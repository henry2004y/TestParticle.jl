# GPU Boris solver using KernelAbstractions.jl

"""
    adapt_field_to_gpu(field::Field, backend::KA.Backend)

Adapt interpolation fields to GPU memory using Adapt.jl.
Analytic functions are returned unchanged.
"""
function adapt_field_to_gpu(field::Field, backend::Backend)
    backend isa CPU && return field

    # Adapt the inner function (FieldInterpolator or analytic)
    adapted_func = Adapt.adapt(backend, field.field_function)

    return Field{is_time_dependent(field), typeof(adapted_func)}(adapted_func)
end

# Fallback for ZeroField
adapt_field_to_gpu(field::ZeroField, backend::Backend) = field

# The solvers the device can run. A kernel takes one step at a time under a step
# size fixed for the whole ensemble, so only those are candidates; the adaptive
# ones decide their own step on the host, once per step for every trajectory,
# which is the SciML loop's job.
const GPUBorisAlgorithm = Union{Boris, MultistepBoris}


# What reaches the device is a parameter container, not the fields on their own,
# so that the kernel can call the same step functions the CPU solvers call and
# reach the fields the same way, through get_q2m, get_EField and get_BField. The
# layout is the one those accessors assume, `(q2m, m, E, B, ...)`, and it is built
# once on the host with the adapted fields.

@inline function boris_update_xv!(i, xv_in, xv_out, p, dt, t, alg)
    r = SVector(xv_in[1, i], xv_in[2, i], xv_in[3, i])
    v_half = SVector(xv_in[4, i], xv_in[5, i], xv_in[6, i])

    r_new, v_half_new = boris_advance(v_half, r, dt, t, p, alg)

    # Scalar write for GPU compatibility
    xv_out[1, i] = r_new[1]
    xv_out[2, i] = r_new[2]
    xv_out[3, i] = r_new[3]
    xv_out[4, i] = v_half_new[1]
    xv_out[5, i] = v_half_new[2]
    xv_out[6, i] = v_half_new[3]

    return
end

@inline function boris_retard_v!(i, xv_in, xv_out, p, dt, t, alg)
    r = SVector(xv_in[1, i], xv_in[2, i], xv_in[3, i])
    v = SVector(xv_in[4, i], xv_in[5, i], xv_in[6, i])

    v_half = boris_half_velocity(v, r, dt, t, p, alg)

    xv_out[4, i] = v_half[1]
    xv_out[5, i] = v_half[2]
    xv_out[6, i] = v_half[3]

    return
end

@kernel function boris_velocity_kernel!(
        xv_out, @Const(xv_in), p, @Const(dt), @Const(t), @Const(offset), @Const(alg)
    )
    i = @index(Global) + offset
    boris_retard_v!(i, xv_in, xv_out, p, dt, t, alg)
end

@kernel function boris_update_kernel!(
        @Const(xv_in), xv_out, p, @Const(dt), @Const(t), @Const(offset), @Const(alg)
    )
    i = @index(Global) + offset
    boris_update_xv!(i, xv_in, xv_out, p, dt, t, alg)
end

@inline function boris_step!(
        backend::Backend, xv_in, xv_out, p, dt, t, irange, workgroup_size, alg
    )
    offset = irange.start - 1
    n_particles = length(irange)
    kernel! = boris_update_kernel!(backend, workgroup_size)
    kernel!(xv_in, xv_out, p, dt, t, offset, alg; ndrange = n_particles)
    synchronize(backend)
    return
end

@inline function boris_step!(::CPU, xv_in, xv_out, p, dt, t, irange, workgroup_size, alg)
    @inbounds for i in irange
        boris_update_xv!(i, xv_in, xv_out, p, dt, t, alg)
    end
    return
end

@inline function boris_velocity_step!(
        backend::Backend, xv_in, xv_out, p, dt, t, irange, workgroup_size, alg
    )
    offset = irange.start - 1
    n_particles = length(irange)
    kernel! = boris_velocity_kernel!(backend, workgroup_size)
    kernel!(xv_out, xv_in, p, dt, t, offset, alg; ndrange = n_particles)
    synchronize(backend)
    return
end

@inline function boris_velocity_step!(
        ::CPU, xv_in, xv_out, p, dt, t, irange, workgroup_size, alg
    )
    @inbounds for i in irange
        boris_retard_v!(i, xv_in, xv_out, p, dt, t, alg)
    end
    return
end

"""
    boris_output_state(xv, p, dt, t, alg)

The state worth reporting for column `xv` of a device array: the position as it
stands together with the velocity synchronised to it. Synchronising costs a
field evaluation, which is why the driver asks for it only at the times it saves
rather than every step.
"""
@inline function boris_output_state(xv, p, dt, t, alg)
    T = eltype(xv)
    r = SVector{3, T}(xv[1], xv[2], xv[3])
    v_half = SVector{3, T}(xv[4], xv[5], xv[6])

    return vcat(r, boris_node_velocity(v_half, r, dt, t, p, alg))
end


@inbounds function _solve_serial(
        prob::TraceProblem, backend::Backend, irange;
        dt::AbstractFloat, plan, save_start::Bool,
        save_end::Bool, save_everystep::Bool, workgroup_size::Int,
        xv_current, xv_next, xv_cpu_buffer, is_cpu_accessible,
        p_gpu, p_host, alg, nout, nt
    )
    (; tspan) = prob
    T = eltype(xv_current)
    n_particles = length(irange)

    sols = Vector{
        typeof(build_solution(prob, :boris, [tspan[1]], [SVector{6, T}(prob.u0)])),
    }(undef, n_particles)

    saved_data = [Vector{SVector{6, T}}(undef, nout) for _ in 1:n_particles]
    saved_times = [Vector{typeof(tspan[1] + dt)}(undef, nout) for _ in 1:n_particles]
    iout_counters = zeros(Int, n_particles)

    nsave = length(plan.times)
    isave = 1
    # Interpolating inside a step needs both ends of it. The device keeps the
    # state at the start of the step in `xv_next` after each swap, so the two are
    # staged through host buffers of their own rather than relying on
    # `xv_cpu_buffer`, which tracks only the array it was aliased to. Like
    # `xv_cpu_buffer` they are sized by the whole ensemble and indexed by the
    # global particle number, because a thread only owns a slice of it.
    if use_saveat(plan)
        xv_cpu_prev = zeros(T, size(xv_current))
        xv_cpu_cur = similar(xv_cpu_prev)
    else
        xv_cpu_prev = Matrix{T}(undef, 0, 0)
        xv_cpu_cur = Matrix{T}(undef, 0, 0)
    end

    if save_start
        if !is_cpu_accessible
            copyto!(xv_cpu_buffer, xv_current)
        end
        for (local_i, i) in enumerate(irange)
            iout_counters[local_i] += 1
            saved_data[local_i][iout_counters[local_i]] =
                SVector{6, T}(
                xv_cpu_buffer[1, i], xv_cpu_buffer[2, i], xv_cpu_buffer[3, i],
                xv_cpu_buffer[4, i], xv_cpu_buffer[5, i], xv_cpu_buffer[6, i]
            )
            saved_times[local_i][iout_counters[local_i]] = tspan[1]
        end
    end

    boris_velocity_step!(
        backend, xv_current, xv_current, p_gpu, dt, tspan[1], irange, workgroup_size, alg
    )

    for it in 1:nt
        t = tspan[1] + (it - 0.5) * dt

        boris_step!(
            backend, xv_current, xv_next, p_gpu, dt, t, irange, workgroup_size, alg
        )

        xv_current, xv_next = xv_next, xv_current

        if use_saveat(plan)
            t_current = tspan[1] + it * dt
            if isave <= nsave && _saveat_reached(plan.times[isave], t_current, plan.dir)
                # The device copies are made once per event rather than per step,
                # so their cost stays proportional to the number of samples.
                copyto!(xv_cpu_prev, xv_next)
                copyto!(xv_cpu_cur, xv_current)

                t_prev = t_current - dt

                while isave <= nsave &&
                        _saveat_reached(plan.times[isave], t_current, plan.dir)
                    t_target = plan.times[isave]
                    for (local_i, i) in enumerate(irange)
                        if iout_counters[local_i] < nout
                            iout_counters[local_i] += 1
                            y_prev = boris_output_state(
                                @view(xv_cpu_prev[:, i]), p_host, dt, t_prev, alg
                            )
                            y_cur = boris_output_state(
                                @view(xv_cpu_cur[:, i]), p_host, dt, t_current, alg
                            )
                            saved_data[local_i][iout_counters[local_i]] =
                                _saveat_interpolate(
                                t_prev, y_prev, t_current, y_cur, t_target
                            )
                            saved_times[local_i][iout_counters[local_i]] = t_target
                        end
                    end
                    isave += 1
                end
            end
        elseif save_everystep
            if !is_cpu_accessible
                copyto!(xv_cpu_buffer, xv_current)
            end

            t_current = tspan[1] + it * dt

            for (local_i, i) in enumerate(irange)
                if iout_counters[local_i] < nout
                    iout_counters[local_i] += 1
                    saved_data[local_i][iout_counters[local_i]] = boris_output_state(
                        @view(xv_cpu_buffer[:, i]), p_host, dt, t_current, alg
                    )
                    saved_times[local_i][iout_counters[local_i]] = t_current
                end
            end
        end
    end

    if save_end
        if !is_cpu_accessible
            copyto!(xv_cpu_buffer, xv_current)
        end
        t_current = tspan[2]

        for (local_i, i) in enumerate(irange)
            if iout_counters[local_i] < nout
                iout_counters[local_i] += 1
                saved_data[local_i][iout_counters[local_i]] = boris_output_state(
                    @view(xv_cpu_buffer[:, i]), p_host, dt, t_current, alg
                )
                saved_times[local_i][iout_counters[local_i]] = t_current
            end
        end
    end

    for local_i in 1:n_particles
        actual_len = iout_counters[local_i]
        if actual_len < nout
            resize!(saved_data[local_i], actual_len)
            resize!(saved_times[local_i], actual_len)
            retcode = ReturnCode.Terminated
        else
            retcode = ReturnCode.Success
        end

        interp = LinearInterpolation(saved_times[local_i], saved_data[local_i])
        sols[local_i] = build_solution(
            prob, :boris, saved_times[local_i], saved_data[local_i];
            interp, retcode, stats = nothing
        )
    end

    return sols
end

function _prepare_boris_solve(
        prob::TraceProblem, backend::Backend, trajectories::Int, dt::AbstractFloat,
        plan, save_start::Bool, save_end::Bool, save_everystep::Bool, maxiters::Int
    )
    (; tspan, p, u0) = prob
    q2m, m, Efunc, Bfunc, _ = p
    T = eltype(u0)

    if abs(dt) < 10 * eps(typeof(dt))
        throw(ArgumentError("time step dt is too small, violating min_dt = 10 * eps(typeof(dt))"))
    end

    Efunc_gpu = adapt_field_to_gpu(Efunc, backend)
    Bfunc_gpu = adapt_field_to_gpu(Bfunc, backend)

    ttotal = tspan[2] - tspan[1]
    nt = round(Int, abs(ttotal / dt))

    if nt > maxiters
        throw(ArgumentError("number of iterations nt ($nt) exceeds maxiters ($maxiters)"))
    end

    nout = save_start ? 1 : 0

    if use_saveat(plan)
        # One slot per requested time, plus the end of the run.
        nout += length(plan.times) + (save_end ? 1 : 0)
    elseif save_everystep
        last_is_step = nt > 0
        nout += nt
        if !save_end && last_is_step
            nout -= 1
        end
        if save_end && !last_is_step
            nout += 1
        end
    elseif save_end
        nout += 1
    end

    n_particles = trajectories
    xv_current = KA.zeros(backend, T, 6, n_particles)
    xv_next = KA.zeros(backend, T, 6, n_particles)

    is_cpu_accessible = xv_current isa Array

    if is_cpu_accessible
        xv_init = xv_current
    else
        xv_init = zeros(T, 6, n_particles)
    end

    for i in 1:n_particles
        new_prob = prob.prob_func(prob, (sim_id = i, repeat = false))
        u0_i = new_prob.u0
        xv_init[:, i] .= u0_i
    end

    if !is_cpu_accessible
        copyto!(xv_current, xv_init)
    end

    if is_cpu_accessible
        xv_cpu_buffer = xv_current
    else
        xv_cpu_buffer = zeros(T, 6, n_particles)
    end

    # One container per side: the device reads the adapted fields, the host the
    # original ones, since a state is synchronised for output after it is copied
    # back. Both answer get_q2m, get_EField and get_BField.
    p_gpu = (q2m, m, Efunc_gpu, Bfunc_gpu)
    p_host = (q2m, m, Efunc, Bfunc)

    return (;
        nt, nout, xv_current, xv_next, xv_cpu_buffer, is_cpu_accessible,
        p_gpu, p_host, tspan, u0, T,
    )
end

@inbounds function solve(
        prob::TraceProblem, alg::GPUBorisAlgorithm, backend::Backend, ::EnsembleSerial;
        dt::AbstractFloat, trajectories::Int = 1,
        saveat = (),
        save_start::Bool = true, save_end::Bool = true, save_everystep::Bool = true,
        workgroup_size::Int = 256, maxiters::Int = 1_000_000
    )
    plan = SavingPlan(
        saveat, prob.tspan, _span_direction(prob.tspan),
        typeof(prob.tspan[1] + dt)
    )
    (;
        nt, nout, xv_current, xv_next, xv_cpu_buffer, is_cpu_accessible,
        p_gpu, p_host,
    ) = _prepare_boris_solve(
        prob, backend, trajectories, dt, plan,
        save_start, save_end, save_everystep, maxiters
    )

    elapsed_time = @elapsed sols = _solve_serial(
        prob, backend, 1:trajectories;
        dt, plan, save_start, save_end, save_everystep, workgroup_size,
        xv_current, xv_next, xv_cpu_buffer, is_cpu_accessible,
        p_gpu, p_host, alg, nout, nt
    )

    return EnsembleSolution(sols, elapsed_time, true)
end

@inbounds function solve(
        prob::TraceProblem, alg::GPUBorisAlgorithm, backend::Backend, ::EnsembleThreads;
        dt::AbstractFloat, trajectories::Int = 1,
        saveat = (),
        save_start::Bool = true, save_end::Bool = true, save_everystep::Bool = true,
        workgroup_size::Int = 256, maxiters::Int = 1_000_000
    )
    plan = SavingPlan(
        saveat, prob.tspan, _span_direction(prob.tspan),
        typeof(prob.tspan[1] + dt)
    )
    (;
        nt, nout, xv_current, xv_next, xv_cpu_buffer, is_cpu_accessible,
        p_gpu, p_host, tspan, u0, T,
    ) = _prepare_boris_solve(
        prob, backend, trajectories, dt, plan,
        save_start, save_end, save_everystep, maxiters
    )

    sols = Vector{
        typeof(build_solution(prob, :boris, [tspan[1]], [SVector{6, T}(u0)])),
    }(undef, trajectories)

    nchunks = Threads.nthreads()
    elapsed_time = @elapsed Threads.@threads for irange in index_chunks(1:trajectories; n = nchunks)
        chunk_sols = _solve_serial(
            prob, backend, irange;
            dt, plan, save_start, save_end, save_everystep, workgroup_size,
            xv_current, xv_next, xv_cpu_buffer, is_cpu_accessible,
            p_gpu, p_host, alg, nout, nt
        )
        for (local_i, i) in enumerate(irange)
            sols[i] = chunk_sols[local_i]
        end
    end

    return EnsembleSolution(sols, elapsed_time, true)
end

@inbounds function solve(
        prob::TraceProblem, alg::GPUBorisAlgorithm, backend::Backend,
        ensemblealg::BasicEnsembleAlgorithm = EnsembleSerial();
        dt::AbstractFloat, trajectories::Int = 1,
        saveat = (),
        save_start::Bool = true, save_end::Bool = true, save_everystep::Bool = true,
        workgroup_size::Int = 256, maxiters::Int = 1_000_000
    )
    return solve(
        prob, alg, backend, ensemblealg;
        dt, trajectories, saveat, save_start, save_end,
        save_everystep, workgroup_size, maxiters
    )
end

function solve(
        prob::TraceProblem, alg::AbstractBoris, backend::Backend, args...; kwargs...
    )
    supported = "Boris() and MultistepBoris{N}(; n)"
    throw(
        ArgumentError(
            "$alg has no GPU path, because it chooses its own time step. Solve it on " *
            "the CPU, `solve(prob, alg)`, or pick a fixed step solver, one of $supported."
        )
    )
end
