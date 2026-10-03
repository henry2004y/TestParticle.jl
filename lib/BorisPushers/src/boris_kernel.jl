# GPU Boris solver using KernelAbstractions.jl

"""
    adapt_field_to_gpu(field, backend::Backend)

Adapt interpolation fields or parameters to GPU memory using Adapt.jl.
Analytic functions and CPU backend return the input unchanged.
"""
adapt_field_to_gpu(field, ::CPU) = field
adapt_field_to_gpu(field, backend::Backend) = Adapt.adapt(backend, field)

const GPUBorisAlgorithm = Union{Boris, MultistepBoris}

@inline function boris_update_xv!(i, xv_in, xv_out, p, dt, t, alg)
    r = SVector(xv_in[i, 1], xv_in[i, 2], xv_in[i, 3])
    v_half = SVector(xv_in[i, 4], xv_in[i, 5], xv_in[i, 6])

    r_new, v_half_new = advance_boris(v_half, r, dt, t, p, alg)

    xv_out[i, 1] = r_new[1]
    xv_out[i, 2] = r_new[2]
    xv_out[i, 3] = r_new[3]
    xv_out[i, 4] = v_half_new[1]
    xv_out[i, 5] = v_half_new[2]
    xv_out[i, 6] = v_half_new[3]

    return
end

@inline function boris_retard_v!(i, xv_in, xv_out, p, dt, t, alg)
    r = SVector(xv_in[i, 1], xv_in[i, 2], xv_in[i, 3])
    v = SVector(xv_in[i, 4], xv_in[i, 5], xv_in[i, 6])

    v_half = update_velocity_half(v, r, dt, t, p, alg)

    xv_out[i, 4] = v_half[1]
    xv_out[i, 5] = v_half[2]
    xv_out[i, 6] = v_half[3]

    return
end

@kernel function boris_velocity_kernel!(
        xv_out, @Const(xv_in), p, @Const(dt), @Const(t), @Const(offset), @Const(alg)
    )
    idx = @index(Global)
    i = idx + offset
    boris_retard_v!(i, xv_in, xv_out, p, dt, t, alg)
end

@kernel function boris_update_kernel!(
        @Const(xv_in), xv_out, p, @Const(dt), @Const(t), @Const(offset), @Const(alg)
    )
    idx = @index(Global)
    i = idx + offset
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

@kernel function boris_trajectory_kernel!(
        xv, p, @Const(dt), @Const(t0), @Const(nt), @Const(offset), @Const(alg)
    )
    idx = @index(Global)
    i = idx + offset
    r = SVector(xv[i, 1], xv[i, 2], xv[i, 3])
    v = SVector(xv[i, 4], xv[i, 5], xv[i, 6])

    v_half = update_velocity_half(v, r, dt, t0, p, alg)

    for it in 1:nt
        t = t0 + (it - 1) * dt
        r, v_half = advance_boris(v_half, r, dt, t, p, alg)
    end

    t_end = t0 + nt * dt
    v_end = update_velocity_node(v_half, r, dt, t_end, p, alg)

    xv[i, 1] = r[1]
    xv[i, 2] = r[2]
    xv[i, 3] = r[3]
    xv[i, 4] = v_end[1]
    xv[i, 5] = v_end[2]
    xv[i, 6] = v_end[3]
end

@kernel function boris_saveat_kernel!(
        saved_data, @Const(xv), @Const(plan_times), p, @Const(dt),
        @Const(t0), @Const(nt), @Const(nsave), @Const(dir),
        @Const(save_start), @Const(save_end), @Const(offset), @Const(alg)
    )
    idx = @index(Global)
    i = idx + offset
    r = SVector(xv[i, 1], xv[i, 2], xv[i, 3])
    v = SVector(xv[i, 4], xv[i, 5], xv[i, 6])

    iout = 0
    if save_start
        iout += 1
        saved_data[i, 1, iout] = r[1]
        saved_data[i, 2, iout] = r[2]
        saved_data[i, 3, iout] = r[3]
        saved_data[i, 4, iout] = v[1]
        saved_data[i, 5, iout] = v[2]
        saved_data[i, 6, iout] = v[3]
    end

    v_half = update_velocity_half(v, r, dt, t0, p, alg)
    isave = 1

    for it in 1:nt
        t_prev = t0 + (it - 1) * dt
        r_prev = r
        v_half_prev = v_half

        r, v_half = advance_boris(v_half, r, dt, t_prev, p, alg)
        t_current = t0 + it * dt

        if isave <= nsave && _saveat_reached(plan_times[isave], t_current, dir)
            v_prev = update_velocity_node(v_half_prev, r_prev, dt, t_prev, p, alg)
            v_cur = update_velocity_node(v_half, r, dt, t_current, p, alg)

            y_prev = SVector(
                r_prev[1], r_prev[2], r_prev[3],
                v_prev[1], v_prev[2], v_prev[3]
            )
            y_cur = SVector(
                r[1], r[2], r[3],
                v_cur[1], v_cur[2], v_cur[3]
            )

            while isave <= nsave && _saveat_reached(plan_times[isave], t_current, dir)
                t_target = plan_times[isave]
                y_target = _saveat_interpolate(t_prev, y_prev, t_current, y_cur, t_target)
                iout += 1
                saved_data[i, 1, iout] = y_target[1]
                saved_data[i, 2, iout] = y_target[2]
                saved_data[i, 3, iout] = y_target[3]
                saved_data[i, 4, iout] = y_target[4]
                saved_data[i, 5, iout] = y_target[5]
                saved_data[i, 6, iout] = y_target[6]
                isave += 1
            end
        end
    end

    if save_end
        t_end = t0 + nt * dt
        v_end = update_velocity_node(v_half, r, dt, t_end, p, alg)
        iout += 1
        saved_data[i, 1, iout] = r[1]
        saved_data[i, 2, iout] = r[2]
        saved_data[i, 3, iout] = r[3]
        saved_data[i, 4, iout] = v_end[1]
        saved_data[i, 5, iout] = v_end[2]
        saved_data[i, 6, iout] = v_end[3]
    end
end

@kernel function boris_everystep_kernel!(
        saved_data, @Const(xv), p, @Const(dt),
        @Const(t0), @Const(nt), @Const(save_start), @Const(save_end),
        @Const(offset), @Const(alg)
    )
    idx = @index(Global)
    i = idx + offset
    r = SVector(xv[i, 1], xv[i, 2], xv[i, 3])
    v = SVector(xv[i, 4], xv[i, 5], xv[i, 6])

    iout = 0
    if save_start
        iout += 1
        saved_data[i, 1, iout] = r[1]
        saved_data[i, 2, iout] = r[2]
        saved_data[i, 3, iout] = r[3]
        saved_data[i, 4, iout] = v[1]
        saved_data[i, 5, iout] = v[2]
        saved_data[i, 6, iout] = v[3]
    end

    v_half = update_velocity_half(v, r, dt, t0, p, alg)

    for it in 1:nt
        t_prev = t0 + (it - 1) * dt
        r, v_half = advance_boris(v_half, r, dt, t_prev, p, alg)
        t_current = t0 + it * dt

        if it < nt || save_end
            v_node = update_velocity_node(v_half, r, dt, t_current, p, alg)
            iout += 1
            saved_data[i, 1, iout] = r[1]
            saved_data[i, 2, iout] = r[2]
            saved_data[i, 3, iout] = r[3]
            saved_data[i, 4, iout] = v_node[1]
            saved_data[i, 5, iout] = v_node[2]
            saved_data[i, 6, iout] = v_node[3]
        end
    end
end

"""
    boris_output_state(xv, p, dt, t, alg)

The state worth reporting for column `xv` of a device array: the position as it
stands together with the velocity synchronised to it.
"""
@inline function boris_output_state(xv, p, dt, t, alg)
    T = eltype(xv)
    r = SVector{3, T}(xv[1], xv[2], xv[3])
    v_half = SVector{3, T}(xv[4], xv[5], xv[6])

    return vcat(r, update_velocity_node(v_half, r, dt, t, p, alg))
end

function _build_saved_times(
        tspan, dt, nt, plan, save_start::Bool, save_end::Bool,
        save_everystep::Bool, ::Type{time_type}
    ) where {time_type}
    saved_times = time_type[]
    if save_start
        push!(saved_times, time_type(tspan[1]))
    end
    if use_saveat(plan)
        append!(saved_times, plan.times)
        if save_end
            push!(saved_times, time_type(tspan[2]))
        end
    elseif save_everystep
        for it in 1:nt
            if it < nt || save_end
                push!(saved_times, time_type(tspan[1] + it * dt))
            end
        end
        if nt == 0 && save_end
            push!(saved_times, time_type(tspan[2]))
        end
    elseif save_end
        push!(saved_times, time_type(tspan[2]))
    end
    return saved_times
end

@inbounds function _solve_serial(
        prob::AbstractODEProblem, backend::Backend, irange;
        dt::AbstractFloat, plan, save_start::Bool,
        save_end::Bool, save_everystep::Bool, workgroup_size::Int,
        xv_current, xv_init, is_cpu_accessible,
        p_gpu, p_host, alg, nout, nt
    )
    (; tspan) = prob
    T = eltype(xv_current)
    n_particles = length(irange)
    offset = irange.start - 1
    time_type = typeof(tspan[1] + dt)
    saved_times = _build_saved_times(
        tspan, dt, nt, plan, save_start, save_end, save_everystep, time_type
    )

    sols = Vector{
        typeof(build_solution(prob, alg, saved_times, [SVector{6, T}(prob.u0)])),
    }(undef, n_particles)

    if !use_saveat(plan) && !save_everystep
        traj_kernel! = boris_trajectory_kernel!(backend, workgroup_size)
        traj_kernel!(
            xv_current, p_gpu, dt, tspan[1], nt, offset, alg;
            ndrange = n_particles
        )
        synchronize(backend)

        xv_cpu_end = is_cpu_accessible ? xv_current : zeros(T, size(xv_current, 1), 6)
        if !is_cpu_accessible
            copyto!(xv_cpu_end, xv_current)
        end

        for (local_i, i) in enumerate(irange)
            traj = Vector{SVector{6, T}}(undef, nout)
            idx = 0
            if save_start
                idx += 1
                traj[idx] = SVector{6, T}(
                    xv_init[i, 1], xv_init[i, 2], xv_init[i, 3],
                    xv_init[i, 4], xv_init[i, 5], xv_init[i, 6]
                )
            end
            if save_end
                idx += 1
                traj[idx] = SVector{6, T}(
                    xv_cpu_end[i, 1], xv_cpu_end[i, 2], xv_cpu_end[i, 3],
                    xv_cpu_end[i, 4], xv_cpu_end[i, 5], xv_cpu_end[i, 6]
                )
            end
            interp = LinearInterpolation(saved_times, traj)
            sols[local_i] = build_solution(
                prob, alg, saved_times, traj;
                interp, retcode = ReturnCode.Success, stats = nothing
            )
        end
    else
        n_total = size(xv_current, 1)
        saved_data_gpu = KA.zeros(backend, T, n_total, 6, nout)

        if use_saveat(plan)
            plan_times_gpu = adapt_field_to_gpu(plan.times, backend)
            nsave = length(plan.times)
            kernel! = boris_saveat_kernel!(backend, workgroup_size)
            kernel!(
                saved_data_gpu, xv_current, plan_times_gpu, p_gpu, dt,
                tspan[1], nt, nsave, plan.dir, save_start, save_end,
                offset, alg; ndrange = n_particles
            )
        else
            kernel! = boris_everystep_kernel!(backend, workgroup_size)
            kernel!(
                saved_data_gpu, xv_current, p_gpu, dt,
                tspan[1], nt, save_start, save_end,
                offset, alg; ndrange = n_particles
            )
        end
        synchronize(backend)

        saved_data_buf = if is_cpu_accessible
            saved_data_gpu
        else
            saved_cpu = zeros(T, n_total, 6, nout)
            copyto!(saved_cpu, saved_data_gpu)
            saved_cpu
        end

        for (local_i, i) in enumerate(irange)
            traj = Vector{SVector{6, T}}(undef, nout)
            for j in 1:nout
                traj[j] = SVector{6, T}(
                    saved_data_buf[i, 1, j],
                    saved_data_buf[i, 2, j],
                    saved_data_buf[i, 3, j],
                    saved_data_buf[i, 4, j],
                    saved_data_buf[i, 5, j],
                    saved_data_buf[i, 6, j]
                )
            end
            interp = LinearInterpolation(saved_times, traj)
            sols[local_i] = build_solution(
                prob, alg, saved_times, traj;
                interp, retcode = ReturnCode.Success, stats = nothing
            )
        end
    end

    return sols
end

function _prepare_boris_solve(
        prob::AbstractODEProblem, backend::Backend, trajectories::Int, dt::AbstractFloat,
        plan, save_start::Bool, save_end::Bool, save_everystep::Bool, maxiters::Int
    )
    (; tspan, p, u0) = prob
    q2m, m = get_q2m(p), p[2]
    Efunc, Bfunc = get_EField(p), get_BField(p)
    T = eltype(u0)

    if abs(dt) < 10 * eps(typeof(dt))
        throw(
            ArgumentError(
                "time step dt is too small, violating min_dt = 10 * eps(typeof(dt))"
            )
        )
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
    xv_current = KA.zeros(backend, T, n_particles, 6)
    is_cpu_accessible = xv_current isa Array

    xv_init = zeros(T, n_particles, 6)
    prob_func = hasproperty(prob, :prob_func) ? prob.prob_func : ((p, ctx) -> p)

    for i in 1:n_particles
        u0_i = if n_particles == 1
            prob.u0
        else
            new_prob = prob_func(prob, (sim_id = i, repeat = false))
            new_prob.u0
        end
        for c in 1:6
            xv_init[i, c] = u0_i[c]
        end
    end

    copyto!(xv_current, xv_init)

    p_gpu = (q2m, m, Efunc_gpu, Bfunc_gpu)
    p_host = (q2m, m, Efunc, Bfunc)

    return (;
        nt, nout, xv_current, xv_init, is_cpu_accessible,
        p_gpu, p_host, tspan, u0, T,
    )
end

@inbounds function SciMLBase.solve(
        prob::AbstractODEProblem, alg::GPUBorisAlgorithm, backend::Backend,
        ::EnsembleSerial;
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
        nt, nout, xv_current, xv_init, is_cpu_accessible,
        p_gpu, p_host,
    ) = _prepare_boris_solve(
        prob, backend, trajectories, dt, plan,
        save_start, save_end, save_everystep, maxiters
    )

    elapsed_time = @elapsed sols = _solve_serial(
        prob, backend, 1:trajectories;
        dt, plan, save_start, save_end, save_everystep, workgroup_size,
        xv_current, xv_init, is_cpu_accessible,
        p_gpu, p_host, alg, nout, nt
    )

    return EnsembleSolution(sols, elapsed_time, true)
end

@inbounds function SciMLBase.solve(
        prob::AbstractODEProblem, alg::GPUBorisAlgorithm, backend::Backend,
        ::EnsembleThreads;
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
        nt, nout, xv_current, xv_init, is_cpu_accessible,
        p_gpu, p_host, tspan, u0, T,
    ) = _prepare_boris_solve(
        prob, backend, trajectories, dt, plan,
        save_start, save_end, save_everystep, maxiters
    )

    if !(backend isa CPU)
        elapsed_time = @elapsed sols = _solve_serial(
            prob, backend, 1:trajectories;
            dt, plan, save_start, save_end, save_everystep, workgroup_size,
            xv_current, xv_init, is_cpu_accessible,
            p_gpu, p_host, alg, nout, nt
        )
        return EnsembleSolution(sols, elapsed_time, true)
    end

    time_type = typeof(tspan[1] + dt)
    saved_times = _build_saved_times(
        tspan, dt, nt, plan, save_start, save_end, save_everystep, time_type
    )
    sols = Vector{
        typeof(build_solution(prob, alg, saved_times, [SVector{6, T}(u0)])),
    }(undef, trajectories)

    nchunks = Threads.nthreads()
    chunks = index_chunks(1:trajectories; n = nchunks)
    elapsed_time = @elapsed Threads.@threads for irange in chunks
        chunk_sols = _solve_serial(
            prob, backend, irange;
            dt, plan, save_start, save_end, save_everystep, workgroup_size,
            xv_current, xv_init, is_cpu_accessible,
            p_gpu, p_host, alg, nout, nt
        )
        for (local_i, i) in enumerate(irange)
            sols[i] = chunk_sols[local_i]
        end
    end

    return EnsembleSolution(sols, elapsed_time, true)
end

@inbounds function SciMLBase.solve(
        prob::AbstractODEProblem, alg::GPUBorisAlgorithm, backend::Backend,
        ensemblealg::BasicEnsembleAlgorithm = EnsembleSerial();
        dt::AbstractFloat, trajectories::Int = 1,
        saveat = (),
        save_start::Bool = true, save_end::Bool = true, save_everystep::Bool = true,
        workgroup_size::Int = 256, maxiters::Int = 1_000_000
    )
    return SciMLBase.solve(
        prob, alg, backend, ensemblealg;
        dt, trajectories, saveat, save_start, save_end,
        save_everystep, workgroup_size, maxiters
    )
end

function SciMLBase.solve(
        prob::AbstractODEProblem, alg::AbstractBoris, backend::Backend, args...;
        kwargs...
    )
    supported = "Boris() and MultistepBoris{N}(; n)"
    throw(
        ArgumentError(
            "$alg has no GPU path, because it chooses its own time step. " *
                "Solve it on the CPU, `solve(prob, alg)`, or pick one of $supported."
        )
    )
end
