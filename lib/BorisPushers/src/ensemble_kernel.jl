# EnsembleKernel solver for hardware backends via KernelAbstractions.jl

"""
    EnsembleKernel(backend=CPU(); workgroup_size=256)

Ensemble algorithm for running particle traces on hardware backends
via KernelAbstractions.jl (such as CUDA, AMDGPU, oneAPI, Metal, or multi-threaded CPU).
"""
struct EnsembleKernel{B <: Backend} <: SciMLBase.EnsembleAlgorithm
    backend::B
    workgroup_size::Int
end

function EnsembleKernel(backend::Backend = CPU(); workgroup_size::Int = 256)
    return EnsembleKernel(backend, workgroup_size)
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

function _foreach_chunked(f, n::Int)
    if n <= 16 || Threads.nthreads() == 1
        for i in 1:n
            f(i)
        end
    else
        chunks = index_chunks(1:n; n = Threads.nthreads())
        Threads.@threads for chunk in chunks
            for i in chunk
                f(i)
            end
        end
    end
    return nothing
end

function _unpack_endpoint_solutions!(
        sols, irange, xv_init, xv_cpu_end, saved_times, nout,
        save_start::Bool, save_end::Bool, prob, alg, ::Type{T}
    ) where {T}
    n_particles = length(irange)
    _foreach_chunked(n_particles) do local_i
        i = irange[local_i]
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
    return sols
end

function _unpack_dense_solutions!(
        sols, irange, saved_data_buf, saved_times, nout, prob, alg, ::Type{T}
    ) where {T}
    n_particles = length(irange)
    _foreach_chunked(n_particles) do local_i
        i = irange[local_i]
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
    return sols
end

@inline function _call_prob_func(prob_func, prob, i::Int, repeat::Int, seed)
    ctx = (sim_id = i, repeat = repeat)
    return prob_func(prob, ctx)
end


function _init_particles!(
        xv_init::AbstractMatrix{T}, prob, prob_func, n_particles::Int, seed
    ) where {T}
    if n_particles == 1
        new_prob = _call_prob_func(prob_func, prob, 1, 0, seed)
        u0 = new_prob.u0
        @inbounds for c in 1:6
            xv_init[1, c] = u0[c]
        end
    else
        _foreach_chunked(n_particles) do i
            new_prob = _call_prob_func(prob_func, prob, i, 0, seed)
            u0_i = new_prob.u0
            @inbounds for c in 1:6
                xv_init[i, c] = u0_i[c]
            end
        end
    end
    return xv_init
end

function _solve_cpu_chunk!(
        sols, irange, xv_current, plan, dt, tspan,
        p_host, alg, nout, nt, save_start, save_end, save_everystep,
        saved_times, prob, ::Type{T}
    ) where {T}
    nsave = length(plan.times)
    for (local_i, i) in enumerate(irange)
        traj = Vector{SVector{6, T}}(undef, nout)
        idx = 0
        r = SVector{3, T}(xv_current[i, 1], xv_current[i, 2], xv_current[i, 3])
        v = SVector{3, T}(xv_current[i, 4], xv_current[i, 5], xv_current[i, 6])

        if save_start
            idx += 1
            traj[idx] = vcat(r, v)
        end

        v_half = update_velocity_half(v, r, dt, tspan[1], p_host, alg)
        isave = 1
        r_prev = r
        v_half_prev = v_half
        t_prev = tspan[1]

        for it in 1:nt
            t_current = tspan[1] + it * dt
            r, v_half = advance_boris(v_half_prev, r_prev, dt, t_prev, p_host, alg)

            if use_saveat(plan)
                if isave <= nsave && saveat_reached(plan.times[isave], t_current, plan.dir)
                    v_prev = update_velocity_node(
                        v_half_prev, r_prev, dt, t_prev, p_host, alg
                    )
                    v_cur = update_velocity_node(v_half, r, dt, t_current, p_host, alg)
                    y_prev = vcat(r_prev, v_prev)
                    y_cur = vcat(r, v_cur)

                    while isave <= nsave &&
                            saveat_reached(plan.times[isave], t_current, plan.dir)
                        t_target = plan.times[isave]
                        idx += 1
                        traj[idx] = saveat_interpolate(
                            t_prev, y_prev, t_current, y_cur, t_target
                        )
                        isave += 1
                    end
                end
            elseif save_everystep
                if it < nt || save_end
                    idx += 1
                    v_node = update_velocity_node(
                        v_half, r, dt, t_current, p_host, alg
                    )
                    traj[idx] = vcat(r, v_node)
                end
            end

            r_prev = r
            v_half_prev = v_half
            t_prev = t_current
        end

        if save_end && (!use_saveat(plan) && !save_everystep)
            idx += 1
            v_node = update_velocity_node(v_half, r, dt, tspan[2], p_host, alg)
            traj[idx] = vcat(r, v_node)
        elseif save_end && use_saveat(plan) && idx < nout
            idx += 1
            v_node = update_velocity_node(v_half, r, dt, tspan[2], p_host, alg)
            traj[idx] = vcat(r, v_node)
        end

        interp = LinearInterpolation(saved_times, traj)
        sols[local_i] = build_solution(
            prob, alg, saved_times, traj;
            interp, retcode = ReturnCode.Success, stats = nothing
        )
    end
    return sols
end

function _solve_kernel_serial(
        prob, backend::Backend, irange;
        dt, plan, save_start::Bool,
        save_end::Bool, save_everystep::Bool, workgroup_size::Int,
        xv_current, xv_init, is_cpu_accessible::Bool,
        p_gpu, p_host, alg, nout::Int, nt::Int
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

    if backend isa CPU
        return _solve_cpu_chunk!(
            sols, irange, xv_current, plan, dt, tspan,
            p_host, alg, nout, nt, save_start, save_end, save_everystep,
            saved_times, prob, T
        )
    end

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

        _unpack_endpoint_solutions!(
            sols, irange, xv_init, xv_cpu_end, saved_times, nout, save_start,
            save_end, prob, alg, T
        )
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

        _unpack_dense_solutions!(
            sols, irange, saved_data_buf, saved_times, nout, prob, alg, T
        )
    end

    return sols
end

function _prepare_ensemble_solve(
        prob, prob_func, backend::Backend, trajectories::Int, dt::Real,
        plan, save_start::Bool, save_end::Bool, save_everystep::Bool, maxiters::Int;
        sort_particles::Bool = false, seed = nothing
    )
    (; tspan, p, u0) = prob
    T = eltype(u0)
    dt_T = T(dt)
    tspan_T = (T(tspan[1]), T(tspan[2]))

    timescale = max(abs(tspan_T[1]), abs(tspan_T[2]), abs(tspan_T[2] - tspan_T[1]))
    min_dt = 10 * eps(T) * timescale
    if abs(dt_T) < min_dt
        throw(
            ArgumentError(
                "time step dt is too small, violating min_dt = 10 * eps(typeof(dt))"
            )
        )
    end

    p_gpu = adapt_params(p, backend)
    p_host = p

    ttotal = tspan_T[2] - tspan_T[1]
    nt = round(Int, abs(ttotal / dt_T))

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

    xv_init = if is_cpu_accessible && (!sort_particles || n_particles <= 1)
        xv_current
    else
        zeros(T, n_particles, 6)
    end

    _init_particles!(xv_init, prob, prob_func, n_particles, seed)

    if sort_particles && n_particles > 1
        perm = morton_sort_particles(xv_init)
        xv_init_sorted = xv_init[perm, :]
        copyto!(xv_current, xv_init_sorted)
        return (;
            nt, nout, xv_current, xv_init = xv_init_sorted, is_cpu_accessible,
            p_gpu, p_host, tspan = tspan_T, u0, T, dt = dt_T, perm
        )
    else
        if xv_current !== xv_init
            copyto!(xv_current, xv_init)
        end
        return (;
            nt, nout, xv_current, xv_init, is_cpu_accessible,
            p_gpu, p_host, tspan = tspan_T, u0, T, dt = dt_T, perm = nothing
        )
    end
end

function _execute_ensemble_kernel(
        base_prob, prob_func, alg::AbstractBoris, ens::EnsembleKernel;
        dt::Real, trajectories::Int = 1,
        saveat = (),
        save_start::Bool = true, save_end::Bool = true, save_everystep::Bool = true,
        maxiters::Int = 1_000_000,
        sort_particles::Bool = false,
        seed = nothing
    )
    if SciMLBase.isadaptive(alg)
        throw(
            ArgumentError(
                "$alg has no GPU path, because it chooses its own time step. " *
                    "Solve it on the CPU, `solve(prob, alg)`, or pick one of " *
                    "Boris() and MultistepBoris{N}(; n)."
            )
        )
    end

    backend = ens.backend
    workgroup_size = ens.workgroup_size
    T = eltype(base_prob.u0)
    dt_T = T(dt)
    tspan_T = (T(base_prob.tspan[1]), T(base_prob.tspan[2]))
    time_type = typeof(tspan_T[1] + dt_T)
    plan = SavingPlan(saveat, tspan_T, time_type)

    (;
        nt, nout, xv_current, xv_init, is_cpu_accessible,
        p_gpu, p_host, tspan, u0, perm
    ) = _prepare_ensemble_solve(
        base_prob, prob_func, backend, trajectories, dt_T, plan,
        save_start, save_end, save_everystep, maxiters;
        sort_particles, seed
    )

    t0 = time_ns()

    sols = if backend isa CPU && Threads.nthreads() > 1 && trajectories > 16
        saved_times = _build_saved_times(
            tspan, dt_T, nt, plan, save_start, save_end, save_everystep, time_type
        )
        res_sols = Vector{
            typeof(build_solution(base_prob, alg, saved_times, [SVector{6, T}(u0)])),
        }(undef, trajectories)

        chunks = index_chunks(1:trajectories; n = Threads.nthreads())
        Threads.@threads for irange in chunks
            chunk_sols = _solve_kernel_serial(
                base_prob, backend, irange;
                dt = dt_T, plan, save_start, save_end, save_everystep, workgroup_size,
                xv_current, xv_init, is_cpu_accessible,
                p_gpu, p_host, alg, nout, nt
            )
            for (local_i, i) in enumerate(irange)
                res_sols[i] = chunk_sols[local_i]
            end
        end
        res_sols
    else
        _solve_kernel_serial(
            base_prob, backend, 1:trajectories;
            dt = dt_T, plan, save_start, save_end, save_everystep, workgroup_size,
            xv_current, xv_init, is_cpu_accessible,
            p_gpu, p_host, alg, nout, nt
        )
    end

    elapsed_time = (time_ns() - t0) * 1.0e-9

    sols_final = perm === nothing ? sols : sols[invperm(perm)]
    return sols_final, elapsed_time
end

function SciMLBase.__solve(
        eprob::SciMLBase.AbstractEnsembleProblem, alg::AbstractBoris, ens::EnsembleKernel;
        dt::Real, trajectories::Int = 1,
        saveat = (),
        save_start::Bool = true, save_end::Bool = true, save_everystep::Bool = true,
        maxiters::Int = 1_000_000,
        sort_particles::Bool = false,
        seed = nothing,
        kwargs...
    )
    base_prob = eprob.prob
    prob_func = eprob.prob_func
    output_func = eprob.output_func
    reduction = eprob.reduction

    sols, elapsed_time = _execute_ensemble_kernel(
        base_prob, prob_func, alg, ens;
        dt, trajectories, saveat, save_start, save_end,
        save_everystep, maxiters, sort_particles, seed
    )

    final_sols = if output_func !== SciMLBase.DEFAULT_OUTPUT_FUNC
        first_out, rerun = output_func(sols[1], 1)
        if rerun
            throw(
                ArgumentError("rerun from output_func is not supported in EnsembleKernel")
            )
        end
        transformed = Vector{typeof(first_out)}(undef, length(sols))
        transformed[1] = first_out
        for i in 2:length(sols)
            out, rerun = output_func(sols[i], i)
            if rerun
                throw(
                    ArgumentError(
                        "rerun from output_func is not supported in EnsembleKernel"
                    )
                )
            end
            transformed[i] = out
        end
        transformed
    else
        sols
    end

    if reduction !== SciMLBase.DEFAULT_REDUCTION
        u, _ = reduction(nothing, final_sols, 1:trajectories)
        return EnsembleSolution(u, elapsed_time, true)
    end

    return EnsembleSolution(final_sols, elapsed_time, true)
end

function SciMLBase.solve(
        eprob::SciMLBase.EnsembleProblem, alg::AbstractBoris, ens::EnsembleKernel;
        kwargs...
    )
    return SciMLBase.__solve(eprob, alg, ens; kwargs...)
end

function SciMLBase.solve(
        eprob::SciMLBase.AbstractEnsembleProblem, alg::AbstractBoris, ens::EnsembleKernel;
        kwargs...
    )
    return SciMLBase.__solve(eprob, alg, ens; kwargs...)
end

function SciMLBase.solve(
        prob::AbstractODEProblem, alg::AbstractBoris, ens::EnsembleKernel;
        trajectories::Int = 1, kwargs...
    )
    prob_func = hasproperty(prob, :prob_func) ? prob.prob_func : ((p, ctx) -> p)
    eprob = EnsembleProblem(prob; prob_func)
    return SciMLBase.__solve(eprob, alg, ens; trajectories, kwargs...)
end

# Backward compatibility overloads
function SciMLBase.solve(
        prob::AbstractODEProblem,
        alg::AbstractBoris, backend::Backend;
        workgroup_size::Int = 256, kwargs...
    )
    return SciMLBase.solve(prob, alg, EnsembleKernel(backend; workgroup_size); kwargs...)
end

function SciMLBase.solve(
        prob::SciMLBase.EnsembleProblem,
        alg::AbstractBoris, backend::Backend;
        workgroup_size::Int = 256, kwargs...
    )
    return SciMLBase.solve(prob, alg, EnsembleKernel(backend; workgroup_size); kwargs...)
end

function SciMLBase.solve(
        prob::AbstractODEProblem,
        alg::AbstractBoris, backend::Backend, ::BasicEnsembleAlgorithm;
        workgroup_size::Int = 256, kwargs...
    )
    return SciMLBase.solve(prob, alg, EnsembleKernel(backend; workgroup_size); kwargs...)
end

function SciMLBase.solve(
        prob::SciMLBase.EnsembleProblem,
        alg::AbstractBoris, backend::Backend, ::BasicEnsembleAlgorithm;
        workgroup_size::Int = 256, kwargs...
    )
    return SciMLBase.solve(prob, alg, EnsembleKernel(backend; workgroup_size); kwargs...)
end
