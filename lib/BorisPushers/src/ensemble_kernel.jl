# EnsembleKernel solver for hardware backends via KernelAbstractions.jl

"""
    EnsembleKernel(backend=CPU(); workgroup_size=256, pin=nothing)

Ensemble algorithm for running particle traces on hardware backends
via KernelAbstractions.jl (such as CUDA, AMDGPU, oneAPI, Metal, or multi-threaded CPU).
When `pin = true`, host staging buffers used for data transfers to and from GPU devices
are pinned (page-locked) to maximize PCIe transfer bandwidth. When `pin === nothing` (default),
pinning is chosen automatically: enabled for GPU runs with `raw_output = true` and
particle count >= 100,000 where transfer bandwidth exceeds pinning latency, and disabled otherwise.
Explicitly passing `pin = true` or `pin = false` overrides the automatic heuristic.
"""
struct EnsembleKernel{B <: Backend} <: SciMLBase.EnsembleAlgorithm
    backend::B
    workgroup_size::Int
    pin::Union{Nothing, Bool}
end

function EnsembleKernel(
        backend::Backend = CPU();
        workgroup_size::Int = 256, pin::Union{Nothing, Bool} = nothing
    )
    return EnsembleKernel(backend, workgroup_size, pin)
end

function EnsembleKernel(backend::B, workgroup_size::Int) where {B <: Backend}
    return EnsembleKernel(backend, workgroup_size, nothing)
end

@inline function _should_pin(
        backend::Backend, pin::Union{Nothing, Bool},
        raw_output::Bool, trajectories::Int
    )
    if pin !== nothing
        return pin
    end
    return !(backend isa CPU) && raw_output && trajectories >= 100_000
end

function _pin_array!(backend::Backend, a::AbstractArray)
    backend isa CPU && return a
    m = parentmodule(typeof(backend))
    while true
        if isdefined(m, :pin)
            try
                m.pin(a)
                return a
            catch
            end
        elseif isdefined(m, :Mem) && isdefined(getfield(m, :Mem), :pin)
            try
                getfield(m, :Mem).pin(a)
                return a
            catch
            end
        end
        parent = parentmodule(m)
        (parent === m || parent === Main || parent === Base || parent === Core) && break
        m = parent
    end
    for pkg in (:CUDA, :AMDGPU, :Metal, :oneAPI)
        if isdefined(Main, pkg)
            mod = getfield(Main, pkg)
            if isdefined(mod, :pin)
                try
                    mod.pin(a)
                    return a
                catch
                end
            elseif isdefined(mod, :Mem) && isdefined(getfield(mod, :Mem), :pin)
                try
                    getfield(mod, :Mem).pin(a)
                    return a
                catch
                end
            end
        end
    end
    return a
end

function _staging_array(backend::Backend, ::Type{T}, dims; pin::Bool = false) where {T}
    a = zeros(T, dims)
    if pin
        _pin_array!(backend, a)
    end
    return a
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

"""
    _term_info(status, nout_out, i, nout) -> (terminated, n_completed, emitted)

Read how the trace of particle `i` ended from the `status` array a kernel
filled: `terminated`, the number of steps that stayed inside the domain, and
how many states were saved for it. A trace that ran the full time span saves
`nout` states and leaves `status` at zero.
"""
@inline function _term_info(
        status::AbstractVector{Int}, nout_out, i::Int, nout::Int
    )
    k = @inbounds status[i]
    k > 0 || return (false, 0, nout)
    emitted = nout_out === nothing ? nout : @inbounds(nout_out[i])
    return (true, k - 1, emitted)
end

"""
    _freeze_trailing!(u, status, nout_out, nout)

Repeat the last saved state of every trace that stopped early into the slots
it did not reach, so raw output never holds zeros for a terminated particle.
"""
function _freeze_trailing!(
        u, status::AbstractVector{Int}, nout_out::AbstractVector{Int}, nout::Int
    )
    @inbounds for i in eachindex(status)
        k = status[i]
        k > 0 || continue
        emitted = nout_out[i]
        emitted == 0 && continue
        emitted >= nout && continue
        for j in (emitted + 1):nout, c in 1:6
            u[i, c, j] = u[i, c, emitted]
        end
    end
    return u
end

_freeze_trailing!(u, ::Nothing, ::Nothing, nout::Int) = u

"""
    _terminated_flags(status, n) -> BitVector

Whether each trace stopped early because it left the domain. Without a boundary
check no trace did.
"""
function _terminated_flags(status, n::Int)
    flags = falses(n)
    status === nothing && return flags
    @inbounds for i in 1:n
        flags[i] = status[i] > 0
    end
    return flags
end

function _unpack_endpoint_solutions!(
        sols, irange, xv_init, xv_cpu_end, saved_times, nout,
        save_start::Bool, save_end::Bool, prob, alg, ::Type{T};
        status = nothing, t0 = nothing, dt = nothing
    ) where {T}
    n_particles = length(irange)
    _foreach_chunked(n_particles) do local_i
        i = irange[local_i]
        terminated, n_completed, _ = status === nothing ?
            (false, 0, 0) : _term_info(status, nothing, i, nout)

        traj = Vector{SVector{6, T}}(undef, nout)
        idx = 0
        # A trace rejected on its first step has no in-domain state beyond the
        # start one, so its end slot is dropped instead of repeating the start.
        drop_end = terminated && n_completed == 0 && save_start
        if save_start
            idx += 1
            traj[idx] = SVector{6, T}(
                xv_init[i, 1], xv_init[i, 2], xv_init[i, 3],
                xv_init[i, 4], xv_init[i, 5], xv_init[i, 6]
            )
        end
        if save_end && !drop_end
            idx += 1
            traj[idx] = SVector{6, T}(
                xv_cpu_end[i, 1], xv_cpu_end[i, 2], xv_cpu_end[i, 3],
                xv_cpu_end[i, 4], xv_cpu_end[i, 5], xv_cpu_end[i, 6]
            )
        end
        t_vec = saved_times
        retcode = ReturnCode.Success
        if terminated
            retcode = ReturnCode.Terminated
            t_vec = saved_times[1:idx]
            traj = traj[1:idx]
            save_end && !drop_end && (t_vec[end] = t0 + n_completed * dt)
        end
        interp = LinearInterpolation(t_vec, traj)
        sols[local_i] = build_solution(
            prob, alg, t_vec, traj;
            interp, retcode, stats = nothing
        )
    end
    return sols
end

function _unpack_dense_solutions!(
        sols, irange, saved_data_buf, saved_times, nout, prob, alg, ::Type{T};
        status = nothing, nout_out = nothing, t0 = nothing, dt = nothing,
        save_end::Bool = false
    ) where {T}
    n_particles = length(irange)
    _foreach_chunked(n_particles) do local_i
        i = irange[local_i]
        emitted = nout
        terminated = false
        n_completed = 0
        if status !== nothing
            terminated, n_completed, emitted = _term_info(status, nout_out, i, nout)
        end

        traj = Vector{SVector{6, T}}(undef, emitted)
        for j in 1:emitted
            traj[j] = SVector{6, T}(
                saved_data_buf[i, 1, j],
                saved_data_buf[i, 2, j],
                saved_data_buf[i, 3, j],
                saved_data_buf[i, 4, j],
                saved_data_buf[i, 5, j],
                saved_data_buf[i, 6, j]
            )
        end
        t_vec = saved_times
        retcode = ReturnCode.Success
        if terminated
            retcode = ReturnCode.Terminated
            t_vec = saved_times[1:emitted]
            # The last state of a trace that stopped early sits at the last node
            # inside the domain, not at the nominal end of the time span.
            if save_end && !isempty(t_vec)
                t_vec[end] = t0 + n_completed * dt
            end
        end
        interp = LinearInterpolation(t_vec, traj)
        sols[local_i] = build_solution(
            prob, alg, t_vec, traj;
            interp, retcode, stats = nothing
        )
    end
    return sols
end

@inline function _call_prob_func(prob_func, prob, i::Int, repeat::Int, seed)
    ctx = (sim_id = i, repeat = repeat)
    return prob_func(prob, ctx)
end


@inline function _init_particle!(
        xv_init::AbstractMatrix, prob, prob_func, i::Int, seed
    )
    new_prob = _call_prob_func(prob_func, prob, i, 0, seed)
    u0_i = new_prob.u0
    @inbounds for c in 1:6
        xv_init[i, c] = u0_i[c]
    end
    return nothing
end

function _init_particles!(
        xv_init::AbstractMatrix{T}, prob, prob_func, n_particles::Int, seed
    ) where {T}
    if n_particles == 1
        _init_particle!(xv_init, prob, prob_func, 1, seed)
    else
        _foreach_chunked(n_particles) do i
            _init_particle!(xv_init, prob, prob_func, i, seed)
        end
    end
    return xv_init
end

@inline function _store_state!(::Nothing, traj, i, idx, y)
    traj[idx + 1] = y
    return idx + 1
end

@inline function _store_state!(raw_out::AbstractArray, traj, i, idx, y)
    @inbounds for c in 1:6
        raw_out[i, c, idx + 1] = y[c]
    end
    return idx + 1
end

@inline _has_isoutside(::typeof(ODE_DEFAULT_ISOUTOFDOMAIN)) = false
@inline _has_isoutside(::Nothing) = false
@inline _has_isoutside(::Any) = true

function _solve_cpu_chunk_unchecked!(
        sols, irange, xv_current, plan, dt, tspan,
        p_host, alg, nout, nt, save_start, save_end, save_everystep,
        saved_times, prob, ::Type{T}; raw_out = nothing
    ) where {T}
    nsave = length(plan.times)
    for (local_i, i) in enumerate(irange)
        traj = raw_out === nothing ? Vector{SVector{6, T}}(undef, nout) : nothing
        idx = 0
        r = SVector{3, T}(xv_current[i, 1], xv_current[i, 2], xv_current[i, 3])
        v = SVector{3, T}(xv_current[i, 4], xv_current[i, 5], xv_current[i, 6])

        if save_start
            idx = _store_state!(raw_out, traj, i, idx, vcat(r, v))
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
                        idx = _store_state!(
                            raw_out, traj, i, idx,
                            saveat_interpolate(t_prev, y_prev, t_current, y_cur, t_target)
                        )
                        isave += 1
                    end
                end
            elseif save_everystep
                if it < nt || save_end
                    v_node = update_velocity_node(
                        v_half, r, dt, t_current, p_host, alg
                    )
                    idx = _store_state!(raw_out, traj, i, idx, vcat(r, v_node))
                end
            end

            r_prev = r
            v_half_prev = v_half
            t_prev = t_current
        end

        if save_end && (!use_saveat(plan) && !save_everystep)
            v_node = update_velocity_node(v_half, r, dt, tspan[2], p_host, alg)
            idx = _store_state!(raw_out, traj, i, idx, vcat(r, v_node))
        elseif save_end && use_saveat(plan) && idx < nout
            v_node = update_velocity_node(v_half, r, dt, tspan[2], p_host, alg)
            idx = _store_state!(raw_out, traj, i, idx, vcat(r, v_node))
        end

        if raw_out === nothing
            interp = LinearInterpolation(saved_times, traj)
            sols[local_i] = build_solution(
                prob, alg, saved_times, traj;
                interp, retcode = ReturnCode.Success, stats = nothing
            )
        end
    end
    return raw_out === nothing ? sols : raw_out
end

function _solve_cpu_chunk_checked!(
        sols, irange, xv_current, plan, dt, tspan,
        p_host, alg, nout, nt, save_start, save_end, save_everystep,
        saved_times, prob, ::Type{T};
        raw_out = nothing, isoutside, term_out = nothing, nout_out = nothing
    ) where {T}
    nsave = length(plan.times)
    for (local_i, i) in enumerate(irange)
        traj = raw_out === nothing ? Vector{SVector{6, T}}(undef, nout) : nothing
        idx = 0
        r = SVector{3, T}(xv_current[i, 1], xv_current[i, 2], xv_current[i, 3])
        v = SVector{3, T}(xv_current[i, 4], xv_current[i, 5], xv_current[i, 6])

        if save_start
            idx = _store_state!(raw_out, traj, i, idx, vcat(r, v))
        end

        v_half = update_velocity_half(v, r, dt, tspan[1], p_host, alg)
        isave = 1
        r_prev = r
        v_half_prev = v_half
        t_prev = tspan[1]
        rejected = 0

        for it in 1:nt
            t_current = tspan[1] + it * dt
            r, v_half = advance_boris(v_half_prev, r_prev, dt, t_prev, p_host, alg)

            v_node = update_velocity_node(v_half, r, dt, t_current, p_host, alg)
            if isoutside(vcat(r, v_node), p_host, t_current)
                r, v_half = r_prev, v_half_prev
                rejected = it
                break
            end

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
                        idx = _store_state!(
                            raw_out, traj, i, idx,
                            saveat_interpolate(t_prev, y_prev, t_current, y_cur, t_target)
                        )
                        isave += 1
                    end
                end
            elseif save_everystep
                if it < nt || save_end
                    idx = _store_state!(raw_out, traj, i, idx, vcat(r, v_node))
                end
            end

            r_prev = r
            v_half_prev = v_half
            t_prev = t_current
        end

        if rejected > 0 && term_out !== nothing
            @inbounds term_out[i] = rejected
        end

        t_end = rejected > 0 ? t_prev : tspan[2]
        fresh_end = rejected == 0 ||
            idx == 0 || plan.dir * (t_end - saved_times[idx]) > 0

        if save_end && (!use_saveat(plan) && !save_everystep) && fresh_end
            v_node = update_velocity_node(v_half, r, dt, t_end, p_host, alg)
            idx = _store_state!(raw_out, traj, i, idx, vcat(r, v_node))
        elseif save_end && use_saveat(plan) && idx < nout && fresh_end
            v_node = update_velocity_node(v_half, r, dt, t_end, p_host, alg)
            idx = _store_state!(raw_out, traj, i, idx, vcat(r, v_node))
        end

        if rejected > 0 && nout_out !== nothing
            @inbounds nout_out[i] = idx
        end

        if raw_out === nothing
            t_vec = saved_times
            traj_vec = traj
            retcode = ReturnCode.Success
            if rejected > 0
                retcode = ReturnCode.Terminated
                t_vec = saved_times[1:idx]
                traj_vec = traj[1:idx]
                if save_end && !isempty(t_vec)
                    t_vec[end] = t_end
                end
            end
            interp = LinearInterpolation(t_vec, traj_vec)
            sols[local_i] = build_solution(
                prob, alg, t_vec, traj_vec; interp, retcode, stats = nothing
            )
        end
    end
    return raw_out === nothing ? sols : raw_out
end

"""
    _solve_cpu_chunk!(sols, irange, xv_current, plan, dt, tspan, p_host, alg, nout, nt,
        save_start, save_end, save_everystep, saved_times, prob, T; raw_out = nothing,
        isoutside = ODE_DEFAULT_ISOUTOFDOMAIN, term_out = nothing, nout_out = nothing)

Advance every particle in `irange` on the host. When no out-of-domain callback is provided,
the fast unchecked path runs without any boundary overhead.
"""
function _solve_cpu_chunk!(
        sols, irange, xv_current, plan, dt, tspan,
        p_host, alg, nout, nt, save_start, save_end, save_everystep,
        saved_times, prob, ::Type{T};
        raw_out = nothing, isoutside::F = ODE_DEFAULT_ISOUTOFDOMAIN,
        term_out = nothing, nout_out = nothing
    ) where {T, F}
    if _has_isoutside(isoutside)
        return _solve_cpu_chunk_checked!(
            sols, irange, xv_current, plan, dt, tspan,
            p_host, alg, nout, nt, save_start, save_end, save_everystep,
            saved_times, prob, T; raw_out, isoutside, term_out, nout_out
        )
    else
        return _solve_cpu_chunk_unchecked!(
            sols, irange, xv_current, plan, dt, tspan,
            p_host, alg, nout, nt, save_start, save_end, save_everystep,
            saved_times, prob, T; raw_out
        )
    end
end

"""
    _solve_raw!(backend, xv_current, xv_init, irange; dt, tspan, plan, workgroup_size,
        p_gpu, p_host, alg, nout, nt, save_start, save_end, save_everystep, prob)

Advance the ensemble and return `(; u, t, terminated)`, where `u` is a `(N, 6, nout)`
array of saved states and `t` the matching times. This is the path `raw_output = true`
takes: building one `ODESolution` per particle dominates the host cost of large
ensembles, and the saved values, including `saveat` interpolation, are identical to
the solution path.

The returned `terminated` flags which traces stopped early because they left the domain.
Their trailing slots hold the last state inside the domain, since `t` is shared by the
whole ensemble and cannot be shortened per particle.
"""
function _solve_raw!(
        backend::Backend, xv_current, xv_init, irange;
        dt, tspan, plan, workgroup_size, p_gpu, p_host, alg, nout, nt,
        save_start, save_end, save_everystep, prob, pin::Bool = false,
        isoutside::F = ODE_DEFAULT_ISOUTOFDOMAIN
    ) where {F}
    T = eltype(xv_current)
    n_total = size(xv_current, 1)
    offset = irange.start - 1
    n_particles = length(irange)
    check = _has_isoutside(isoutside)
    isoutside_gpu = check ? adapt_field_to_gpu(isoutside, backend) : nothing
    time_type = typeof(tspan[1] + dt)
    saved_times = _build_saved_times(
        tspan, dt, nt, plan, save_start, save_end, save_everystep, time_type
    )
    dense = use_saveat(plan) || save_everystep
    # With exactly one saved state per particle nothing has to be copied: the final
    # state already lives in `xv_current`, so `u` can alias it.
    alias = !dense && !save_start && nout == 1

    if backend isa CPU
        u = alias ? reshape(xv_current, n_total, 6, 1) : zeros(T, n_total, 6, nout)
        status = check ? zeros(Int, n_total) : nothing
        nout_out = check ? zeros(Int, n_total) : nothing
        chunks = if Threads.nthreads() > 1 && n_particles > 16
            index_chunks(irange; n = Threads.nthreads())
        else
            [irange]
        end
        Threads.@threads for chunk in chunks
            _solve_cpu_chunk!(
                nothing, chunk, xv_current, plan, dt, tspan,
                p_host, alg, nout, nt, save_start, save_end, save_everystep,
                saved_times, prob, T; raw_out = u, isoutside,
                term_out = status, nout_out
            )
        end
        _freeze_trailing!(u, status, nout_out, nout)
        return (; u, t = saved_times, terminated = _terminated_flags(status, n_total))
    end

    if !dense
        if check
            status = KA.zeros(backend, Int, n_particles)
            kernel! = boris_trajectory_checked_kernel!(backend, workgroup_size)
            kernel!(
                xv_current, p_gpu, dt, tspan[1], nt, offset, alg,
                isoutside_gpu, status; ndrange = n_particles
            )
        else
            status = nothing
            kernel! = boris_trajectory_kernel!(backend, workgroup_size)
            kernel!(
                xv_current, p_gpu, dt, tspan[1], nt, offset, alg;
                ndrange = n_particles
            )
        end
        synchronize(backend)

        xv_host = _staging_array(backend, T, (n_total, 6); pin)
        copyto!(xv_host, xv_current)
        status_host = _copy_to_host(status, n_particles)
        if alias
            return (;
                u = reshape(xv_host, n_total, 6, 1), t = saved_times,
                terminated = _terminated_flags(status_host, n_total)
            )
        end
        u = _staging_array(backend, T, (n_total, 6, nout); pin)
        save_start && _copy_column!(u, xv_init, irange, 1)
        save_end && _copy_column!(u, xv_host, irange, nout)
        return (;
            u, t = saved_times, terminated = _terminated_flags(status_host, n_total)
        )
    end

    u = _staging_array(backend, T, (n_total, 6, nout); pin)
    saved_gpu = KA.zeros(backend, T, n_total, 6, nout)
    if check
        status = KA.zeros(backend, Int, n_particles)
        nout_out = KA.zeros(backend, Int, n_particles)
        if use_saveat(plan)
            plan_times_gpu = adapt_field_to_gpu(plan.times, backend)
            nsave = length(plan.times)
            kernel! = boris_saveat_checked_kernel!(backend, workgroup_size)
            kernel!(
                saved_gpu, xv_current, plan_times_gpu, p_gpu, dt,
                tspan[1], nt, nsave, plan.dir, save_start, save_end,
                offset, alg, isoutside_gpu, status, nout_out;
                ndrange = n_particles
            )
        else
            kernel! = boris_everystep_checked_kernel!(backend, workgroup_size)
            kernel!(
                saved_gpu, xv_current, p_gpu, dt,
                tspan[1], nt, save_start, save_end,
                offset, alg, isoutside_gpu, status, nout_out;
                ndrange = n_particles
            )
        end
    else
        status = nothing
        nout_out = nothing
        if use_saveat(plan)
            plan_times_gpu = adapt_field_to_gpu(plan.times, backend)
            nsave = length(plan.times)
            kernel! = boris_saveat_kernel!(backend, workgroup_size)
            kernel!(
                saved_gpu, xv_current, plan_times_gpu, p_gpu, dt,
                tspan[1], nt, nsave, plan.dir, save_start, save_end,
                offset, alg; ndrange = n_particles
            )
        else
            kernel! = boris_everystep_kernel!(backend, workgroup_size)
            kernel!(
                saved_gpu, xv_current, p_gpu, dt,
                tspan[1], nt, save_start, save_end,
                offset, alg; ndrange = n_particles
            )
        end
    end
    synchronize(backend)
    copyto!(u, saved_gpu)

    status_host = _copy_to_host(status, n_particles)
    nout_out_host = _copy_to_host(nout_out, n_particles)
    _freeze_trailing!(u, status_host, nout_out_host, nout)

    return (; u, t = saved_times, terminated = _terminated_flags(status_host, n_total))
end

"""
    _copy_to_host(status, n) -> Union{Nothing, Vector{Int}}

Bring a device-side status array back to the host, or `nothing` when no boundary
check was active.
"""
_copy_to_host(::Nothing, ::Int) = nothing
_copy_to_host(status, n::Int, check::Bool) = check ? _copy_to_host(status, n) : nothing

function _copy_to_host(status::AbstractArray, n::Int)
    host = Vector{Int}(undef, n)
    copyto!(host, status)
    return host
end

@inline function _copy_column!(out, src, irange, j)
    @inbounds for i in irange
        for c in 1:6
            out[i, c, j] = src[i, c]
        end
    end
    return out
end

function _solve_kernel_serial(
        prob, backend::Backend, irange;
        dt, plan, save_start::Bool,
        save_end::Bool, save_everystep::Bool, workgroup_size::Int,
        xv_current, xv_init, is_cpu_accessible::Bool,
        p_gpu, p_host, alg, nout::Int, nt::Int,
        pin::Bool = false, isoutside::F = ODE_DEFAULT_ISOUTOFDOMAIN
    ) where {F}
    (; tspan) = prob
    T = eltype(xv_current)
    n_particles = length(irange)
    offset = irange.start - 1
    check = _has_isoutside(isoutside)
    isoutside_gpu = check ? adapt_field_to_gpu(isoutside, backend) : nothing
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
            saved_times, prob, T; isoutside
        )
    end

    if !use_saveat(plan) && !save_everystep
        if check
            status = KA.zeros(backend, Int, n_particles)
            traj_kernel! = boris_trajectory_checked_kernel!(backend, workgroup_size)
            traj_kernel!(
                xv_current, p_gpu, dt, tspan[1], nt, offset, alg,
                isoutside_gpu, status; ndrange = n_particles
            )
        else
            status = nothing
            traj_kernel! = boris_trajectory_kernel!(backend, workgroup_size)
            traj_kernel!(
                xv_current, p_gpu, dt, tspan[1], nt, offset, alg;
                ndrange = n_particles
            )
        end
        synchronize(backend)

        xv_cpu_end = is_cpu_accessible ? xv_current : zeros(T, size(xv_current, 1), 6)
        if !is_cpu_accessible
            copyto!(xv_cpu_end, xv_current)
        end

        _unpack_endpoint_solutions!(
            sols, irange, xv_init, xv_cpu_end, saved_times, nout, save_start,
            save_end, prob, alg, T;
            status = _copy_to_host(status, n_particles),
            t0 = tspan[1], dt
        )
    else
        n_total = size(xv_current, 1)
        saved_data_gpu = KA.zeros(backend, T, n_total, 6, nout)

        if check
            status = KA.zeros(backend, Int, n_particles)
            nout_out = KA.zeros(backend, Int, n_particles)
            if use_saveat(plan)
                plan_times_gpu = adapt_field_to_gpu(plan.times, backend)
                nsave = length(plan.times)
                kernel! = boris_saveat_checked_kernel!(backend, workgroup_size)
                kernel!(
                    saved_data_gpu, xv_current, plan_times_gpu, p_gpu, dt,
                    tspan[1], nt, nsave, plan.dir, save_start, save_end,
                    offset, alg, isoutside_gpu, status, nout_out;
                    ndrange = n_particles
                )
            else
                kernel! = boris_everystep_checked_kernel!(backend, workgroup_size)
                kernel!(
                    saved_data_gpu, xv_current, p_gpu, dt,
                    tspan[1], nt, save_start, save_end,
                    offset, alg, isoutside_gpu, status, nout_out;
                    ndrange = n_particles
                )
            end
        else
            status = nothing
            nout_out = nothing
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
        end
        synchronize(backend)

        saved_data_buf = if is_cpu_accessible
            saved_data_gpu
        else
            saved_cpu = _staging_array(backend, T, (n_total, 6, nout); pin)
            copyto!(saved_cpu, saved_data_gpu)
            saved_cpu
        end

        _unpack_dense_solutions!(
            sols, irange, saved_data_buf, saved_times, nout, prob, alg, T;
            status = _copy_to_host(status, n_particles),
            nout_out = _copy_to_host(nout_out, n_particles),
            t0 = tspan[1], dt, save_end
        )
    end

    return sols
end

function _prepare_ensemble_solve(
        prob, prob_func, backend::Backend, trajectories::Int, dt::Real,
        plan, save_start::Bool, save_end::Bool, save_everystep::Bool, maxiters::Int;
        seed = nothing, u0 = nothing, pin::Bool = false
    )
    (; tspan, p) = prob
    prob_u0 = prob.u0
    T = eltype(prob_u0)
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

    xv_init = if u0 !== nothing
        # One bulk copy instead of one `prob_func` (and one `remake`) per particle.
        size(u0) == (n_particles, 6) ||
            throw(ArgumentError("u0 must be a $n_particles x 6 matrix"))
        mat = u0 isa Matrix{T} ? u0 : Matrix{T}(u0)
        pin && _pin_array!(backend, mat)
        mat
    elseif is_cpu_accessible
        xv_current
    else
        _staging_array(backend, T, (n_particles, 6); pin)
    end

    u0 === nothing && _init_particles!(xv_init, prob, prob_func, n_particles, seed)

    if xv_current !== xv_init
        copyto!(xv_current, xv_init)
    end
    return (;
        nt, nout, xv_current, xv_init, is_cpu_accessible,
        p_gpu, p_host, tspan = tspan_T, u0 = prob_u0, T, dt = dt_T,
    )
end

function _execute_ensemble_kernel(
        base_prob, prob_func, alg::AbstractBoris, ens::EnsembleKernel;
        dt::Real, trajectories::Int = 1,
        saveat = (),
        save_start::Bool = true, save_end::Bool = true, save_everystep::Bool = true,
        maxiters::Int = 1_000_000,
        seed = nothing,
        u0 = nothing,
        raw_output::Bool = false,
        pin::Union{Nothing, Bool} = nothing,
        isoutside::F = ODE_DEFAULT_ISOUTOFDOMAIN
    ) where {F}
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
    effective_pin = _should_pin(
        backend, pin !== nothing ? pin : ens.pin, raw_output, trajectories
    )
    T = eltype(base_prob.u0)
    dt_T = T(dt)
    tspan_T = (T(base_prob.tspan[1]), T(base_prob.tspan[2]))
    time_type = typeof(tspan_T[1] + dt_T)
    plan = SavingPlan(saveat, tspan_T, time_type)

    (;
        nt, nout, xv_current, xv_init, is_cpu_accessible,
        p_gpu, p_host, tspan,
    ) = _prepare_ensemble_solve(
        base_prob, prob_func, backend, trajectories, dt_T, plan,
        save_start, save_end, save_everystep, maxiters;
        seed, u0, pin = effective_pin
    )

    t0 = time_ns()

    sols = if raw_output
        _solve_raw!(
            backend, xv_current, xv_init, 1:trajectories;
            dt = dt_T, tspan, plan, workgroup_size, p_gpu, p_host, alg,
            nout, nt, save_start, save_end, save_everystep, prob = base_prob,
            pin = effective_pin, isoutside
        )
    elseif backend isa CPU && Threads.nthreads() > 1 && trajectories > 16
        saved_times = _build_saved_times(
            tspan, dt_T, nt, plan, save_start, save_end, save_everystep, time_type
        )
        sample_sol = build_solution(
            base_prob, alg, saved_times, [SVector{6, T}(base_prob.u0)]
        )
        res_sols = Vector{typeof(sample_sol)}(undef, trajectories)

        chunks = index_chunks(1:trajectories; n = Threads.nthreads())
        Threads.@threads for irange in chunks
            chunk_sols = _solve_kernel_serial(
                base_prob, backend, irange;
                dt = dt_T, plan, save_start, save_end, save_everystep, workgroup_size,
                xv_current, xv_init, is_cpu_accessible,
                p_gpu, p_host, alg, nout, nt, pin = effective_pin, isoutside
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
            p_gpu, p_host, alg, nout, nt, pin = effective_pin, isoutside
        )
    end

    elapsed_time = (time_ns() - t0) * 1.0e-9

    return sols, elapsed_time
end

"""
    solve(prob, alg, backend; trajectories, dt, u0 = nothing, raw_output = false)

Trace an ensemble on a hardware backend.

  - `u0`: a `(trajectories, 6)` matrix of initial states. When given, `prob_func` is
    skipped entirely, which avoids one problem copy per particle.
  - `raw_output`: return a NamedTuple `(; u, t, terminated)` instead of an
    `EnsembleSolution`, where `u` is a `(trajectories, 6, nout)` array of the saved
    states and `t` the matching times. No `ODESolution` is built, which is what makes
    large ensembles affordable. `terminated` flags the traces that stopped early
    because they left the domain, whose trailing states are the last inside it.
  - `isoutside`: boundary check `(u, p, t)`. A step that ends outside the domain is
    rolled back and the trace stops at the last state inside it, reported with
    `ReturnCode.Terminated`. On a device the check runs inside the kernel, with the
    parameter container the device fields use.
"""
function SciMLBase.__solve(
        eprob::SciMLBase.AbstractEnsembleProblem, alg::AbstractBoris, ens::EnsembleKernel;
        dt::Real, trajectories::Int = 1,
        saveat = (),
        save_start::Bool = true, save_end::Bool = true, save_everystep::Bool = true,
        maxiters::Int = 1_000_000,
        seed = nothing,
        u0 = nothing,
        raw_output::Bool = false,
        pin::Union{Nothing, Bool} = nothing,
        isoutside::F = ODE_DEFAULT_ISOUTOFDOMAIN,
        kwargs...
    ) where {F}
    base_prob = eprob.prob
    prob_func = eprob.prob_func
    output_func = eprob.output_func
    reduction = eprob.reduction

    sols, elapsed_time = _execute_ensemble_kernel(
        base_prob, prob_func, alg, ens;
        dt, trajectories, saveat, save_start, save_end,
        save_everystep, maxiters, seed, u0, raw_output, pin, isoutside
    )

    raw_output && return sols

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
        workgroup_size::Int = 256, pin::Union{Nothing, Bool} = nothing, kwargs...
    )
    return SciMLBase.solve(
        prob, alg, EnsembleKernel(backend; workgroup_size, pin); kwargs...
    )
end

function SciMLBase.solve(
        prob::SciMLBase.EnsembleProblem,
        alg::AbstractBoris, backend::Backend;
        workgroup_size::Int = 256, pin::Union{Nothing, Bool} = nothing, kwargs...
    )
    return SciMLBase.solve(
        prob, alg, EnsembleKernel(backend; workgroup_size, pin); kwargs...
    )
end

function SciMLBase.solve(
        prob::AbstractODEProblem,
        alg::AbstractBoris, backend::Backend, ::BasicEnsembleAlgorithm;
        workgroup_size::Int = 256, pin::Union{Nothing, Bool} = nothing, kwargs...
    )
    return SciMLBase.solve(
        prob, alg, EnsembleKernel(backend; workgroup_size, pin); kwargs...
    )
end

function SciMLBase.solve(
        prob::SciMLBase.EnsembleProblem,
        alg::AbstractBoris, backend::Backend, ::BasicEnsembleAlgorithm;
        workgroup_size::Int = 256, pin::Union{Nothing, Bool} = nothing, kwargs...
    )
    return SciMLBase.solve(
        prob, alg, EnsembleKernel(backend; workgroup_size, pin); kwargs...
    )
end
