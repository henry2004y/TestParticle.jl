# Native GC solver

"""
    TraceGCProblem{uType, tType, isinplace, P, F <: AbstractODEFunction, PF} <: AbstractODEProblem{uType, tType, isinplace}

Problem type for tracing guiding centers.
"""
struct TraceGCProblem{uType, tType, isinplace, P, F <: AbstractODEFunction, PF} <:
    AbstractODEProblem{uType, tType, isinplace}
    f::F
    "initial condition"
    u0::uType
    "time span"
    tspan::tType
    "(q2m, E, B)"
    p::P
    "function for setting initial conditions"
    prob_func::PF
end

function TraceGCProblem(u0, tspan, p; prob_func = DEFAULT_PROB_FUNC)
    _f = ODEFunction{true, DEFAULT_SPECIALIZATION}(x -> nothing) # dummy func
    return TraceGCProblem{
        typeof(u0), typeof(tspan), true, typeof(p), typeof(_f),
        typeof(prob_func),
    }(
        _f,
        u0,
        tspan,
        p,
        prob_func
    )
end

# For remake
function TraceGCProblem{iip}(; f, u0, tspan, p, prob_func) where {iip}
    return TraceGCProblem{
        typeof(u0), typeof(tspan), iip, typeof(p), typeof(f),
        typeof(prob_func),
    }(
        f,
        u0,
        tspan,
        p,
        prob_func
    )
end

struct RK4Method{T, N}
    # intermediate variables used in the solver
    k1::MVector{N, T}
    k2::MVector{N, T}
    k3::MVector{N, T}
    k4::MVector{N, T}
    y_tmp::MVector{N, T}
end

function RK4Method(u0::AbstractArray{T}) where {T}
    N = length(u0)
    k1 = MVector{N, T}(undef)
    k2 = MVector{N, T}(undef)
    k3 = MVector{N, T}(undef)
    k4 = MVector{N, T}(undef)
    y_tmp = MVector{N, T}(undef)

    return RK4Method{T, N}(k1, k2, k3, k4, y_tmp)
end

"""
    update_rk4!(xv, method, param, dt, t)

Update state using the RK4 method.
"""
@muladd function update_rk4!(dy, y, method, param, dt, t)
    (; k1, k2, k3, k4, y_tmp) = method

    # k1 = f(t, y)
    trace_gc!(k1, y, param, t)

    # k2 = f(t + dt/2, y + dt/2 * k1)
    @. y_tmp = y + 0.5 * dt * k1
    trace_gc!(k2, y_tmp, param, t + 0.5 * dt)

    # k3 = f(t + dt/2, y + dt/2 * k2)
    @. y_tmp = y + 0.5 * dt * k2
    trace_gc!(k3, y_tmp, param, t + 0.5 * dt)

    # k4 = f(t + dt, y + dt * k3)
    @. y_tmp = y + dt * k3
    trace_gc!(k4, y_tmp, param, t + dt)

    # y_{n+1} = y_n + dt/6 * (k1 + 2k2 + 2k3 + k4)
    @. dy = (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4)

    return
end

"""
    solve(prob::TraceGCProblem; trajectories::Int=1, dt::AbstractFloat,
    savestepinterval::Int=1, isoutofdomain::Function=ODE_DEFAULT_ISOUTOFDOMAIN)

Trace guiding centers using the RK4 method with specified `prob`.
"""
function solve(
        prob::TraceGCProblem, ensemblealg::BasicEnsembleAlgorithm = EnsembleSerial();
        trajectories::Int = 1, savestepinterval::Int = 1, dt::AbstractFloat,
        isoutofdomain::Function = ODE_DEFAULT_ISOUTOFDOMAIN,
        save_start::Bool = true, save_end::Bool = true, save_everystep::Bool = true,
        alg::Symbol = :rk4
    )
    if alg != :rk4
        @warn "Only :rk4 algorithm is supported for native TraceGCProblem currently. Using :rk4."
    end

    return _solve(
        ensemblealg, prob, trajectories, dt, savestepinterval, isoutofdomain,
        save_start, save_end, save_everystep
    )
end

function _solve(
        ::EnsembleSerial, prob::TraceGCProblem, trajectories, dt, savestepinterval, isoutofdomain,
        save_start, save_end, save_everystep
    )
    sols, nt,
        nout = _prepare_gc(
        prob, trajectories, dt, savestepinterval,
        save_start, save_end, save_everystep
    )
    irange = 1:trajectories
    _rk4!(
        sols, prob, irange, savestepinterval, dt, nt, nout, isoutofdomain,
        save_start, save_end, save_everystep
    )

    return sols
end

function _solve(
        ::EnsembleThreads, prob::TraceGCProblem, trajectories, dt, savestepinterval, isoutofdomain,
        save_start, save_end, save_everystep
    )
    sols, nt,
        nout = _prepare_gc(
        prob, trajectories, dt, savestepinterval,
        save_start, save_end, save_everystep
    )

    nchunks = Threads.nthreads()
    Threads.@threads for irange in index_chunks(1:trajectories; n = nchunks)
        _rk4!(
            sols, prob, irange, savestepinterval, dt, nt, nout, isoutofdomain,
            save_start, save_end, save_everystep
        )
    end

    return sols
end

function _get_sol_type(prob::TraceGCProblem, dt)
    u0 = prob.u0
    tspan = prob.tspan
    T_t = typeof(tspan[1] + dt)
    t = Vector{T_t}(undef, 0)
    # Force u to be Vector{SVector{4, T}}
    T = eltype(u0)
    u = Vector{SVector{4, T}}(undef, 0)
    interp = LinearInterpolation(t, u)
    alg = :rk4

    sol = build_solution(prob, alg, t, u; interp = interp)
    return typeof(sol)
end

function _prepare_gc(
        prob::TraceGCProblem, trajectories, dt, savestepinterval,
        save_start, save_end, save_everystep
    )
    ttotal = prob.tspan[2] - prob.tspan[1]
    nt = round(Int, ttotal / dt) |> abs

    nout = 0
    if save_start
        nout += 1
    end

    if save_everystep
        steps = nt ÷ savestepinterval
        last_is_step = (nt > 0) && (nt % savestepinterval == 0)
        nout += steps
        if !save_end && last_is_step
            nout -= 1
        end
        if save_end && !last_is_step
            nout += 1
        end
    elseif save_end
        nout += 1
    end

    sol_type = _get_sol_type(prob, dt)
    sols = Vector{sol_type}(undef, trajectories)

    return sols, nt, nout
end

"""
Apply RK4 method for particles with index in `irange`.
"""
function _rk4!(
        sols, prob, irange, savestepinterval, dt, nt, nout, isoutofdomain,
        save_start, save_end, save_everystep
    )
    (; tspan, p, u0) = prob
    T = eltype(u0)
    xv = MVector{4, T}(undef)
    dx = MVector{4, T}(undef)

    # We need a reusable method struct for each thread/trajectory?
    # For serial/chunks, we can reuse if defined inside loop or per thread.
    # To be safe and thread-local, we create it here or inside the loop.
    # Since specific types might depend on u0, let's create inside.

    @fastmath @inbounds for i in irange
        traj = Vector{SVector{4, T}}(undef, nout)
        tsave = Vector{typeof(tspan[1] + dt)}(undef, nout)

        # New method instance for each trajectory to avoid race conditions if shared (though strictly local here)
        # Using u0 type for RK4Method
        paramRK4 = RK4Method(MVector{4, T}(undef))

        # set initial conditions for each trajectory i
        iout = 0
        new_prob = prob.prob_func(prob, i, false)
        xv .= new_prob.u0

        if save_start
            iout += 1
            traj[iout] = SVector{4, T}(xv)
            tsave[iout] = tspan[1]
        end

        it = 1
        while it <= nt
            t = tspan[1] + (it - 1) * dt

            update_rk4!(dx, xv, paramRK4, p, dt, t)
            @. xv += dx

            if save_everystep && (it % savestepinterval == 0)
                iout += 1
                if iout <= nout
                    traj[iout] = SVector{4, T}(xv)
                    tsave[iout] = t + dt
                end
            end

            if isoutofdomain(xv, p, t + dt)
                break
            end
            it += 1
        end

        # Handle save_end logic
        final_step = min(it, nt)
        should_save_final = false
        if save_end
            should_save_final = true
        elseif save_everystep && (final_step > 0) && (final_step % savestepinterval == 0)
            should_save_final = true
        end

        if iout < nout && should_save_final && iout > 0 && tsave[iout] != tspan[2]
            # If we haven't saved the very last step yet (and we should)
            # Note: if it ran to completion, t + dt could be tspan[2].
            # If we broke early, we save the last state.
            iout += 1
            traj[iout] = SVector{4, T}(xv)
            tsave[iout] = (it > nt) ? tspan[2] : (tspan[1] + it * dt)
        end

        if iout < nout
            resize!(traj, iout)
            resize!(tsave, iout)
        end

        alg = :rk4
        t = tsave
        interp = LinearInterpolation(t, traj)
        retcode = ReturnCode.Default
        stats = nothing

        sols[i] = build_solution(
            prob, alg, t, traj; interp = interp, retcode = retcode, stats = stats
        )
    end

    return
end
