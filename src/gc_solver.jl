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

"""
    update_dp5(y, param, dt, t)::(dy, E)

Update state using the Dormand-Prince 5(4) method.
Returns dy as the update and E as the error estimate, both as SArrays.
"""
@muladd function update_dp5(y, param, dt, t)
    # Coefficients for DP5
    c2, c3, c4, c5, c6 = 1 / 5, 3 / 10, 4 / 5, 8 / 9, 1.0
    a21 = 1 / 5
    a31, a32 = 3 / 40, 9 / 40
    a41, a42, a43 = 44 / 45, -56 / 15, 32 / 9
    a51, a52, a53, a54 = 19372 / 6561, -25360 / 2187, 64448 / 6561, -212 / 729
    a61, a62, a63, a64, a65 = 9017 / 3168, -355 / 33, 46732 / 5247, 49 / 176, -5103 / 18656
    a71, a72, a73, a74, a75, a76 = 35 / 384, 0.0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84

    # Coefficients for 5th order solution and error estimate
    # b2 = b7 = e2 = 0.0
    b1, b3, b4, b5, b6 = 35 / 384, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84
    e1, e3, e4, e5, e6, e7 = 71 / 57600, -71 / 16695, 71 / 1920, -17253 / 339200, 22 / 525, -1 / 40

    k1 = trace_gc(y, param, t)

    y_tmp = y + dt * (a21 * k1)
    k2 = trace_gc(y_tmp, param, t + c2 * dt)

    y_tmp = y + dt * (a31 * k1 + a32 * k2)
    k3 = trace_gc(y_tmp, param, t + c3 * dt)

    y_tmp = y + dt * (a41 * k1 + a42 * k2 + a43 * k3)
    k4 = trace_gc(y_tmp, param, t + c4 * dt)

    y_tmp = y + dt * (a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4)
    k5 = trace_gc(y_tmp, param, t + c5 * dt)

    y_tmp = y + dt * (a61 * k1 + a62 * k2 + a63 * k3 + a64 * k4 + a65 * k5)
    k6 = trace_gc(y_tmp, param, t + c6 * dt)

    y_tmp = y + dt * (a71 * k1 + a72 * k2 + a73 * k3 + a74 * k4 + a75 * k5 + a76 * k6)
    k7 = trace_gc(y_tmp, param, t + dt)

    # y_{n+1} update
    dy = dt * (b1 * k1 + b3 * k3 + b4 * k4 + b5 * k5 + b6 * k6)
    E = dt * (e1 * k1 + e3 * k3 + e4 * k4 + e5 * k5 + e6 * k6 + e7 * k7)

    return dy, E
end

"""
    update_rk4(y, param, dt, t)

Update state using the RK4 method.
Returns dy as the update SArray.
"""
@muladd function update_rk4(y, param, dt, t)
    k1 = trace_gc(y, param, t)

    y_tmp = y + 0.5 * dt * k1
    k2 = trace_gc(y_tmp, param, t + 0.5 * dt)

    y_tmp = y + 0.5 * dt * k2
    k3 = trace_gc(y_tmp, param, t + 0.5 * dt)

    y_tmp = y + dt * k3
    k4 = trace_gc(y_tmp, param, t + dt)

    dy = (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4)

    return dy
end


"""
    solve(prob::TraceGCProblem; trajectories::Int=1, dt::AbstractFloat,
    savestepinterval::Int=1, isoutofdomain::Function=ODE_DEFAULT_ISOUTOFDOMAIN,
    alg::Symbol=:rk4, abstol=1e-6, reltol=1e-6, maxiters=10000,
    save_fields::Bool=false, save_work::Bool=false)

Trace guiding centers using the RK4 method with specified `prob`.
If `alg` is `:rk45`, uses adaptive time stepping.
"""
function solve(
        prob::TraceGCProblem, ensemblealg::BasicEnsembleAlgorithm = EnsembleSerial();
        trajectories::Int = 1, savestepinterval::Int = 1,
        dt::Union{AbstractFloat, Nothing} = nothing,
        isoutofdomain::Function = ODE_DEFAULT_ISOUTOFDOMAIN,
        save_start::Bool = true, save_end::Bool = true, save_everystep::Bool = true,
        alg::Symbol = :rk4, abstol = 1.0e-6, reltol = 1.0e-6, maxiters = 10000,
        save_fields::Bool = false, save_work::Bool = false
    )
    if alg != :rk4 && alg != :rk45
        @warn "Only :rk4 and :rk45 are supported for native TraceGCProblem currently. Using :rk4."
    end

    return if save_fields
        if save_work
            return _solve(
                ensemblealg, prob, trajectories, dt, savestepinterval, isoutofdomain,
                save_start, save_end, save_everystep, alg, abstol, reltol, maxiters,
                Val(true), Val(true)
            )
        else
            return _solve(
                ensemblealg, prob, trajectories, dt, savestepinterval, isoutofdomain,
                save_start, save_end, save_everystep, alg, abstol, reltol, maxiters,
                Val(true), Val(false)
            )
        end
    else
        if save_work
            return _solve(
                ensemblealg, prob, trajectories, dt, savestepinterval, isoutofdomain,
                save_start, save_end, save_everystep, alg, abstol, reltol, maxiters,
                Val(false), Val(true)
            )
        else
            return _solve(
                ensemblealg, prob, trajectories, dt, savestepinterval, isoutofdomain,
                save_start, save_end, save_everystep, alg, abstol, reltol, maxiters,
                Val(false), Val(false)
            )
        end
    end
end

function _solve(
        ::EnsembleSerial, prob::TraceGCProblem, trajectories, dt, savestepinterval,
        isoutofdomain, save_start, save_end, save_everystep, alg, abstol, reltol, maxiters,
        ::Val{SaveFields}, ::Val{SaveWork}
    ) where {SaveFields, SaveWork}
    sols, nt, nout = _prepare_gc(
        prob, trajectories, dt, savestepinterval,
        save_start, save_end, save_everystep, alg, Val(SaveFields), Val(SaveWork)
    )
    irange = 1:trajectories

    if alg == :rk45
        _rk45!(
            sols, prob, irange, dt, isoutofdomain,
            save_start, save_end, save_everystep, abstol, reltol, maxiters,
            Val(SaveFields), Val(SaveWork)
        )
    else
        _rk4!(
            sols, prob, irange, savestepinterval, dt, nt, nout, isoutofdomain,
            save_start, save_end, save_everystep,
            Val(SaveFields), Val(SaveWork)
        )
    end

    return sols
end

function _solve(
        ::EnsembleThreads, prob::TraceGCProblem, trajectories, dt, savestepinterval,
        isoutofdomain, save_start, save_end, save_everystep, alg, abstol, reltol, maxiters,
        ::Val{SaveFields}, ::Val{SaveWork}
    ) where {SaveFields, SaveWork}
    sols, nt,
        nout = _prepare_gc(
        prob, trajectories, dt, savestepinterval,
        save_start, save_end, save_everystep, alg, Val(SaveFields), Val(SaveWork)
    )

    nchunks = Threads.nthreads()
    Threads.@threads for irange in index_chunks(1:trajectories; n = nchunks)
        if alg == :rk45
            _rk45!(
                sols, prob, irange, dt, isoutofdomain,
                save_start, save_end, save_everystep, abstol, reltol, maxiters,
                Val(SaveFields), Val(SaveWork)
            )
        else
            _rk4!(
                sols, prob, irange, savestepinterval, dt, nt, nout, isoutofdomain,
                save_start, save_end, save_everystep,
                Val(SaveFields), Val(SaveWork)
            )
        end
    end

    return sols
end

function _get_sol_type(prob::TraceGCProblem, dt, alg, ::Val{SaveFields}, ::Val{SaveWork}) where {SaveFields, SaveWork}
    u0 = prob.u0
    tspan = prob.tspan
    dt_guess = isnothing(dt) ? one(eltype(u0)) : dt
    T_t = typeof(tspan[1] + dt_guess)
    t = Vector{T_t}(undef, 0)
    # Force u to be Vector{SVector{4, T}}
    T = eltype(u0)

    n_vars = 4
    if SaveFields
        n_vars += 6
    end
    if SaveWork
        n_vars += 4
    end

    u = Vector{SVector{n_vars, T}}(undef, 0)
    interp = LinearInterpolation(t, u)

    sol = build_solution(prob, alg, t, u; interp = interp)
    return typeof(sol)
end

function _prepare_gc(
        prob::TraceGCProblem, trajectories, dt, savestepinterval,
        save_start, save_end, save_everystep, alg, ::Val{SaveFields}, ::Val{SaveWork}
    ) where {SaveFields, SaveWork}
    ttotal = prob.tspan[2] - prob.tspan[1]
    if isnothing(dt)
        if alg == :rk4
            error("dt must be provided for fixed step solver :rk4")
        end
        nt = 0
    else
        if alg == :rk4 && sign(dt) != sign(ttotal) && ttotal != 0
            error("dt must have the same sign as tspan[2] - tspan[1] for fixed step solver :rk4")
        end
        nt = round(Int, ttotal / dt)
    end

    nout = 0
    if save_start
        nout += 1
    end

    # For :rk45, we initialize nout=0 since we don't know the exact steps.
    if alg == :rk4 && save_everystep
        steps = nt ÷ savestepinterval
        last_is_step = (nt > 0) && (nt % savestepinterval == 0)
        nout += steps
        if !save_end && last_is_step
            nout -= 1
        end
        if save_end && !last_is_step
            nout += 1
        end
    elseif alg == :rk4 && save_end
        nout += 1
    end

    sol_type = _get_sol_type(prob, dt, alg, Val(SaveFields), Val(SaveWork))
    sols = Vector{sol_type}(undef, trajectories)

    return sols, nt, nout
end

@inline function _prepare_saved_data_gc(xv, p, t, ::Val{SaveFields}, ::Val{SaveWork}) where {SaveFields, SaveWork}
    data = xv
    if SaveFields
        # p = (q, q2m, μ, Efunc, Bfunc)
        r = get_x(xv)
        T = eltype(xv)
        E = p[4](r, t)
        B = p[5](r, t)
        data = vcat(data, E, B)
    end
    if SaveWork
        work = get_work_rates_gc(xv, p, t)
        data = vcat(data, work)
    end
    return data
end

"""
Apply RK4 method for particles with index in `irange`.
"""
function _rk4!(
        sols, prob, irange, savestepinterval, dt, nt, nout, isoutofdomain,
        save_start, save_end, save_everystep,
        ::Val{SaveFields}, ::Val{SaveWork}
    ) where {SaveFields, SaveWork}
    (; tspan, p, u0) = prob
    T = eltype(u0)

    vars_dim = 4
    if SaveFields
        vars_dim += 6
    end
    if SaveWork
        vars_dim += 4
    end

    @fastmath @inbounds for i in irange
        traj = SVector{vars_dim, T}[]
        sizehint!(traj, nout)
        tsave = typeof(tspan[1] + dt)[]
        sizehint!(tsave, nout)

        # set initial conditions for each trajectory i
        iout = 0
        new_prob = prob.prob_func(prob, i, false)
        xv = new_prob.u0

        if save_start
            iout += 1
            push!(
                traj,
                _prepare_saved_data_gc(xv, p, tspan[1], Val(SaveFields), Val(SaveWork))
            )
            push!(tsave, tspan[1])
        end

        it = 1
        while it <= nt
            t = tspan[1] + (it - 1) * dt

            dx = update_rk4(xv, p, dt, t)
            xv += dx

            isoutofdomain(xv, p, t + dt) && break

            if save_everystep && (it % savestepinterval == 0)
                iout += 1
                if iout <= nout
                    push!(
                        traj,
                        _prepare_saved_data_gc(
                            xv, p, t + dt, Val(SaveFields), Val(SaveWork)
                        )
                    )
                    push!(tsave, t + dt)
                end
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

        if iout < nout && should_save_final && iout > 0 &&
                (isempty(tsave) || tsave[end] != tspan[2])
            iout += 1
            t_final = (it > nt) ? tspan[2] : (tspan[1] + it * dt)
            push!(
                traj, _prepare_saved_data_gc(
                    xv, p, t_final, Val(SaveFields), Val(SaveWork)
                )
            )
            push!(tsave, t_final)
        end

        alg = :rk4
        t = tsave
        interp = LinearInterpolation(t, traj)
        retcode = ReturnCode.Default
        stats = nothing

        sols[i] = build_solution(prob, alg, t, traj; interp, retcode, stats)
    end

    return
end

"""
Apply RK45 method for particles with index in `irange`.
"""
function _rk45!(
        sols, prob, irange, dt_initial, isoutofdomain,
        save_start, save_end, save_everystep, abstol, reltol, maxiters,
        ::Val{SaveFields}, ::Val{SaveWork},
    ) where {SaveFields, SaveWork}
    (; tspan, p, u0) = prob
    T = eltype(u0)

    safety = 0.9
    max_growth = 5.0
    min_growth = 0.2

    vars_dim = 4
    if SaveFields
        vars_dim += 6
    end
    if SaveWork
        vars_dim += 4
    end

    @fastmath @inbounds for i in irange
        traj = SVector{vars_dim, T}[]
        tsave = typeof(tspan[1] + one(T))[]

        new_prob = prob.prob_func(prob, i, false)
        xv = new_prob.u0

        t = tspan[1]

        if isnothing(dt_initial)
            # Estimate initial step size based on gyroperiod
            # p = (q, q2m, μ, E, B)
            q2m = p[2]
            B_field = p[5]
            R = get_x(xv)
            B_vec = B_field(R, t)
            Bmag = norm(B_vec)

            if Bmag == 0 || q2m == 0
                dt = 1.0e-6
            else
                omega = abs(q2m * Bmag)
                dt = 0.5 * 2π / omega # 0.5 gyroperiod
            end
        else
            dt = dt_initial
        end

        if save_start
            push!(traj, _prepare_saved_data_gc(xv, p, t, Val(SaveFields), Val(SaveWork)))
            push!(tsave, t)
        end

        steps = 0
        while t < tspan[2] && steps < maxiters
            if t + dt > tspan[2]
                dt = tspan[2] - t
            end

            dx, E = update_dp5(xv, p, dt, t)

            error_ratio = 0.0

            y_next = xv + dx
            sum_sq_error = 0.0
            for k in 1:4
                sc = abstol + max(abs(xv[k]), abs(y_next[k])) * reltol
                sum_sq_error += (E[k] / sc)^2
            end
            error_ratio = 0.5 * sqrt(sum_sq_error)

            if error_ratio <= 1.0
                t += dt
                xv = y_next

                isoutofdomain(xv, p, t) && break

                if save_everystep
                    push!(
                        traj, _prepare_saved_data_gc(
                            xv, p, t, Val(SaveFields), Val(SaveWork)
                        )
                    )
                    push!(tsave, t)
                end

                steps += 1
            end

            scale = if error_ratio == 0.0
                max_growth
            else
                safety * (1.0 / error_ratio)^0.2
            end
            scale = max(min_growth, min(scale, max_growth))
            dt *= scale

            dt < 1.0e-14 && break
        end

        if save_end && (isempty(tsave) || tsave[end] != t)
            push!(traj, _prepare_saved_data_gc(xv, p, t, Val(SaveFields), Val(SaveWork)))
            push!(tsave, t)
        end

        alg = :rk45
        t_final = tsave
        u_final = traj
        interp = LinearInterpolation(t_final, u_final)
        retcode = steps >= maxiters ? ReturnCode.MaxIters : ReturnCode.Success
        stats = nothing

        sols[i] = build_solution(prob, alg, t_final, u_final; interp, retcode, stats)
    end
    return
end
