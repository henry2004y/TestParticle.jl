# Native GC solver driven by Vern7.
#
# The adaptive step itself comes from `Vern7`, which lives in OrdinaryDiffEqVerner
# and is deliberately not a dependency of this package: importing it costs about
# 0.4 s, a third of this package's own load time, and most users never trace a
# guiding center. The two hooks below are implemented by
# `TestParticleVernerExt` once that package is loaded.
#
# The parts that stay local keep the native ensemble contract unaffected: the
# initial step estimated from the gyroperiod, `isoutside` domain checks that roll
# the last step back, linear interpolation between bracketing steps for `saveat`,
# the extra saved columns, and the per-trajectory seeded RNG.

@noinline function _verner_missing()
    throw(
        ArgumentError(
            "GCVern7 and the guiding-center mode of AdaptiveHybrid need " *
            "OrdinaryDiffEqVerner. Run `using OrdinaryDiffEqVerner` first, or use " *
            "GCRK4() for a solver that this package implements itself."
        )
    )
end

"""
    init_gc_integrator(xv, tspan, p, dt, abstol, reltol, maxiters)

Integrator for one guiding-center trajectory, set up to take single adaptive steps
without saving anything itself. The state is the out-of-place `SVector{4}` form
[`trace_gc`](@ref) expects and `dt` is only the proposed first step.

Implemented by `TestParticleVernerExt` for `SVector{4}` states; the fallback here
reports the missing dependency.
"""
init_gc_integrator(xv, tspan, p, dt, abstol, reltol, maxiters) = _verner_missing()

"""
    gc_step!(integrator, isoutside, p) -> ReturnCode

Advance `integrator` by one adaptive step and run the domain check on the state it
landed on. Returns `ReturnCode.Success` once the trace has advanced inside the
domain, `Terminated` if `isoutside` rejected the new state, or the failure code the
integrator reported.

`Terminated` is checked first on purpose: when a step blows up, the domain
predicate gets to classify it, so a predicate that also watches for NaN keeps its
say over a solver level failure code.
"""
gc_step!(integrator, isoutside, p) = _verner_missing()

"""
    _gc_initial_dt(p, xv, t, ::Type{Tt}, dt_initial) -> Tt

Step size to start from: the given one, else half a gyroperiod at `xv`, else a
small fallback where there is no rotation to measure.
"""
function _gc_initial_dt(p, xv, t, ::Type{Tt}, dt_initial) where {Tt}
    isnothing(dt_initial) || return Tt(dt_initial)

    q2m = p.q2m
    Bmag = norm(p.B(get_x(xv), t))
    if Bmag == 0 || q2m == 0
        return Tt(1.0e-6)
    end
    return Tt(0.5 * 2π / abs(q2m * Bmag))
end

"""
Apply Vern7 method for particles with index in `irange`.
"""
function _vern7!(
        sols, prob, irange, alg, plan, dt_initial, isoutside,
        save_start, save_end, save_everystep, maxiters,
        ::Val{SaveFields}, ::Val{SaveWork}, seed = nothing
    ) where {SaveFields, SaveWork}
    (; tspan, p, u0) = prob
    (; abstol, reltol) = alg
    T = eltype(u0)
    Tt = typeof(tspan[1] + one(T))

    vars_dim = 4
    if SaveFields
        vars_dim += 6
    end
    if SaveWork
        vars_dim += 4
    end

    @inbounds for i in irange
        traj = SVector{vars_dim, T}[]
        tsave = Tt[]

        rng = isnothing(seed) ? default_rng() : Xoshiro(_rng_seed(seed, i))
        new_prob = prob.prob_func(prob, EnsembleContext(i, 1, 0, nothing, rng, seed))
        xv = SVector{4, T}(new_prob.u0)

        t = Tt(tspan[1])

        integrator = init_gc_integrator(
            xv, (Tt(tspan[1]), Tt(tspan[2])), p,
            _gc_initial_dt(p, xv, t, Tt, dt_initial), abstol, reltol, maxiters
        )

        if save_start
            push!(traj, _prepare_saved_data_gc(xv, p, t, Val(SaveFields), Val(SaveWork)))
            push!(tsave, t)
        end

        steps = 0
        nsave = length(plan.times)
        isave = 1
        retcode = ReturnCode.Success
        while t < tspan[2] && steps < maxiters
            code = gc_step!(integrator, isoutside, p)
            if code != ReturnCode.Success
                retcode = code
                break
            end

            t_prev = t
            xv_prev = xv
            t = Tt(integrator.t)
            xv = SVector{4, T}(integrator.u)

            if use_saveat(plan)
                while isave <= nsave && saveat_reached(plan.times[isave], t, plan.dir)
                    t_target = plan.times[isave]
                    push!(
                        traj, _prepare_saved_data_gc(
                            saveat_interpolate(t_prev, xv_prev, t, xv, t_target),
                            p, t_target, Val(SaveFields), Val(SaveWork)
                        )
                    )
                    push!(tsave, t_target)
                    isave += 1
                end
            elseif save_everystep
                push!(
                    traj, _prepare_saved_data_gc(
                        xv, p, t, Val(SaveFields), Val(SaveWork)
                    )
                )
                push!(tsave, t)
            end

            steps += 1
        end

        # Only reaching the end of the time span counts as a completed solve.
        if retcode == ReturnCode.Success && t < tspan[2]
            retcode = ReturnCode.MaxIters
        end

        if save_end && (isempty(tsave) || tsave[end] != t)
            push!(traj, _prepare_saved_data_gc(xv, p, t, Val(SaveFields), Val(SaveWork)))
            push!(tsave, t)
        end

        t_final = tsave
        u_final = traj
        interp = LinearInterpolation(t_final, u_final)
        stats = nothing

        sols[i] = build_solution(prob, alg, t_final, u_final; interp, retcode, stats)
    end
    return
end
