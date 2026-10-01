# Boris integration through the SciML loop. The algorithms live in
# OrdinaryDiffEqBoris; this file adapts `TraceProblem` and TestParticle's keyword
# surface to the SciML interface. Only one trajectory is adapted here: several of
# them are a SciML `EnsembleProblem`, see the Boris tutorial.

_boris_rhs(u, p, t) = nothing

function _boris_problem(prob::TraceProblem)
    T = eltype(prob.u0)
    return ODEProblem(_boris_rhs, SVector{6, T}(prob.u0), prob.tspan, prob.p)
end

function _boris_initial_dt(prob::TraceProblem, alg)
    q2m = get_q2m(prob.p)
    u0 = prob.u0
    r0 = SVector{3}(u0[1], u0[2], u0[3])
    Bmag = norm(get_BField(prob.p)(r0, prob.tspan[1]))
    tdir = sign(prob.tspan[2] - prob.tspan[1])

    return tdir * (2π * alg.safety) / (abs(q2m) * Bmag)
end

function _boris_check_limits(prob::TraceProblem, dt, alg, maxiters)
    if abs(dt) < 10 * eps(typeof(dt))
        throw(
            ArgumentError(
                "time step dt is too small, violating min_dt = 10 * eps(typeof(dt))"
            )
        )
    end

    if !isadaptive(alg)
        ttotal = prob.tspan[2] - prob.tspan[1]
        nt = abs(round(Int, ttotal / dt))
        if nt > maxiters
            throw(
                ArgumentError("number of iterations nt ($nt) exceeds maxiters ($maxiters)")
            )
        end
    end

    return
end

# `isoutside` takes `(u, p, t)`; a DiscreteCallback condition takes
# `(u, t, integrator)`. The step that leaves the domain is rolled back, so the
# last saved state stays inside it.
function _boris_callback(isoutside::F) where {F}
    isoutside === ODE_DEFAULT_ISOUTOFDOMAIN && return nothing

    condition = (u, t, integrator) -> isoutside(u, integrator.p, t)
    affect! = function (integrator)
        integrator.u = integrator.uprev
        return terminate!(integrator)
    end

    return DiscreteCallback(condition, affect!)
end

# With `saveat` the saved times are already the requested ones, so there is
# nothing to select; only the field and work columns remain to be appended.
_boris_finalize(sol, p, ::Val{false}, ::Val{false}) = sol

function _boris_finalize(
        sol, p, ::Val{SaveFields}, ::Val{SaveWork}
    ) where {SaveFields, SaveWork}
    u = [
        _prepare_saved_data(sol.u[i], p, sol.t[i], Val(SaveFields), Val(SaveWork))
            for i in eachindex(sol.t)
    ]

    return build_solution(
        sol.prob, sol.alg, sol.t, u;
        interp = LinearInterpolation(sol.t, u), retcode = sol.retcode, stats = nothing
    )
end

"""
    solve(prob::TraceProblem, alg::AbstractBoris; kwargs...)

Trace one particle with a Boris method through the SciML loop and return an
`ODESolution`. Trace several by wrapping the problem in a SciML
`EnsembleProblem`, which `solve(prob, alg, ensemblealg)` does for convenience.
`trajectories` belongs to the ensemble; it is not a keyword here.

# Keywords
  - `dt`: time step. Optional for adaptive methods, which otherwise start from
    `safety * 2π / |q B / m|`.
  - `saveat`: times to save at, as a collection or as an interval. The state is
    interpolated linearly inside a step, the dense output these methods admit.
  - `isoutside`: boundary check `(u, p, t)`; the trace terminates when it holds.
  - `save_start::Bool=true`, `save_end::Bool=true`, `save_everystep::Bool=true`.
    With `saveat`, `save_start` and `save_end` add the ends of the time span to
    the requested times.
  - `save_fields::Bool=false`: append E and B to every saved state.
  - `save_work::Bool=false`: append the work rates to every saved state.
  - `maxiters::Int=1_000_000`: maximum number of steps.
  - `rng`, `seed`: accepted and ignored; a SciML ensemble passes them to every
    trajectory it solves.
"""
function solve(
        prob::TraceProblem, alg::AbstractBoris;
        saveat = (),
        dt::Union{Nothing, AbstractFloat} = nothing,
        isoutside::F = ODE_DEFAULT_ISOUTOFDOMAIN,
        save_start::Bool = true,
        save_end::Bool = true,
        save_everystep::Bool = true,
        save_fields::Bool = false,
        save_work::Bool = false,
        maxiters::Int = 1_000_000,
        rng = nothing,
        seed = nothing,
        kwargs...
    ) where {F}
    step = dt === nothing ? _boris_initial_dt(prob, alg) : dt
    _boris_check_limits(prob, step, alg, maxiters)

    ode_prob = _boris_problem(prob)
    callback = _boris_callback(isoutside)

    # `save_everystep` defaults to `isempty(saveat)` in SciML, so it is left out
    # when `saveat` is given: setting it would add every step to the requested
    # times instead of selecting from them. The branch is written out rather than
    # splatted so that each call is specialized on its own keywords.
    sol = if isempty(saveat)
        solve(
            ode_prob, alg; dt = step, save_start, save_end, save_everystep,
            maxiters, callback, dense = false, kwargs...
        )
    else
        solve(
            ode_prob, alg; dt = step, saveat, save_start, save_end,
            maxiters, callback, dense = false, kwargs...
        )
    end

    return _boris_finalize(sol, prob.p, Val(save_fields), Val(save_work))
end

"""
    solve(prob::TraceProblem, alg::AbstractBoris, ensemblealg; kwargs...)

Trace `trajectories` particles with a Boris method and return an
`EnsembleSolution`. A shorthand for

```julia
solve(EnsembleProblem(prob; prob_func = prob.prob_func), alg, ensemblealg;
    trajectories, kwargs...)
```

which is the form to use for anything it does not cover, such as an
`output_func`, a `reduction`, or a `prob_func` other than the one carried by
`prob`.

# Keywords
  - `trajectories::Int=1`: number of particles.
  - `seed`: master seed. SciML derives one reproducible generator per trajectory
    from it and hands it to `prob_func` as `ctx.rng`.
  - `batch_size`, `pmap_batch_size`: SciML's own ensemble keywords, controlling
    how trajectories are grouped before the reduction and how many are sent to a
    worker at a time.
  - Every keyword of the single trajectory `solve` above.
"""
function solve(
        prob::TraceProblem, alg::AbstractBoris,
        ensemblealg::BasicEnsembleAlgorithm;
        trajectories::Int = 1,
        seed::Union{Nothing, Integer} = nothing,
        kwargs...
    )
    # Copy the problem per trajectory only when `prob_func` may mutate it, the
    # rule SciML applies to an `EnsembleProblem` built from a custom `prob_func`.
    ensemble_prob = EnsembleProblem(
        prob; prob_func = prob.prob_func,
        safetycopy = prob.prob_func !== DEFAULT_PROB_FUNC
    )

    return solve(ensemble_prob, alg, ensemblealg; trajectories, seed, kwargs...)
end
