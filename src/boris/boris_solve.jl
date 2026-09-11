# Boris integration through the SciML loop.
#
# The algorithms themselves live in OrdinaryDiffEqBoris. This file only adapts
# the TraceProblem container and TestParticle's keyword surface to the SciML
# interface.
#
# Two ways of choosing the output times are supported. Unless `saveat` is given,
# the integrator saves every accepted step and TestParticle's own rules
# (`save_start`, `save_end`, `save_everystep`) select from that series. With
# `saveat` the times go to SciML, which interpolates inside a step; the Boris
# methods declare linear dense output for that, because a Boris step forms no
# derivative stages to interpolate with.

# The Boris methods never evaluate the right-hand side, so a placeholder is
# enough. Building the problem out of place keeps static initial conditions
# working without an in-place cache.
_boris_rhs(u, p, t) = nothing

function _boris_problem(prob::TraceProblem)
    T = eltype(prob.u0)
    return ODEProblem(_boris_rhs, SVector{6, T}(prob.u0), prob.tspan, prob.p)
end

_boris_isadaptive(alg) = hasfield(typeof(alg), :safety)

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

    if !_boris_isadaptive(alg)
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

# `isoutside` follows TestParticle's `(u, p, t)` convention, whereas a
# DiscreteCallback condition receives `(u, t, integrator)`.
#
# The condition is only tested once a step has been taken, so the offending step
# is rolled back before terminating. That keeps the last saved state inside the
# domain, matching the native loop, which drops the step that leaves it.
function _boris_callback(isoutside::F) where {F}
    isoutside === ODE_DEFAULT_ISOUTOFDOMAIN && return nothing

    condition = (u, t, integrator) -> isoutside(u, integrator.p, t)
    affect! = function (integrator)
        integrator.u = integrator.uprev
        return terminate!(integrator)
    end

    return DiscreteCallback(condition, affect!)
end

function _boris_prob_func(prob::TraceProblem)
    return function (p, ctx)
        new_prob = prob.prob_func(p, ctx)
        T = eltype(new_prob.u0)

        return ODEProblem(
            _boris_rhs, SVector{6, T}(new_prob.u0), new_prob.tspan, new_prob.p
        )
    end
end

"""
    _boris_saved_indices(n, save_start, save_end, save_everystep)

Indices of the every-step series `t[1:n]` that TestParticle's saving rules keep.
"""
function _boris_saved_indices(n, save_start, save_end, save_everystep)
    idxs = Int[]
    save_start && push!(idxs, 1)
    save_everystep && append!(idxs, 2:(n - 1))
    save_end && push!(idxs, n)

    return sort!(unique!(idxs))
end

function _boris_build(
        sol::AbstractODESolution, p, idxs, ::Val{SaveFields}, ::Val{SaveWork}
    ) where {SaveFields, SaveWork}
    tsave = sol.t[idxs]
    u = [
        _prepare_saved_data(sol.u[i], p, sol.t[i], Val(SaveFields), Val(SaveWork))
            for i in idxs
    ]

    return build_solution(
        sol.prob, sol.alg, tsave, u;
        interp = LinearInterpolation(tsave, u), retcode = sol.retcode, stats = nothing
    )
end

function _boris_finalize(
        sol::AbstractODESolution, p, save_start, save_end,
        save_everystep, ::Val{SaveFields}, ::Val{SaveWork}
    ) where {SaveFields, SaveWork}
    idxs = _boris_saved_indices(length(sol.t), save_start, save_end, save_everystep)

    return _boris_build(sol, p, idxs, Val(SaveFields), Val(SaveWork))
end

# With `saveat` the output times are already the ones SciML saved, so there is
# nothing to select; only the field and work columns remain to be appended.
function _boris_finalize(
        sol::AbstractODESolution, p, ::Val{SaveFields}, ::Val{SaveWork}
    ) where {SaveFields, SaveWork}
    return _boris_build(sol, p, eachindex(sol.t), Val(SaveFields), Val(SaveWork))
end

"""
    solve(prob::TraceProblem, alg::AbstractBoris, ensemblealg=EnsembleSerial(); kwargs...)

Trace particles with a Boris method, integrated by the SciML loop.

# Keywords
  - `dt`: time step. Optional for adaptive methods, which otherwise start from
    `safety * 2π / |q B / m|`.
  - `saveat`: times to save at, as a collection or as an interval. The solution
    is interpolated linearly inside a step, which is the dense output these
    methods admit.
  - `isoutside`: boundary check `(u, p, t)`; the trace terminates when it holds.
  - `save_start::Bool=true`, `save_end::Bool=true`, `save_everystep::Bool=true`.
    With `saveat`, `save_start` and `save_end` add the ends of the time span to
    the requested times.
  - `save_fields::Bool=false`: append E and B to every saved state.
  - `save_work::Bool=false`: append the work rates to every saved state.
  - `maxiters::Int=1_000_000`: maximum number of steps.
  - `trajectories::Int=1`, `batch_size::Int`, `seed`: ensemble controls.
"""
@inline function solve(
        prob::TraceProblem, alg::AbstractBoris,
        ensemblealg::BasicEnsembleAlgorithm = EnsembleSerial();
        trajectories::Int = 1,
        saveat = (),
        dt::Union{Nothing, AbstractFloat} = nothing,
        isoutside::F = ODE_DEFAULT_ISOUTOFDOMAIN,
        save_start::Bool = true,
        save_end::Bool = true,
        save_everystep::Bool = true,
        save_fields::Bool = false,
        save_work::Bool = false,
        maxiters::Int = 1_000_000,
        batch_size::Int = _default_batch_size(ensemblealg, trajectories),
        seed::Union{Nothing, Integer} = nothing,
    ) where {F}
    step = dt === nothing ? _boris_initial_dt(prob, alg) : dt
    _boris_check_limits(prob, step, alg, maxiters)

    solve_kwargs = if isempty(saveat)
        (
            dt = step,
            save_start = true,
            save_end = true,
            save_everystep = true,
            maxiters,
            callback = _boris_callback(isoutside),
            dense = false,
        )
    else
        (
            dt = step,
            saveat,
            save_start,
            save_end,
            maxiters,
            callback = _boris_callback(isoutside),
            dense = false,
        )
    end

    ensemble_prob = EnsembleProblem(
        _boris_problem(prob); prob_func = _boris_prob_func(prob)
    )
    ensemble_kwargs = if ensemblealg isa Union{EnsembleSplitThreads, EnsembleDistributed}
        (; batch_size)
    else
        NamedTuple()
    end
    ensemble_kwargs = isnothing(seed) ? ensemble_kwargs : merge(ensemble_kwargs, (; seed))

    elapsed_time = @elapsed esol = SciMLBase.solve(
        ensemble_prob, alg, ensemblealg;
        trajectories, ensemble_kwargs..., solve_kwargs...
    )

    sols = if isempty(saveat)
        [
            _boris_finalize(
                sol, prob.p, save_start, save_end, save_everystep,
                Val(save_fields), Val(save_work)
            ) for sol in esol.u
        ]
    else
        [
            _boris_finalize(sol, prob.p, Val(save_fields), Val(save_work))
                for sol in esol.u
        ]
    end

    return EnsembleSolution(sols, elapsed_time, true)
end
