# Boris integration through the SciML loop.
#
# The algorithms themselves live in OrdinaryDiffEqBoris. This file only adapts
# the TraceProblem container and TestParticle's keyword surface to the SciML
# interface.
#
# The integrator always saves every accepted step. TestParticle's saving rules
# (including `savestepinterval`) are then applied by selecting from that series
# and rebuilding the solution. This avoids asking SciML to interpolate at
# `saveat` points, which would require derivative stages that a Boris step never
# forms.

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
    _boris_saved_indices(n, savestepinterval, save_start, save_end, save_everystep)

Indices of the every-step series `t[1:n]` that TestParticle's saving rules keep.
"""
function _boris_saved_indices(n, savestepinterval, save_start, save_end, save_everystep)
    idxs = Int[]
    save_start && push!(idxs, 1)

    if save_everystep
        if savestepinterval == 1
            append!(idxs, 2:(n - 1))
        else
            # The last step is not eligible, matching a fixed-step count.
            for j in 1:div(n - 2, savestepinterval)
                push!(idxs, 1 + j * savestepinterval)
            end
        end
    end

    save_end && push!(idxs, n)

    return sort!(unique!(idxs))
end

function _boris_finalize(
        sol::AbstractODESolution, p, savestepinterval, save_start, save_end,
        save_everystep, ::Val{SaveFields}, ::Val{SaveWork}
    ) where {SaveFields, SaveWork}
    idxs = _boris_saved_indices(
        length(sol.t), savestepinterval, save_start, save_end, save_everystep
    )

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

"""
    solve(prob::TraceProblem, alg::AbstractBoris, ensemblealg=EnsembleSerial(); kwargs...)

Trace particles with a Boris method, integrated by the SciML loop.

# Keywords
  - `dt`: time step. Optional for adaptive methods, which otherwise start from
    `safety * 2π / |q B / m|`.
  - `savestepinterval::Int=1`: save every `savestepinterval`-th step. Deprecated
    in favour of `saveat`.
  - `isoutside`: boundary check `(u, p, t)`; the trace terminates when it holds.
  - `save_start::Bool=true`, `save_end::Bool=true`, `save_everystep::Bool=true`.
  - `save_fields::Bool=false`: append E and B to every saved state.
  - `save_work::Bool=false`: append the work rates to every saved state.
  - `maxiters::Int=1_000_000`: maximum number of steps.
  - `trajectories::Int=1`, `batch_size::Int`, `seed`: ensemble controls.
"""
@inline function solve(
        prob::TraceProblem, alg::AbstractBoris,
        ensemblealg::BasicEnsembleAlgorithm = EnsembleSerial();
        trajectories::Int = 1,
        savestepinterval::Int = 1,
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

    solve_kwargs = (
        dt = step,
        save_start = true,
        save_end = true,
        save_everystep = true,
        maxiters,
        callback = _boris_callback(isoutside),
        dense = false,
    )

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

    sols = [
        _boris_finalize(
            sol, prob.p, savestepinterval, save_start, save_end, save_everystep,
            Val(save_fields), Val(save_work)
        ) for sol in esol.u
    ]

    return EnsembleSolution(sols, elapsed_time, true)
end
