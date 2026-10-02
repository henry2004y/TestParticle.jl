# TraceProblem constructors and SciML interface integration

"""
    TraceProblem(u0, tspan, p; prob_func = DEFAULT_PROB_FUNC)
    TraceProblem(f, u0, tspan, p; prob_func = DEFAULT_PROB_FUNC)

Problem specification for particle tracing. Supports both Boris solvers and standard ODE
solvers (e.g. Tsit5, Vern7). When `f` is not provided, defaults to `trace!` (in-place)
for mutable arrays or `trace` (out-of-place) for `StaticArray`s.
"""
function TraceProblem(u0, tspan, p; prob_func = DEFAULT_PROB_FUNC)
    isinplace = !(u0 isa StaticArray)
    _func = isinplace ? trace! : trace
    _f = ODEFunction{isinplace, DEFAULT_SPECIALIZATION}(_func)
    return TraceProblem{
        typeof(u0), typeof(tspan), isinplace, typeof(p), typeof(_f), typeof(prob_func),
    }(_f, u0, tspan, p, prob_func)
end

function TraceProblem(f::Function, u0, tspan, p; prob_func = DEFAULT_PROB_FUNC)
    isinplace = !(u0 isa StaticArray)
    _f = ODEFunction{isinplace, DEFAULT_SPECIALIZATION}(f)
    return TraceProblem{
        typeof(u0), typeof(tspan), isinplace, typeof(p), typeof(_f), typeof(prob_func),
    }(_f, u0, tspan, p, prob_func)
end

function TraceProblem(
        f::AbstractODEFunction{iip}, u0, tspan, p;
        prob_func = DEFAULT_PROB_FUNC
    ) where {iip}
    return TraceProblem{
        typeof(u0), typeof(tspan), iip, typeof(p), typeof(f), typeof(prob_func),
    }(f, u0, tspan, p, prob_func)
end

function TraceProblem{iip}(; f, u0, tspan, p, prob_func = DEFAULT_PROB_FUNC) where {iip}
    return TraceProblem{
        typeof(u0), typeof(tspan), iip, typeof(p), typeof(f), typeof(prob_func),
    }(f, u0, tspan, p, prob_func)
end

function SciMLBase.remake(
        prob::TraceProblem;
        f = prob.f,
        u0 = prob.u0,
        tspan = prob.tspan,
        p = prob.p,
        prob_func = prob.prob_func,
        kwargs...
    )
    isinplace = !(u0 isa StaticArray)
    _f = if f === prob.f && (prob.f.f === trace! || prob.f.f === trace)
        isinplace ? ODEFunction{true, DEFAULT_SPECIALIZATION}(trace!) :
            ODEFunction{false, DEFAULT_SPECIALIZATION}(trace)
    else
        f
    end
    return TraceProblem{
        typeof(u0), typeof(tspan), isinplace, typeof(p), typeof(_f), typeof(prob_func),
    }(_f, u0, tspan, p, prob_func)
end

function solve(
        prob::TraceProblem, alg::SciMLBase.AbstractODEAlgorithm,
        ensemblealg::BasicEnsembleAlgorithm;
        trajectories::Int = 1,
        seed::Union{Nothing, Integer} = nothing,
        kwargs...
    )
    ensemble_prob = EnsembleProblem(
        prob; prob_func = prob.prob_func,
        safetycopy = prob.prob_func !== DEFAULT_PROB_FUNC
    )

    return solve(ensemble_prob, alg, ensemblealg; trajectories, seed, kwargs...)
end
