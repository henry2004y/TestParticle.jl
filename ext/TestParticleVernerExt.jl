module TestParticleVernerExt

# Vern7 backing for the adaptive guiding-center solver. Kept out of the package
# proper so that `using TestParticle` does not pay for OrdinaryDiffEqVerner.

using TestParticle
import TestParticle: init_gc_integrator, gc_step!, trace_gc

using OrdinaryDiffEqVerner: Vern7
using SciMLBase: ODEProblem, ReturnCode
import SciMLBase
using StaticArrays: SVector
using PrecompileTools: @setup_workload, @compile_workload

function init_gc_integrator(
        xv::SVector{4, T}, tspan, p, dt, abstol, reltol, maxiters
    ) where {T}
    prob = ODEProblem(trace_gc, xv, tspan, p)
    return SciMLBase.init(
        prob, Vern7();
        dt, abstol, reltol, maxiters,
        save_start = false, save_end = false, save_everystep = false, dense = false
    )
end

function gc_step!(
        integrator::SciMLBase.DEIntegrator, isoutside::F, p
    ) where {F}
    t_before = integrator.t
    SciMLBase.step!(integrator)

    y = SVector{4}(integrator.u)
    isoutside(y, p, integrator.t) && return ReturnCode.Terminated

    code = SciMLBase.check_error(integrator)
    code != ReturnCode.Success && return code
    # Rejected steps are retried internally, so staying put means the solver
    # stopped without recording why.
    return integrator.t == t_before ? ReturnCode.DtLessThanMin : ReturnCode.Success
end

@setup_workload begin
    x = range(-10, 10, length = 4)
    y = range(-10, 10, length = 6)
    z = range(-10, 10, length = 8)
    B = fill(0.0, 3, length(x), length(y), length(z))
    E = fill(0.0, 3, length(x), length(y), length(z))
    B[3, :, :, :] .= 10.0e-9

    gc_init = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0]
    tspan = (0.0, 1.0)

    @compile_workload begin
        stateinit_gc, param_gc = prepare_gc(gc_init, x, y, z; E, B)
        solve(TraceGCProblem(stateinit_gc, tspan, param_gc), GCVern7())
    end
end

end # module TestParticleVernerExt
