module TestParticleVernerExt

# Vern7 backing for the default guiding-center solver in AdaptiveHybrid and the
# deprecated GCRK45 stub. Kept out of the package proper so that `using TestParticle`
# does not pay for OrdinaryDiffEqVerner.

using TestParticle
import TestParticle: default_gc_alg, _needs_fsal_init

using OrdinaryDiffEqVerner: Vern7
using PrecompileTools: @setup_workload, @compile_workload

default_gc_alg(::Val{:default}) = Vern7()
_needs_fsal_init(::Vern7) = false

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
        stateinit_gc, param_gc = prepare_gc(
            gc_init;
            E = build_interpolator(E, (x, y, z); bc = TestParticle.FillExtrap(NaN)),
            B = build_interpolator(B, (x, y, z); bc = TestParticle.FillExtrap(NaN))
        )
        solve(TraceGCProblem(stateinit_gc, tspan, param_gc), Vern7())
    end
end

end # module TestParticleVernerExt

