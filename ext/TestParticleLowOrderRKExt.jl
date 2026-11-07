module TestParticleLowOrderRKExt

# RK4 backing for the deprecated GCRK4 stub.

using TestParticle
import TestParticle: _rk4_alg

using OrdinaryDiffEqLowOrderRK: RK4

_rk4_alg(::Val{:default}) = RK4()

end # module TestParticleLowOrderRKExt

