# #  Unit support

# This example shows how to trace charged particles with Unitful units.

# Tracing with units is convenient although not as performant as tracing in dimensionless units.

using TestParticle, StaticArrays
using Unitful
using OrdinaryDiffEqVerner

const Bmag = 1e-8u"T"
const Ω = abs(Unitful.q) * Bmag / Unitful.mp
const t_max = 1 / Ω |> u"s"

B_field(x) = SA[0.0u"T", 0.0u"T", Bmag]
E_field(x) = SA[1e-8u"V/m", 0.0u"V/m", 0.0u"V/m"]

x0 = [0.0, 0.0, 0.0] * u"m" # [m]
v0 = [0.0, 0.01, 0.0] * Unitful.c0 # [m/s]
u0 = [x0..., v0...]
tspan = (0.0u"s", 2π * t_max) # [s]

param = prepare(E_field, B_field, species = User, q = Unitful.q, m = Unitful.mp)
prob = ODEProblem(trace!, u0, tspan, param)
solve(prob, Vern9())
