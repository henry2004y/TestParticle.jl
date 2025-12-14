# Unit support

This example shows how to trace charged particles with Unitful units.

TL;DR: Tracing with units is convenient although not as performant as tracing in dimensionless units.

```@example unit
using TestParticle, StaticArrays
using Unitful
using OrdinaryDiffEqVerner
import TestParticle as TP

const Bmag = 1e-8u"T"
const Ω = abs(Unitful.q) * Bmag / Unitful.mp
const t_max = 1 / Ω |> u"s"
const Emag = 1e-8u"V/m"

B(x) = SA[0.0u"T", 0.0u"T", Bmag]
E(x) = SA[Emag, 0.0u"V/m", 0.0u"V/m"]

x0 = [0.0, 0.0, 0.0] * u"m" # [m]
v0 = [0.0, 0.01, 0.0] * Unitful.c0 # [m/s]
u0 = [x0..., v0...]
tspan = (0.0u"s", 2π * t_max) # [s]

param = prepare(E, B, q = Unitful.q, m = Unitful.mp)
prob = ODEProblem(trace!, u0, tspan, param)
sol = solve(prob, Vern9())
```

## Heterogeneous arrays with `ArrayPartition`

`ArrayPartitions` in DiffEq can be used for heterogeneous arrays.
See https://docs.sciml.ai/DiffEqDocs/stable/features/diffeq_arrays for more details.

```@example unit
using RecursiveArrayTools
u0_p = ArrayPartition(x0, v0)
prob_p = ODEProblem(trace!, u0_p, tspan, param)
sol_p = solve(prob_p, Vern9())
```

## Tracing in standard SI units

```@example unit
const Bmag_SI = Bmag / u"T"
const Emag_SI = Emag / u"V/m"
### Solving in SI units
B_SI(x) = SA[0, 0, Bmag_SI]
E_SI(x) = SA[Emag_SI, 0.0, 0.0]

## Initial conditions
x0_SI = [0.0, 0.0, 0.0]
v0_SI = [0.0, 0.01TP.c, 0.0]
u0_SI = [x0_SI..., v0_SI...]
tspan_SI = tspan ./ u"s" # [s]

param_SI = prepare(E_SI, B_SI, species = Proton)
prob_SI = ODEProblem(trace!, u0_SI, tspan_SI, param_SI)
sol_SI = solve(prob_SI, Vern9())
```

## Compare performance

```@example unit
using Chairmarks, Test
@assert sol_p[end][4] ≈ sol[end][4]
println("SI units (unitless): ")
display(@b solve(prob_SI, Vern9()))
println("Partitioned (unitful): ")
display(@b solve(prob_p, Vern9()))
println("Basic (unitful): ")
display(@b solve(prob, Vern9()))
```


## Related API

```@docs; canonical=false
trace!
```