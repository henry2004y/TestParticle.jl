```@meta
CurrentModule = TestParticle
```

# TestParticle.jl

TestParticle.jl is a test particle tracer in a static electromagnetic field.

This package supports charged particle tracing in analytic/numerical relativistic/non-relativistic

* electric and magnetic field;
* body force field.

All tracing are performed in 3D, as is the nature for the fields. For a numerical field, the mesh is constructed with [Meshes.jl](https://github.com/JuliaGeometry/Meshes.jl), and the field is interpolated with the aid of [Interpolations.jl](https://github.com/JuliaMath/Interpolations.jl).
For an analytical field, the user is responsible for providing the function for calculating the field at a given spatial location.
The actual tracing is done through [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl), thanks to the ODE system of the equations of motion.

For all the tracing methods, we provide both an inplace version (with `!` at the end of the function name) and a non-inplace version using StaticArrays. The non-inplace version requires the initial conditions to be static a static vector. Use them at your convenience.

The single particle motions are the basics in understanding the test particle method. Check out [Single-Particle Motions](https://henry2004y.github.io/KeyNotes/contents/single.html) for more complete maths.

## Installation

```julia
julia> ]
pkg> add TestParticle
```

## Usage

It would be better to understand the basic workflow of [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl) before digging into TestParticle.jl. All we are doing here can be concluded as contructing the ODE system from Newton's 2nd law and preparing the field/particle data. Check more in [examples](@ref).

## Acknowledgement

Nothing can be done such easily without the support of the Julia community. We appreciate all the contributions from developers around the world.
