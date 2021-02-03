```@meta
CurrentModule = TestParticle
```

# TestParticle

Test particle tracing in a static electromagnetic field.

This package supports charged particle tracing in
* an analytic electric and magnetic field;
* a numeric electric and magnetic field.

All tracing are done in 3D. For a numerical field, the mesh is constructed with [Meshes.jl](https://github.com/JuliaGeometry/Meshes.jl), and the field is interpolated with the aid of [Interpolations.jl](https://github.com/JuliaMath/Interpolations.jl).
For an analytical field, the user is responsible for providing the function for calculating the field at a given spatial location.
The actual tracing is done through [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl).

Nothing can be done such easily without the support of the Julia community. I appreciate all the contributions from developers around the world.