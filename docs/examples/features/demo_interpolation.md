# Field Interpolation

A robust field interpolation is the prerequisite for pushing particles.
This example demonstrates the construction of scalar/vector field interpolators for Cartesian/Spherical grids.
If the field is analytic, you can directly pass the generated function to [`prepare`](@ref).

```@example interp
import TestParticle as TP
using StaticArrays
using Chairmarks

function setup_spherical_field()
   r = logrange(0.1, 10.0, length = 16)
   r_uniform = range(0.1, 10.0, length = 16)
   θ = range(0, π, length = 16)
   ϕ = range(0, 2π, length = 16)

   B₀ = 1e-8 # [nT]
   B = zeros(3, length(r), length(θ), length(ϕ)) # vector
   A = zeros(length(r), length(θ), length(ϕ)) # scalar

   for (iθ, θ_val) in enumerate(θ)
      sinθ, cosθ = sincos(θ_val)
      B[1, :, iθ, :] .= B₀ * cosθ
      B[2, :, iθ, :] .= -B₀ * sinθ
      A[:, iθ, :] .= B₀ * sinθ
   end

   B_field_nu = TP.getinterp(TP.StructuredGrid, B, r, θ, ϕ)
   A_field_nu = TP.getinterp_scalar(TP.StructuredGrid, A, r, θ, ϕ)
   B_field = TP.getinterp(TP.StructuredGrid, B, r_uniform, θ, ϕ)
   A_field = TP.getinterp_scalar(TP.StructuredGrid, A, r_uniform, θ, ϕ)

   return B_field_nu, A_field_nu, B_field, A_field
end

function setup_cartesian_field()
   x = range(-10, 10, length = 16)
   y = range(-10, 10, length = 16)
   z = range(-10, 10, length = 16)
   B = zeros(3, length(x), length(y), length(z)) # vector
   B[3, :, :, :] .= 10e-9
   A = zeros(length(x), length(y), length(z)) # scalar
   A[:, :, :] .= 10e-9

   B_field = TP.getinterp(B, x, y, z)
   A_field = TP.getinterp_scalar(A, x, y, z)

   return B_field, A_field
end

function setup_cartesian_nonuniform_field()
   x = logrange(0.1, 10.0, length = 16)
   y = range(-10, 10, length = 16)
   z = range(-10, 10, length = 16)
   B = zeros(3, length(x), length(y), length(z)) # vector
   B[3, :, :, :] .= 10e-9
   A = zeros(length(x), length(y), length(z)) # scalar
   A[:, :, :] .= 10e-9

   B_field = TP.getinterp(TP.RectilinearGrid, B, x, y, z)
   A_field = TP.getinterp_scalar(TP.RectilinearGrid, A, x, y, z)

   return B_field, A_field
end

function setup_time_dependent_field()
   x = range(-10, 10, length = 16)
   y = range(-10, 10, length = 16)
   z = range(-10, 10, length = 16)

   # Create two time snapshots
   B0 = zeros(3, length(x), length(y), length(z))
   B0[3, :, :, :] .= 1.0 # Bz = 1 at t=0

   B1 = zeros(3, length(x), length(y), length(z))
   B1[3, :, :, :] .= 2.0 # Bz = 2 at t=1

   times = [0.0, 1.0]

   function loader(i)
       if i == 1
           # For demonstration, we assume we load from disk here
           return TP.getinterp(TP.CartesianGrid, B0, x, y, z)
       elseif i == 2
           return TP.getinterp(TP.CartesianGrid, B1, x, y, z)
       else
           error("Index out of bounds")
       end
   end

   # B_field_t(x, t)
   B_field_t = TP.LazyTimeInterpolator(times, loader)

   return B_field_t
end

B_sph_nu, A_sph_nu, B_sph, A_sph = setup_spherical_field();
B_car, A_car = setup_cartesian_field();
B_car_nu, A_car_nu = setup_cartesian_nonuniform_field();
B_td = setup_time_dependent_field();
loc = SA[1.0, 1.0, 1.0];
```

## Gridded spherical interpolation

```@repl interp
@be B_sph_nu($loc)
@be A_sph_nu($loc)
```

## Uniform spherical interpolation

```@repl interp
@be B_sph($loc)
@be A_sph($loc)
```

## Uniform Cartesian interpolation

```@repl interp
@be B_car($loc)
@be A_car($loc)
```

## Non-uniform Cartesian interpolation

```@repl interp
@be B_car_nu($loc)
@be A_car_nu($loc)
```

Based on the benchmarks, for the same grid size, gridded interpolation (`StructuredGrid` with non-uniform ranges, `RectilinearGrid`) is 2x slower than uniform mesh interpolation (`StructuredGrid` with uniform ranges, `CartesianGrid`).

## Time-dependent field interpolation

For time-dependent fields, we can use [`LazyTimeInterpolator`](@ref). It takes a list of time points and a loader function that returns a spatial interpolator for a given time index. The interpolator will linearly interpolate between the two nearest time points.

```@repl interp
@be B_td($loc, 0.5)
```

## Related API

```@docs; canonical=false
TestParticle.getinterp
TestParticle.getinterp_scalar
TestParticle.prepare
```
