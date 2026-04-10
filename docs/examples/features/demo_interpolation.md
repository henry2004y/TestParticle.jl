# Field Interpolation

A robust field interpolation is the prerequisite for pushing particles.
This example demonstrates the construction of scalar/vector field interpolators for Cartesian/Spherical grids.
If the field is analytic, you can directly pass the generated function to [`prepare`](@ref).

```@example interp
import TestParticle as TP
using StaticArrays
using Chairmarks

function setup_spherical_field(ns = 16)
   r = logrange(0.1, 10.0, length = ns)
   r_uniform = range(0.1, 10.0, length = ns)
   θ = range(0, π, length = ns)
   ϕ = range(0, 2π, length = ns)

   B₀ = 1e-8 # [nT]
   B = zeros(3, length(r), length(θ), length(ϕ)) # vector
   A = zeros(length(r), length(θ), length(ϕ)) # scalar

   for (iθ, θ_val) in enumerate(θ)
      sinθ, cosθ = sincos(θ_val)
      B[1, :, iθ, :] .= B₀ * cosθ
      B[2, :, iθ, :] .= -B₀ * sinθ
      A[:, iθ, :] .= B₀ * sinθ
   end

   B_field_nu = TP.build_interpolator(TP.StructuredGrid, B, r, θ, ϕ)
   A_field_nu = TP.build_interpolator(TP.StructuredGrid, A, r, θ, ϕ)
   B_field = TP.build_interpolator(TP.StructuredGrid, B, r_uniform, θ, ϕ)
   A_field = TP.build_interpolator(TP.StructuredGrid, A, r_uniform, θ, ϕ)

   return B_field_nu, A_field_nu, B_field, A_field
end

function setup_cartesian_field(ns = 16)
   x = range(-10, 10, length = ns)
   y = range(-10, 10, length = ns)
   z = range(-10, 10, length = ns)
   B = zeros(3, length(x), length(y), length(z)) # vector
   B[3, :, :, :] .= 10e-9
   A = zeros(length(x), length(y), length(z)) # scalar
   A[:, :, :] .= 10e-9

   B_field = TP.build_interpolator(B, x, y, z)
   A_field = TP.build_interpolator(A, x, y, z)

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

   B_field = TP.build_interpolator(TP.RectilinearGrid, B, x, y, z)
   A_field = TP.build_interpolator(TP.RectilinearGrid, A, x, y, z)

   return B_field, A_field
end

function setup_time_dependent_field(ns = 16)
   x = range(-10, 10, length = ns)
   y = range(-10, 10, length = ns)
   z = range(-10, 10, length = ns)

   # Create two time snapshots
   B0 = zeros(3, length(x), length(y), length(z))
   B0[3, :, :, :] .= 1.0 # Bz = 1 at t=0

   B1 = zeros(3, length(x), length(y), length(z))
   B1[3, :, :, :] .= 2.0 # Bz = 2 at t=1

   times = [0.0, 1.0]

   function loader(i)
       if i == 1
           # For demonstration, we assume we load from disk here
           return TP.build_interpolator(TP.CartesianGrid, B0, x, y, z)
       elseif i == 2
           return TP.build_interpolator(TP.CartesianGrid, B1, x, y, z)
       else
           error("Index out of bounds")
       end
   end

   # B_field_t(x, t)
   B_field_t = TP.LazyTimeInterpolator(times, loader)

   return B_field_t
end

function setup_mixed_precision_field(ns = 11, order = 1, bc = TP.FillExtrap(NaN); coeffs = OnTheFly())
   x = range(0.0f0, 10.0f0, length = ns)
   y = range(0.0f0, 10.0f0, length = ns)
   z = range(0.0f0, 10.0f0, length = ns)
   B = fill(0.0f0, 3, ns, ns, ns)
   B[3, :, :, :] .= 1.0f-8

   itp = TP.build_interpolator(B, x, y, z, order, bc; coeffs)
   return itp
end

B_sph_nu, A_sph_nu, B_sph, A_sph = setup_spherical_field();
B_car, A_car = setup_cartesian_field();
B_car_nu, A_car_nu = setup_cartesian_nonuniform_field();
B_td = setup_time_dependent_field();
itp_f32 = setup_mixed_precision_field();

loc = SA[1.0, 1.0, 1.0];
loc_f32 = SA[1.0f0, 1.0f0, 1.0f0];
loc_f64 = SA[1.0, 1.0, 1.0];
```

## Gridded spherical interpolation

!!! note "Input Location"
    For spherical data, the input location is still in Cartesian coordinates!

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

## Mixed precision interpolation

Numerical field data from files is often stored in `Float32`. TestParticle supports constructing interpolators from `Float32` data and ranges, which can then be queried with both `Float32` and `Float64` location vectors.

```@repl interp
@be itp_f32($loc_f32)
@be itp_f32($loc_f64)
```

## Memory usage analysis

Large numerical field arrays can consume significant amounts of memory. It's important that interpolators are memory-efficient during both construction and storage.

For **linear interpolation (`order=1`)**, construction is near zero-allocation because it creates a wrapper around a reinterpreted view of your existing array.

For **higher-order interpolation (`order = 3`)**, `FastInterpolations.jl` must precompute interpolation coefficients to ensure $O(1)$ lookup performance. This process **allocates an additional array of the same size** as your input data.

We can measure this difference using `@be`.

```@repl interp
# Order 1: Minimal allocations (Uses a view)
@be setup_mixed_precision_field(11, 1)

# Order 3: Significant allocations (Computes coefficients)
@be setup_mixed_precision_field(11, 3)
```

Comparing the ratios relative to the original array size illustrates the overhead:

```@repl interp
# Original 4D field array size
B = fill(0.0f0, 3, 11, 11, 11);
size_B = Base.summarysize(B)

# Total size for order=1 (Essentially the input array size)
size_itp1 = Base.summarysize(itp_f32)

# Total size for order=3 (Input array + Coefficient array)
itp_f32_q = setup_mixed_precision_field(11, 3);
size_itp3 = Base.summarysize(itp_f32_q)

# Ratios relative to raw data
size_itp1 / size_B
size_itp3 / size_B
```

As a rule of thumb, linear interpolation has nearly zero memory overhead (ratio ≈ 1.0), while cubic interpolation doubles the memory footprint (ratio ≈ 8.0) to store the coefficients.

## On-the-fly vs Precomputed coefficients

Cubic interpolation (`order = 3`) requires high-order coefficients. By default, TestParticle uses `OnTheFly()` coefficients, which are calculated at query time. This saves memory but increases evaluation time. For maximum performance, you can use `PreCompute()`, which stores the coefficients in an additional array.

```@repl interp
using FastInterpolations: OnTheFly, PreCompute
# Benchmark evaluation time
itp_fly = setup_mixed_precision_field(11, 3; coeffs = OnTheFly());
itp_pre = setup_mixed_precision_field(11, 3; coeffs = PreCompute());

@be itp_fly($loc_f64)
@be itp_pre($loc_f64)

# Compare total object size
Base.summarysize(itp_fly)
Base.summarysize(itp_pre)
```

As shown, `PreCompute()` is roughly 3-4x faster but consumes significantly more memory (comparable to the input data array).

## Time-dependent field interpolation

For time-dependent fields, we can use [`LazyTimeInterpolator`](@ref). It takes a list of time points and a loader function that returns a spatial interpolator for a given time index. The interpolator will linearly interpolate between the two nearest time points.

```@repl interp
@be B_td($loc, 0.5)
```

## Related API

```@docs; canonical=false
TestParticle.build_interpolator
TestParticle.prepare
```
