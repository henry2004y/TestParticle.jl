# Field Interpolation

A robust field interpolation is the prerequisite for pushing particles.
This example demonstrates the construction of scalar/vector field interpolators for Cartesian/Spherical grids.
If the field is analytic, you can directly pass the generated function to [`prepare`](@ref).

```@example interp
using TestParticle
using Meshes
using StaticArrays
using Chairmarks
import TestParticle as TP

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

   B_field_nu = build_interpolator(StructuredGrid, B, r, θ, ϕ)
   A_field_nu = build_interpolator(StructuredGrid, A, r, θ, ϕ)
   B_field = build_interpolator(StructuredGrid, B, r_uniform, θ, ϕ)
   A_field = build_interpolator(StructuredGrid, A, r_uniform, θ, ϕ)

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

   B_field = build_interpolator(B, x, y, z)
   A_field = build_interpolator(A, x, y, z)

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

   B_field = build_interpolator(RectilinearGrid, B, x, y, z)
   A_field = build_interpolator(RectilinearGrid, A, x, y, z)

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
           return build_interpolator(CartesianGrid, B0, x, y, z)
       elseif i == 2
           return build_interpolator(CartesianGrid, B1, x, y, z)
       else
           error("Index out of bounds")
       end
   end

   # B_field_t(x, t)
   B_field_t = LazyTimeInterpolator(times, loader)

   return B_field_t
end

function setup_mixed_precision_field(
      ns = 11, order = 1, bc = FillExtrap(NaN);
      coeffs = OnTheFly(), store = StorePolicy()
   )
   x = range(0.0f0, 10.0f0, length = ns)
   y = range(0.0f0, 10.0f0, length = ns)
   z = range(0.0f0, 10.0f0, length = ns)
   B = fill(0.0f0, 3, ns, ns, ns)
   B[3, :, :, :] .= 1.0f-8

   itp = build_interpolator(B, x, y, z, order, bc; coeffs, store)
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

Large numerical field arrays can consume significant amounts of memory during
interpolator construction. `FastInterpolations.jl` uses `StorePolicy()` by default,
which copies grids and field data so that the interpolator owns a stable snapshot.
This protects it from later caller mutations, but constructing an interpolator for a
45 GB field requires another approximately 45 GB allocation.

TestParticle forwards the `store` keyword to `FastInterpolations.jl`. Use
`StorePolicy(copy = false)` to alias compatible grids and field data instead of
copying them. This removes the field-sized construction allocation for linear
interpolation and for cubic interpolation backed by dense scalar or `SVector` arrays.

```@repl interp
x_memory = range(0.0f0, 1.0f0, length = 32);
B_memory = rand(Float32, 3, 32, 32, 32);
B_memory_svector = fill(SA[1.0f0, 1.0f0, 1.0f0], 32, 32, 32);
copy_store = StorePolicy();
zero_copy_store = StorePolicy(copy = false);

# Default owned-copy construction
@be build_interpolator($B_memory, $x_memory, $x_memory, $x_memory, 1; store = $copy_store)
@be build_interpolator($B_memory_svector, $x_memory, $x_memory, $x_memory, 3; store = $copy_store)

# Zero-copy construction
@be build_interpolator($B_memory, $x_memory, $x_memory, $x_memory, 1; store = $zero_copy_store)
@be build_interpolator($B_memory_svector, $x_memory, $x_memory, $x_memory, 3; store = $zero_copy_store)
```

With zero-copy storage, the interpolator may alias both the field and grid arrays.
They must remain alive and must not be mutated or resized while the interpolator is
in use. If only the large field should be aliased, use
`StorePolicy(copy_values = false)` to retain owned grid storage. Inputs that require
an element-type conversion may still be copied.

Component-first vector fields with shape `(3, nx, ny, nz)` are exposed internally as
a `ReinterpretArray`. FastInterpolations v0.4 currently cannot evaluate ND cardinal
interpolation from that aliased wrapper, so TestParticle warns and falls back to
owned storage for `order = 3`. Store the field directly as a dense 3D array of
`SVector`s to use zero-copy cubic interpolation.

`Base.summarysize` reports all objects reachable through an interpolator, including
aliased inputs, so construction allocation measurements are a better way to verify
that the copy was eliminated.

## On-the-fly vs Precomputed coefficients

Cubic interpolation (`order = 3`) requires high-order coefficients. TestParticle uses
`OnTheFly()` coefficients, which are calculated during each query. Zero-copy storage
reduces construction time and memory but does not change interpolation throughput.
`PreCompute()` could trade additional memory for faster queries, but it is not yet
supported for TestParticle's ND cardinal interpolation.

!!! note "Status in FastInterpolations v0.4.19"
    `PreCompute()` is available for the global natural cubic spline (`CubicInterp`) in ND. However, it is **not yet supported** for the local Hermite cubic spline (`CardinalInterp`, i.e. `order = 3` in TestParticle) in ND. Constructing such an interpolator with `coeffs = PreCompute()` raises an `ArgumentError`, so `OnTheFly()` remains the only option for cubic interpolation in TestParticle.

```@repl interp
# Benchmark evaluation time
itp_fly = setup_mixed_precision_field(11, 3; coeffs = OnTheFly());
# itp_pre = setup_mixed_precision_field(11, 3; coeffs = PreCompute()); # unsupported for ND cardinal cubic

@be itp_fly($loc_f64)

# Compare total object size
Base.summarysize(itp_fly)
```

As shown, `OnTheFly()` preserves memory efficiency while providing higher-order accuracy. Once `PreCompute()` support for ND local Hermite cubic interpolation lands, it will offer a faster alternative for memory-abundant systems.

## Time-dependent field interpolation

For time-dependent fields, we can use [`LazyTimeInterpolator`](@ref). It takes a list of time points and a loader function that returns a spatial interpolator for a given time index. The interpolator will linearly interpolate between the two nearest time points.

```@repl interp
@be B_td($loc, 0.5)
```

## GPU interpolation

The interpolators above are built on FastInterpolations.jl and live on the host. Pushing a
particle on a device asks for something different: a plain struct with `Int32` indices whose
axes and data can be copied to device memory and evaluated inside a kernel. TestParticle
provides `GPUSphericalGrid` for spherical grids, and `GPUGrid3D`, `GPUGrid2D` and `GPUGrid1D`
for Cartesian ones.

`prepare` performs the conversion by itself once the parameters are handed to a device
backend, see [GPU Ensemble Tracing](@ref). The very same grids can be built on the `CPU()`
backend, which is what the following does with `TP.CPU()`, so that the device interpolator can
be inspected and benchmarked where no device is available:

```@example interp
B_gpu = TP._to_gpu_spherical_grid(B_sph.itp, TP.CPU());
A_gpu = TP._to_gpu_spherical_grid(A_sph.itp, TP.CPU());
B_nu_gpu = TP._to_gpu_spherical_grid(B_sph_nu.itp, TP.CPU());
```

A uniform grid vector becomes a `GPUUniformAxis`, which locates a cell with a multiply, and a
non-uniform one a `GPUNonUniformAxis`, which brackets it with a binary search:

```@repl interp
typeof(B_gpu)
B_gpu.axis_r
typeof(B_nu_gpu.axis_r)
```

Queries are unchanged, only the backend behind them differs. The location is still Cartesian,
and a vector result is rotated back into the Cartesian basis:

```@repl interp
@be B_gpu($loc)
@be B_sph($loc)
B_gpu($loc)
A_gpu($loc)
```

Outside the grid the boundary condition decides, filling with `NaN` in `r` and in `θ` and
wrapping periodically in `ϕ`:

```@repl interp
B_gpu(0.01, 0.0, 0.0)
B_gpu(20.0, 0.0, 0.0)
```

One caveat comes with the device grids: they always interpolate linearly, whatever `order`
the host interpolator was built with. A field prepared with `order = 3` is therefore not
reproduced on the device, and `order = 1`, the default, is the one to use when the host and
the device are meant to agree.

## Related API

```@docs; canonical=false
build_interpolator
prepare
```
