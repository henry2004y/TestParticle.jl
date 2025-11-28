
# Internal {#Internal}

## Public APIs {#Public-APIs}
<details class='jldocstring custom-block' open>
<summary><a id='TestParticle.BiMaxwellian' href='#TestParticle.BiMaxwellian'><span class="jlbinding">TestParticle.BiMaxwellian</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



Type for BiMaxwellian velocity distributions with respect to the magnetic field.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/henry2004y/TestParticle.jl/blob/4c31acb4391a45c5eba4a7610d42ce329e4172dc/src/sampler.jl#L31-L33" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TestParticle.BiMaxwellian-Union{Tuple{U}, Tuple{T}, Tuple{AbstractVector{U}, AbstractVector{T}, Any, Any, Any}} where {T<:AbstractFloat, U<:AbstractFloat}' href='#TestParticle.BiMaxwellian-Union{Tuple{U}, Tuple{T}, Tuple{AbstractVector{U}, AbstractVector{T}, Any, Any, Any}} where {T<:AbstractFloat, U<:AbstractFloat}'><span class="jlbinding">TestParticle.BiMaxwellian</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
 BiMaxwellian(B::Vector{U}, u0::Vector{T}, ppar, pperp, n; m=mᵢ)
```


Construct a BiMaxwellian distribution with magnetic field `B`, bulk velocity `u0`, parallel thermal pressure `ppar`, perpendicular thermal pressure `pperp`, and number density `n` in SI units. The default particle is proton.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/henry2004y/TestParticle.jl/blob/4c31acb4391a45c5eba4a7610d42ce329e4172dc/src/sampler.jl#L45-L51" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TestParticle.Cartesian' href='#TestParticle.Cartesian'><span class="jlbinding">TestParticle.Cartesian</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



Cartesian grid.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/henry2004y/TestParticle.jl/blob/4c31acb4391a45c5eba4a7610d42ce329e4172dc/src/utility/interpolation.jl#L7-L9" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TestParticle.Maxwellian' href='#TestParticle.Maxwellian'><span class="jlbinding">TestParticle.Maxwellian</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



Type for Maxwellian velocity distributions.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/henry2004y/TestParticle.jl/blob/4c31acb4391a45c5eba4a7610d42ce329e4172dc/src/sampler.jl#L8-L10" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TestParticle.Maxwellian-Union{Tuple{T}, Tuple{AbstractVector{T}, Any, Any}} where T' href='#TestParticle.Maxwellian-Union{Tuple{T}, Tuple{AbstractVector{T}, Any, Any}} where T'><span class="jlbinding">TestParticle.Maxwellian</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
 Maxwellian(u0::AbstractVector{T}, p, n; m=mᵢ)
```


Construct a Maxwellian distribution with bulk velocity `u0`, thermal pressure `p`, and number density `n` in SI units. The default particle is proton.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/henry2004y/TestParticle.jl/blob/4c31acb4391a45c5eba4a7610d42ce329e4172dc/src/sampler.jl#L18-L23" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TestParticle.Spherical' href='#TestParticle.Spherical'><span class="jlbinding">TestParticle.Spherical</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



Spherical grid with uniform r, θ and ϕ.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/henry2004y/TestParticle.jl/blob/4c31acb4391a45c5eba4a7610d42ce329e4172dc/src/utility/interpolation.jl#L11-L13" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TestParticle.SphericalNonUniformR' href='#TestParticle.SphericalNonUniformR'><span class="jlbinding">TestParticle.SphericalNonUniformR</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



Spherical grid with non-uniform r and uniform θ, ϕ.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/henry2004y/TestParticle.jl/blob/4c31acb4391a45c5eba4a7610d42ce329e4172dc/src/utility/interpolation.jl#L15-L17" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TestParticle.energy2velocity-Tuple{Any}' href='#TestParticle.energy2velocity-Tuple{Any}'><span class="jlbinding">TestParticle.energy2velocity</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Return velocity magnitude from energy in [eV].


<Badge type="info" class="source-link" text="source"><a href="https://github.com/henry2004y/TestParticle.jl/blob/4c31acb4391a45c5eba4a7610d42ce329e4172dc/src/utility/utility.jl#L211-L213" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TestParticle.get_energy-Tuple{Any}' href='#TestParticle.get_energy-Tuple{Any}'><span class="jlbinding">TestParticle.get_energy</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Calculate the energy [eV] of a relativistic particle from γv.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/henry2004y/TestParticle.jl/blob/4c31acb4391a45c5eba4a7610d42ce329e4172dc/src/utility/utility.jl#L200-L202" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TestParticle.get_energy-Tuple{SciMLBase.AbstractODESolution}' href='#TestParticle.get_energy-Tuple{SciMLBase.AbstractODESolution}'><span class="jlbinding">TestParticle.get_energy</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Return the energy [eV] from relativistic `sol`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/henry2004y/TestParticle.jl/blob/4c31acb4391a45c5eba4a7610d42ce329e4172dc/src/utility/utility.jl#L184-L186" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TestParticle.get_gc-Tuple{Any, Any}' href='#TestParticle.get_gc-Tuple{Any, Any}'><span class="jlbinding">TestParticle.get_gc</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
 get_gc(xu, param)
 get_gc(x, y, z, vx, vy, vz, bx, by, bz, q2m)
```


Calculate the coordinates of the guiding center according to the phase space coordinates of a particle. Reference: [wiki](https://en.wikipedia.org/wiki/Guiding_center)

Nonrelativistic definition:

$$\mathbf{X}=\mathbf{x}-m\frac{\mathbf{b}\times\mathbf{v}}{qB}$$


<Badge type="info" class="source-link" text="source"><a href="https://github.com/henry2004y/TestParticle.jl/blob/4c31acb4391a45c5eba4a7610d42ce329e4172dc/src/gc.jl#L76-L88" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TestParticle.get_gc_func-Tuple{Any}' href='#TestParticle.get_gc_func-Tuple{Any}'><span class="jlbinding">TestParticle.get_gc_func</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
 get_gc_func(param)
```


Return the function for plotting the orbit of guiding center.

**Example**

```julia
param = prepare(E, B; species = Proton)
# The definitions of stateinit, tspan, E and B are ignored.
prob = ODEProblem(trace!, stateinit, tspan, param)
sol = solve(prob, Vern7(); dt = 2e-11)

f = Figure(fontsize = 18)
ax = Axis3(f[1, 1], aspect = :data)
gc = param |> get_gc_func
gc_plot(x, y, z, vx, vy, vz) = (gc(SA[x, y, z, vx, vy, vz])...,)
lines!(ax, sol, idxs = (gc_plot, 1, 2, 3, 4, 5, 6))
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/henry2004y/TestParticle.jl/blob/4c31acb4391a45c5eba4a7610d42ce329e4172dc/src/gc.jl#L136-L155" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TestParticle.get_gyrofrequency' href='#TestParticle.get_gyrofrequency'><span class="jlbinding">TestParticle.get_gyrofrequency</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



Return the gyrofrequency.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/henry2004y/TestParticle.jl/blob/4c31acb4391a45c5eba4a7610d42ce329e4172dc/src/utility/utility.jl#L130-L132" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TestParticle.get_gyroperiod' href='#TestParticle.get_gyroperiod'><span class="jlbinding">TestParticle.get_gyroperiod</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



Return the gyroperiod.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/henry2004y/TestParticle.jl/blob/4c31acb4391a45c5eba4a7610d42ce329e4172dc/src/utility/utility.jl#L154-L156" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TestParticle.get_gyroradius-Tuple{AbstractFloat, AbstractFloat}' href='#TestParticle.get_gyroradius-Tuple{AbstractFloat, AbstractFloat}'><span class="jlbinding">TestParticle.get_gyroradius</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Return the gyroradius.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/henry2004y/TestParticle.jl/blob/4c31acb4391a45c5eba4a7610d42ce329e4172dc/src/utility/utility.jl#L141-L143" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TestParticle.get_velocity-Tuple{Any}' href='#TestParticle.get_velocity-Tuple{Any}'><span class="jlbinding">TestParticle.get_velocity</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Return velocity from relativistic γv in `sol`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/henry2004y/TestParticle.jl/blob/4c31acb4391a45c5eba4a7610d42ce329e4172dc/src/utility/utility.jl#L166-L168" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TestParticle.prepare' href='#TestParticle.prepare'><span class="jlbinding">TestParticle.prepare</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
prepare(args...; kwargs...) -> (q2m, m, E, B, F)
prepare(E, B, F = ZeroField(); kwargs...)
prepare(grid::CartesianGrid, E, B, F = ZeroField(); kwargs...)
prepare(x, E, B, F = ZeroField(); dir = 1, kwargs...)
prepare(x, y, E, B, F = ZeroField(); kwargs...)
prepare(x, y, z, E, B, F = ZeroField(); kwargs...)
```


Return a tuple consists of particle charge-mass ratio for a prescribed `species` of charge `q` and mass `m`, mass `m` for a prescribed `species`, analytic/interpolated EM field functions, and external force `F`.

Prescribed `species` are `Electron` and `Proton`; other species can be manually specified with `species=Ion/User`, `q` and `m`.

Direct range input for uniform grid in 1/2/3D is supported. For 1D grid, an additional keyword `dir` is used for specifying the spatial direction, 1 -&gt; x, 2 -&gt; y, 3 -&gt; z. For 3D grid, the default grid type is `Cartesian`. To use `Spherical` grid, an additional keyword `gridtype` is needed. For `Spherical` grid, dimensions of field arrays should be `(Br, Bθ, Bϕ)`.

**Keywords**
- `order::Int=1`: order of interpolation in [1,2,3].
  
- `bc::Int=1`: type of boundary conditions, 1 -&gt; NaN, 2 -&gt; periodic.
  
- `species::Species=Proton`: particle species.
  
- `q=1.0`: particle charge. Only works when `Species=User`.
  
- `m=1.0`: particle mass. Only works when `Species=User`.
  
- `gridtype::Grid=Cartesian()`: type of grid in `Cartesian()`, `Spherical()`, `SphericalNonUniformR`.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/henry2004y/TestParticle.jl/blob/4c31acb4391a45c5eba4a7610d42ce329e4172dc/src/prepare.jl#L81-L108" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TestParticle.sample-Union{Tuple{Maxwellian{U, T}}, Tuple{T}, Tuple{U}} where {U, T}' href='#TestParticle.sample-Union{Tuple{Maxwellian{U, T}}, Tuple{T}, Tuple{U}} where {U, T}'><span class="jlbinding">TestParticle.sample</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
 sample(vdf::Maxwellian)
```


Sample a 3D velocity from a [`Maxwellian`](/api#TestParticle.Maxwellian) distribution `vdf` using the Box-Muller method.

```julia
 sample(vdf::BiMaxwellian)
```


Sample a 3D velocity from a [`BiMaxwellian`](/api#TestParticle.BiMaxwellian) distribution `vdf` using the Box-Muller method.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/henry2004y/TestParticle.jl/blob/4c31acb4391a45c5eba4a7610d42ce329e4172dc/src/sampler.jl#L69-L77" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TestParticle.trace!-NTuple{4, Any}' href='#TestParticle.trace!-NTuple{4, Any}'><span class="jlbinding">TestParticle.trace!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
trace!(dy, y, p, t)
```


ODE equations for charged particle moving in static EM field and external force field with in-place form.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/henry2004y/TestParticle.jl/blob/4c31acb4391a45c5eba4a7610d42ce329e4172dc/src/equations.jl#L25-L29" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TestParticle.trace-Tuple{Any, Any, Any}' href='#TestParticle.trace-Tuple{Any, Any, Any}'><span class="jlbinding">TestParticle.trace</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
trace(y, p, t) -> SVector{6, Float64}
```


ODE equations for charged particle moving in static EM field and external force field with out-of-place form.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/henry2004y/TestParticle.jl/blob/4c31acb4391a45c5eba4a7610d42ce329e4172dc/src/equations.jl#L38-L42" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TestParticle.trace_gc!-Tuple{Any, Any, Tuple{Float64, Float64, Float64, TestParticle.AbstractField, TestParticle.AbstractField}, Any}' href='#TestParticle.trace_gc!-Tuple{Any, Any, Tuple{Float64, Float64, Float64, TestParticle.AbstractField, TestParticle.AbstractField}, Any}'><span class="jlbinding">TestParticle.trace_gc!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
 trace_gc!(dy, y, p, t)
```


Guiding center equations for nonrelativistic charged particle moving in static EM field with in-place form. Variable `y = (x, y, z, u)`, where `u` is the velocity along the magnetic field at (x,y,z).


<Badge type="info" class="source-link" text="source"><a href="https://github.com/henry2004y/TestParticle.jl/blob/4c31acb4391a45c5eba4a7610d42ce329e4172dc/src/equations.jl#L154-L159" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TestParticle.trace_gc_1st!-Tuple{Any, Any, Tuple{Float64, Float64, Float64, TestParticle.AbstractField, TestParticle.AbstractField}, Any}' href='#TestParticle.trace_gc_1st!-Tuple{Any, Any, Tuple{Float64, Float64, Float64, TestParticle.AbstractField, TestParticle.AbstractField}, Any}'><span class="jlbinding">TestParticle.trace_gc_1st!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



1st order approximation of guiding center equations.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/henry2004y/TestParticle.jl/blob/4c31acb4391a45c5eba4a7610d42ce329e4172dc/src/equations.jl#L196-L198" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TestParticle.trace_gc_drifts!-NTuple{4, Any}' href='#TestParticle.trace_gc_drifts!-NTuple{4, Any}'><span class="jlbinding">TestParticle.trace_gc_drifts!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
 trace_gc_drifts!(dx, x, p, t)
```


Equations for tracing the guiding center using analytical drifts, including the grad-B drift, curvature drift, and ExB drift. Parallel velocity is also added. This expression requires the full particle trajectory `p.sol`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/henry2004y/TestParticle.jl/blob/4c31acb4391a45c5eba4a7610d42ce329e4172dc/src/equations.jl#L127-L132" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TestParticle.trace_normalized!-NTuple{4, Any}' href='#TestParticle.trace_normalized!-NTuple{4, Any}'><span class="jlbinding">TestParticle.trace_normalized!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
 trace_normalized!(dy, y, p, t)
```


Normalized ODE equations for charged particle moving in static EM field with in-place form. If the field is in 2D X-Y plane, periodic boundary should be applied for the field in z via the extrapolation function provided by Interpolations.jl.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/henry2004y/TestParticle.jl/blob/4c31acb4391a45c5eba4a7610d42ce329e4172dc/src/equations.jl#L76-L82" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TestParticle.trace_relativistic!-NTuple{4, Any}' href='#TestParticle.trace_relativistic!-NTuple{4, Any}'><span class="jlbinding">TestParticle.trace_relativistic!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
 trace_relativistic!(dy, y, p, t)
```


ODE equations for relativistic charged particle (x, γv) moving in static EM field with in-place form.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/henry2004y/TestParticle.jl/blob/4c31acb4391a45c5eba4a7610d42ce329e4172dc/src/equations.jl#L49-L53" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TestParticle.trace_relativistic-Tuple{Any, Any, Any}' href='#TestParticle.trace_relativistic-Tuple{Any, Any, Any}'><span class="jlbinding">TestParticle.trace_relativistic</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
 trace_relativistic(y, p, t) -> SVector{6}
```


ODE equations for relativistic charged particle (x, γv) moving in static EM field with out-of-place form.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/henry2004y/TestParticle.jl/blob/4c31acb4391a45c5eba4a7610d42ce329e4172dc/src/equations.jl#L63-L67" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TestParticle.trace_relativistic_normalized!-NTuple{4, Any}' href='#TestParticle.trace_relativistic_normalized!-NTuple{4, Any}'><span class="jlbinding">TestParticle.trace_relativistic_normalized!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
 trace_relativistic_normalized!(dy, y, p, t)
```


Normalized ODE equations for relativistic charged particle (x, γv) moving in static EM field with in-place form.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/henry2004y/TestParticle.jl/blob/4c31acb4391a45c5eba4a7610d42ce329e4172dc/src/equations.jl#L94-L98" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TestParticle.trace_relativistic_normalized-Tuple{Any, Any, Any}' href='#TestParticle.trace_relativistic_normalized-Tuple{Any, Any, Any}'><span class="jlbinding">TestParticle.trace_relativistic_normalized</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
 trace_relativistic_normalized(y, p, t)
```


Normalized ODE equations for relativistic charged particle (x, γv) moving in static EM field with out-of-place form.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/henry2004y/TestParticle.jl/blob/4c31acb4391a45c5eba4a7610d42ce329e4172dc/src/equations.jl#L111-L115" target="_blank" rel="noreferrer">source</a></Badge>

</details>


## Private types and methods {#Private-types-and-methods}
<details class='jldocstring custom-block' open>
<summary><a id='TestParticle.AbstractTraceSolution' href='#TestParticle.AbstractTraceSolution'><span class="jlbinding">TestParticle.AbstractTraceSolution</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



Abstract type for tracing solutions.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/henry2004y/TestParticle.jl/blob/4c31acb4391a45c5eba4a7610d42ce329e4172dc/src/TestParticle.jl#L36-L38" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TestParticle.Field' href='#TestParticle.Field'><span class="jlbinding">TestParticle.Field</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



```julia
 Field{itd, F} <: AbstractField{itd}
```


A representation of a field function `f`, defined by:

time-independent field

$$\mathbf{F} = F(\mathbf{x}),$$

time-dependent field

$$\mathbf{F} = F(\mathbf{x}, t).$$

**Arguments**
- `field_function::Function`: the function of field.
  
- `itd::Bool`: whether the field function is time dependent.
  
- `F`: the type of `field_function`.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/henry2004y/TestParticle.jl/blob/4c31acb4391a45c5eba4a7610d42ce329e4172dc/src/prepare.jl#L12-L34" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TestParticle.GCTuple' href='#TestParticle.GCTuple'><span class="jlbinding">TestParticle.GCTuple</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



The type of parameter tuple for guiding center problem.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/henry2004y/TestParticle.jl/blob/4c31acb4391a45c5eba4a7610d42ce329e4172dc/src/prepare.jl#L56-L58" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TestParticle.Grid' href='#TestParticle.Grid'><span class="jlbinding">TestParticle.Grid</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



Type for grid.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/henry2004y/TestParticle.jl/blob/4c31acb4391a45c5eba4a7610d42ce329e4172dc/src/utility/interpolation.jl#L3-L5" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TestParticle.Species' href='#TestParticle.Species'><span class="jlbinding">TestParticle.Species</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



Type for the particles, `Proton`, `Electron`, `Ion`, or `User`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/henry2004y/TestParticle.jl/blob/4c31acb4391a45c5eba4a7610d42ce329e4172dc/src/TestParticle.jl#L31-L33" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TestParticle.TraceSolution-Union{Tuple{Any}, Tuple{deriv}, Tuple{Any, Type{deriv}}} where deriv' href='#TestParticle.TraceSolution-Union{Tuple{Any}, Tuple{deriv}, Tuple{Any, Type{deriv}}} where deriv'><span class="jlbinding">TestParticle.TraceSolution</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Interpolate solution at time `x`. Forward tracing only.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/henry2004y/TestParticle.jl/blob/4c31acb4391a45c5eba4a7610d42ce329e4172dc/src/pusher.jl#L52-L54" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TestParticle.VDF' href='#TestParticle.VDF'><span class="jlbinding">TestParticle.VDF</span></a> <Badge type="info" class="jlObjectType jlType" text="Type" /></summary>



Abstract type for velocity distribution functions.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/henry2004y/TestParticle.jl/blob/4c31acb4391a45c5eba4a7610d42ce329e4172dc/src/sampler.jl#L3-L5" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TestParticle._boris!-NTuple{8, Any}' href='#TestParticle._boris!-NTuple{8, Any}'><span class="jlbinding">TestParticle._boris!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Apply Boris method for particles with index in `irange`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/henry2004y/TestParticle.jl/blob/4c31acb4391a45c5eba4a7610d42ce329e4172dc/src/pusher.jl#L228-L230" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TestParticle._prepare-Tuple{TraceProblem, Any, Any, Any}' href='#TestParticle._prepare-Tuple{TraceProblem, Any, Any, Any}'><span class="jlbinding">TestParticle._prepare</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Prepare for advancing.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/henry2004y/TestParticle.jl/blob/4c31acb4391a45c5eba4a7610d42ce329e4172dc/src/pusher.jl#L216-L218" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TestParticle.cart2sph-Tuple{Any, Any, Any}' href='#TestParticle.cart2sph-Tuple{Any, Any, Any}'><span class="jlbinding">TestParticle.cart2sph</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Convert from Cartesian to spherical coordinates vector.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/henry2004y/TestParticle.jl/blob/4c31acb4391a45c5eba4a7610d42ce329e4172dc/src/utility/utility.jl#L14-L16" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TestParticle.cross!-Tuple{Any, Any, Any}' href='#TestParticle.cross!-Tuple{Any, Any, Any}'><span class="jlbinding">TestParticle.cross!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



In-place cross product.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/henry2004y/TestParticle.jl/blob/4c31acb4391a45c5eba4a7610d42ce329e4172dc/src/pusher.jl#L167-L169" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TestParticle.dipole-Tuple{Any, Any}' href='#TestParticle.dipole-Tuple{Any, Any}'><span class="jlbinding">TestParticle.dipole</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Calculates the magnetic field from a dipole with magnetic moment `M` at `r`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/henry2004y/TestParticle.jl/blob/4c31acb4391a45c5eba4a7610d42ce329e4172dc/src/utility/dipole.jl#L16-L18" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TestParticle.dipole_fieldline' href='#TestParticle.dipole_fieldline'><span class="jlbinding">TestParticle.dipole_fieldline</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
 dipole_fieldline(ϕ, L=2.5, nP=100)
```


Creates `nP` points on one field line of the magnetic field from a dipole. In a centered dipole magnetic field model, the path along a given L shell can be described as r = L*cos²λ, where r is the radial distance (in planetary radii) to a point on the line, λ is its co-latitude, and L is the L-shell of interest.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/henry2004y/TestParticle.jl/blob/4c31acb4391a45c5eba4a7610d42ce329e4172dc/src/utility/dipole.jl#L29-L36" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TestParticle.getB_CS_harris' href='#TestParticle.getB_CS_harris'><span class="jlbinding">TestParticle.getB_CS_harris</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
 getB_CS_harris(B₀, L)
```


Return the magnetic field at location `r` near a current sheet with magnetic strength `B₀` and sheet length `L`. The current sheet is assumed to lie in the z = 0 plane.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/henry2004y/TestParticle.jl/blob/4c31acb4391a45c5eba4a7610d42ce329e4172dc/src/utility/current_sheet.jl#L3-L8" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TestParticle.getB_bottle-NTuple{8, Any}' href='#TestParticle.getB_bottle-NTuple{8, Any}'><span class="jlbinding">TestParticle.getB_bottle</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
 getB_bottle(x, y, z, distance, a, b, I1, I2) -> StaticVector{Float64, 3}
```


Get magnetic field from a magnetic bottle. Reference: [wiki](https://en.wikipedia.org/wiki/Magnetic_mirror#Magnetic_bottles)

**Arguments**
- `x,y,z::Float`: particle coordinates in [m].
  
- `distance::Float`: distance between solenoids in [m].
  
- `a::Float`: radius of each side coil in [m].
  
- `b::Float`: radius of central coil in [m].
  
- `I1::Float`: current in the solenoid times number of windings in side coils in [A].
  
- `I2::Float`: current in the central solenoid times number of windings in the central loop in [A].
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/henry2004y/TestParticle.jl/blob/4c31acb4391a45c5eba4a7610d42ce329e4172dc/src/utility/confinement.jl#L70-L85" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TestParticle.getB_current_loop-Tuple{Any, Any, Any, TestParticle.Currentloop}' href='#TestParticle.getB_current_loop-Tuple{Any, Any, Any, TestParticle.Currentloop}'><span class="jlbinding">TestParticle.getB_current_loop</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
 getB_current_loop(x, y, z, cl::Currentloop) -> StaticVector{Float64, 3}
```


Get magnetic field at `[x, y, z]` from a magnetic mirror generated from two coils.

**Arguments**
- `x,y,z::Float`: particle coordinates in [m].
  
- `distance::Float`: distance between solenoids in [m].
  
- `a::Float`: radius of each side coil in [m].
  
- `I1::Float`: current in the solenoid times number of windings in side coils.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/henry2004y/TestParticle.jl/blob/4c31acb4391a45c5eba4a7610d42ce329e4172dc/src/utility/confinement.jl#L18-L29" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TestParticle.getB_dipole-Tuple{Any}' href='#TestParticle.getB_dipole-Tuple{Any}'><span class="jlbinding">TestParticle.getB_dipole</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Analytic magnetic field function for testing. Return in SI unit.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/henry2004y/TestParticle.jl/blob/4c31acb4391a45c5eba4a7610d42ce329e4172dc/src/utility/dipole.jl#L8-L10" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TestParticle.getB_mirror-NTuple{6, Any}' href='#TestParticle.getB_mirror-NTuple{6, Any}'><span class="jlbinding">TestParticle.getB_mirror</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
 getB_mirror(x, y, z, distance, a, I1) -> StaticVector{Float64, 3}
```


Get magnetic field at `[x, y, z]` from a magnetic mirror generated from two coils.

**Arguments**
- `x,y,z::Float`: particle coordinates in [m].
  
- `distance::Float`: distance between solenoids in [m].
  
- `a::Float`: radius of each side coil in [m].
  
- `I1::Float`: current in the solenoid times number of windings in side coils.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/henry2004y/TestParticle.jl/blob/4c31acb4391a45c5eba4a7610d42ce329e4172dc/src/utility/confinement.jl#L43-L54" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TestParticle.getB_tokamak_coil-NTuple{7, Any}' href='#TestParticle.getB_tokamak_coil-NTuple{7, Any}'><span class="jlbinding">TestParticle.getB_tokamak_coil</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
 getB_tokamak_coil(x, y, z, a, b, ICoils, IPlasma) -> StaticVector{Float64, 3}
```


Get the magnetic field from a Tokamak topology consists of 16 coils. Original: [Tokamak-Fusion-Reactor](https://github.com/BoschSamuel/Simulation-of-a-Tokamak-Fusion-Reactor/blob/master/Simulation2.m)

**Arguments**
- `x,y,z::Float`: location in [m].
  
- `a::Float`: radius of each coil in [m].
  
- `b::Float`: radius of central region in [m].
  
- `ICoil::Float`: current in the coil times number of windings in [A].
  
- `IPlasma::Float`: current of the plasma in [A].
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/henry2004y/TestParticle.jl/blob/4c31acb4391a45c5eba4a7610d42ce329e4172dc/src/utility/confinement.jl#L101-L114" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TestParticle.getB_tokamak_profile-Tuple{AbstractFloat, AbstractFloat, AbstractFloat, Any, AbstractFloat, AbstractFloat, AbstractFloat}' href='#TestParticle.getB_tokamak_profile-Tuple{AbstractFloat, AbstractFloat, AbstractFloat, Any, AbstractFloat, AbstractFloat, AbstractFloat}'><span class="jlbinding">TestParticle.getB_tokamak_profile</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
 getB_tokamak_profile(x, y, z, q_profile, a, R₀, Bζ0) -> StaticVector{Float64, 3}
```


Reconstruct the magnetic field distribution from a safe factor(q) profile. Reference: Tokamak, 4th Edition, John Wesson.

**Arguments**
- `x,y,z::Float`: location in [m].
  
- `q_profile::Function`: profile of q. The variable of this function must be the normalized radius.
  
- `a::Float`: minor radius [m].
  
- `R₀::Float`: major radius [m].
  
- `Bζ0::Float`: toroidal magnetic field on axis [T].
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/henry2004y/TestParticle.jl/blob/4c31acb4391a45c5eba4a7610d42ce329e4172dc/src/utility/confinement.jl#L186-L199" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TestParticle.getE_dipole-Tuple{Any}' href='#TestParticle.getE_dipole-Tuple{Any}'><span class="jlbinding">TestParticle.getE_dipole</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Analytic electric field function for testing.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/henry2004y/TestParticle.jl/blob/4c31acb4391a45c5eba4a7610d42ce329e4172dc/src/utility/dipole.jl#L3-L5" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TestParticle.get_interpolator-Union{Tuple{T}, Tuple{Cartesian, AbstractArray{T, 4}, Any, Any, Any}, Tuple{Cartesian, AbstractArray{T, 4}, Any, Any, Any, Int64}, Tuple{Cartesian, AbstractArray{T, 4}, Any, Any, Any, Int64, Int64}} where T' href='#TestParticle.get_interpolator-Union{Tuple{T}, Tuple{Cartesian, AbstractArray{T, 4}, Any, Any, Any}, Tuple{Cartesian, AbstractArray{T, 4}, Any, Any, Any, Int64}, Tuple{Cartesian, AbstractArray{T, 4}, Any, Any, Any, Int64, Int64}} where T'><span class="jlbinding">TestParticle.get_interpolator</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
 get_interpolator(A, gridx, gridy, gridz, order::Int=1, bc::Int=1)
 get_interpolator(gridtype, A, grid1, grid2, grid3, order::Int=1, bc::Int=1)
```


Return a function for interpolating field array `A` on the grid given by `gridx`, `gridy`, and `gridz`.

**Arguments**
- `gridtype`: `Cartesian`, `Spherical` or `SphericalNonUniformR`.
  
- `A`: field array. For vector field, the first dimension should be 3.
  
- `order::Int=1`: order of interpolation in [1,2,3].
  
- `bc::Int=1`: type of boundary conditions, 1 -&gt; NaN, 2 -&gt; periodic, 3 -&gt; Flat.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/henry2004y/TestParticle.jl/blob/4c31acb4391a45c5eba4a7610d42ce329e4172dc/src/utility/interpolation.jl#L163-L176" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TestParticle.get_rotation_matrix-Tuple{AbstractVector{<:AbstractFloat}, Real}' href='#TestParticle.get_rotation_matrix-Tuple{AbstractVector{<:AbstractFloat}, Real}'><span class="jlbinding">TestParticle.get_rotation_matrix</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
 get_rotation_matrix(axis::AbstractVector, angle::Real) --> SMatrix{3,3}
```


Create a rotation matrix for rotating a 3D vector around a unit `axis` by an `angle` in radians. Reference: [Rotation matrix from axis and angle](https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle)

**Example**

```julia
using LinearAlgebra
v = [-0.5, 1.0, 1.0]
v̂ = normalize(v)
θ = deg2rad(-74)
R = get_rotation_matrix(v̂, θ)
```



<Badge type="info" class="source-link" text="source"><a href="https://github.com/henry2004y/TestParticle.jl/blob/4c31acb4391a45c5eba4a7610d42ce329e4172dc/src/utility/utility.jl#L105-L121" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TestParticle.getchargemass-Tuple{TestParticle.Species, Real, Real}' href='#TestParticle.getchargemass-Tuple{TestParticle.Species, Real, Real}'><span class="jlbinding">TestParticle.getchargemass</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
 getchargemass(species::Species, q, m)
```


Return charge and mass for `species`. For `species = Ion`, `q` and `m` are charge and mass numbers. For `species = User`, the input `q` and `m` are returned as is.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/henry2004y/TestParticle.jl/blob/4c31acb4391a45c5eba4a7610d42ce329e4172dc/src/utility/utility.jl#L46-L51" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TestParticle.getinterp' href='#TestParticle.getinterp'><span class="jlbinding">TestParticle.getinterp</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
 getinterp(::Grid, A, gridx, gridy, gridz, order::Int=1, bc::Int=1)
```


Return a function for interpolating field array `A` on the grid given by `gridx`, `gridy`, and `gridz`.

**Arguments**
- `order::Int=1`: order of interpolation in [1,2,3].
  
- `bc::Int=1`: type of boundary conditions, 1 -&gt; NaN, 2 -&gt; periodic, 3 -&gt; Flat.
  
- `dir::Int`: 1/2/3, representing x/y/z direction.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/henry2004y/TestParticle.jl/blob/4c31acb4391a45c5eba4a7610d42ce329e4172dc/src/utility/interpolation.jl#L30-L41" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TestParticle.getinterp_scalar' href='#TestParticle.getinterp_scalar'><span class="jlbinding">TestParticle.getinterp_scalar</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
 getinterp_scalar(::Grid, A, gridx, gridy, gridz, order::Int=1, bc::Int=1)
```


Return a function for interpolating scalar array `A` on the grid given by `gridx`, `gridy`, and `gridz`. Currently only 3D arrays are supported.

**Arguments**
- `order::Int=1`: order of interpolation in [1,2,3].
  
- `bc::Int=1`: type of boundary conditions, 1 -&gt; NaN, 2 -&gt; periodic, 3 -&gt; Flat.
  
- `dir::Int`: 1/2/3, representing x/y/z direction.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/henry2004y/TestParticle.jl/blob/4c31acb4391a45c5eba4a7610d42ce329e4172dc/src/utility/interpolation.jl#L141-L152" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TestParticle.is_time_dependent-Tuple{Function}' href='#TestParticle.is_time_dependent-Tuple{Function}'><span class="jlbinding">TestParticle.is_time_dependent</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Judge whether the field function is time dependent.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/henry2004y/TestParticle.jl/blob/4c31acb4391a45c5eba4a7610d42ce329e4172dc/src/prepare.jl#L7-L9" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TestParticle.makegrid-Tuple{Meshes.CartesianGrid{M, C, N} where {M<:Meshes.𝔼, C<:CoordRefSystems.Cartesian, N}}' href='#TestParticle.makegrid-Tuple{Meshes.CartesianGrid{M, C, N} where {M<:Meshes.𝔼, C<:CoordRefSystems.Cartesian, N}}'><span class="jlbinding">TestParticle.makegrid</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Return uniform range from 2D/3D CartesianGrid.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/henry2004y/TestParticle.jl/blob/4c31acb4391a45c5eba4a7610d42ce329e4172dc/src/utility/utility.jl#L66-L68" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TestParticle.set_axes_equal-Tuple{Any}' href='#TestParticle.set_axes_equal-Tuple{Any}'><span class="jlbinding">TestParticle.set_axes_equal</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
 set_axes_equal(ax)
```


Set 3D plot axes to equal scale for Matplotlib. Make axes of 3D plot have equal scale so that spheres appear as spheres and cubes as cubes. Required since `ax.axis('equal')` and `ax.set_aspect('equal')` don&#39;t work on 3D.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/henry2004y/TestParticle.jl/blob/4c31acb4391a45c5eba4a7610d42ce329e4172dc/src/utility/utility.jl#L85-L91" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TestParticle.solve' href='#TestParticle.solve'><span class="jlbinding">TestParticle.solve</span></a> <Badge type="info" class="jlObjectType jlFunction" text="Function" /></summary>



```julia
 solve(prob::TraceProblem; trajectories::Int=1, dt::AbstractFloat,
	 savestepinterval::Int=1, isoutofdomain::Function=ODE_DEFAULT_ISOUTOFDOMAIN)
```


Trace particles using the Boris method with specified `prob`.

**keywords**
- `trajectories::Int`: number of trajectories to trace.
  
- `dt::AbstractFloat`: time step.
  
- `savestepinterval::Int`: saving output interval.
  
- `isoutofdomain::Function`: a function with input of position and velocity vector `xv` that determines whether to stop tracing.
  


<Badge type="info" class="source-link" text="source"><a href="https://github.com/henry2004y/TestParticle.jl/blob/4c31acb4391a45c5eba4a7610d42ce329e4172dc/src/pusher.jl#L178-L190" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TestParticle.sph2cart-Tuple{Any, Any, Any}' href='#TestParticle.sph2cart-Tuple{Any, Any, Any}'><span class="jlbinding">TestParticle.sph2cart</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Convert from spherical to Cartesian coordinates vector.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/henry2004y/TestParticle.jl/blob/4c31acb4391a45c5eba4a7610d42ce329e4172dc/src/utility/utility.jl#L3-L5" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TestParticle.sph_to_cart_vector-NTuple{5, Any}' href='#TestParticle.sph_to_cart_vector-NTuple{5, Any}'><span class="jlbinding">TestParticle.sph_to_cart_vector</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Convert a vector from spherical to Cartesian.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/henry2004y/TestParticle.jl/blob/4c31acb4391a45c5eba4a7610d42ce329e4172dc/src/utility/utility.jl#L29-L31" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TestParticle.update_location!-Tuple{Any, Any}' href='#TestParticle.update_location!-Tuple{Any, Any}'><span class="jlbinding">TestParticle.update_location!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



Update location in one timestep `dt`.


<Badge type="info" class="source-link" text="source"><a href="https://github.com/henry2004y/TestParticle.jl/blob/4c31acb4391a45c5eba4a7610d42ce329e4172dc/src/pusher.jl#L156-L158" target="_blank" rel="noreferrer">source</a></Badge>

</details>

<details class='jldocstring custom-block' open>
<summary><a id='TestParticle.update_velocity!-NTuple{5, Any}' href='#TestParticle.update_velocity!-NTuple{5, Any}'><span class="jlbinding">TestParticle.update_velocity!</span></a> <Badge type="info" class="jlObjectType jlMethod" text="Method" /></summary>



```julia
 update_velocity!(xv, paramBoris, param, dt, t)
```


Update velocity using the Boris method, Birdsall, Plasma Physics via Computer Simulation. Reference: [DTIC](https://apps.dtic.mil/sti/citations/ADA023511)


<Badge type="info" class="source-link" text="source"><a href="https://github.com/henry2004y/TestParticle.jl/blob/4c31acb4391a45c5eba4a7610d42ce329e4172dc/src/pusher.jl#L114-L119" target="_blank" rel="noreferrer">source</a></Badge>

</details>

