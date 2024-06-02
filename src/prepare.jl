# Construction of tracing parameters.

"Judge whether the field function is time dependent."
function is_time_dependent(f::Function)
   applicable(f, zeros(6), 0.0) ? true : false
end

abstract type AbstractField{itd} end

"""
    Field{itd, F} <: AbstractField{itd}

A representation of a field function `f`, defined by:

time-independent field
```math
\\mathbf{F} = F(\\mathbf{x}),
```

time-dependent field
```math
\\mathbf{F} = F(\\mathbf{x}, t).
```

# Arguments
- `field_function::Function`: the function of field.
- `itd::Bool`: whether the field function is time dependent.
- `F`: the type of `field_function`.
"""
struct Field{itd, F} <: AbstractField{itd}
   field_function::F
   Field{itd, F}(field_function::F) where {itd, F} =
      isa(itd, Bool) ? new(field_function) : throw(ArgumentError("itd must be a boolean."))
end

Field(f::Function) = Field{is_time_dependent(f), typeof(f)}(f)

(f::AbstractField{true})(xu, t) = f.field_function(xu, t)
(f::AbstractField{true})(xu) = throw(ArgumentError("Time-dependent field function must have a time argument."))
(f::AbstractField{false})(xu, t) = f.field_function(xu)
(f::AbstractField{false})(xu) = f.field_function(xu)

"The type of parameter tuple for full test particle problem."
FullTPTuple = Tuple{Float64, Float64, AbstractField, AbstractField, AbstractField}

"The type of parameter tuple for normal test particle problem."
TPTuple = Tuple{Float64, AbstractField, AbstractField}

"The type of parameter tuple for normalized test particle problem."
TPNormalizedTuple = Tuple{AbstractFloat, AbstractField, AbstractField}

"The type of parameter tuple for guiding center problem."
GCTuple = Tuple{Float64, Float64, Float64, AbstractField, AbstractField}


"""
    prepare(grid::CartesianGrid, E, B; kwargs...) -> (q2m, E, B)

Return a tuple consists of particle charge-mass ratio for a prescribed `species` and
interpolated EM field functions.

# keywords
- `order::Int=1`: order of interpolation in [1,2,3].
- `bc::Int=1`: type of boundary conditions, 1 -> NaN, 2 -> periodic.
- `species::Species=Proton`: particle species.
- `q::AbstractFloat=1.0`: particle charge. Only works when `Species=User`.
- `m::AbstractFloat=1.0`: particle mass. Only works when `Species=User`.

    prepare(grid::CartesianGrid, E, B, F; species=Proton, q=1.0, m=1.0) -> (q, m, E, B, F)

Return a tuple consists of particle charge, mass for a prescribed `species` of charge `q`
and mass `m`, interpolated EM field functions, and external force `F`.

    prepare(x::AbstractRange, y::AbstractRange, z::AbstractRange, E, B; kwargs...) -> (q2m, E, B)
    prepare(x, y, E, B; kwargs...) -> (q2m, E, B)
    prepare(x::AbstractRange, E, B; kwargs...) -> (q2m, E, B)

Direct range input for uniform grid in 2/3D is also accepted.

    prepare(E, B; kwargs...) -> (q2m, E, B)

Return a tuple consists of particle charge-mass ratio for a prescribed `species` of charge
`q` and mass `m` and analytic EM field functions. Prescribed `species` are `Electron` and
`Proton`; other species can be manually specified with `species=Ion/User`, `q` and `m`.

    prepare(E, B, F; kwargs...) -> (q, m, E, B, F)

Return a tuple consists of particle charge, mass for a prescribed `species` of charge `q`
and mass `m`, analytic EM field functions, and external force `F`.
"""
function prepare(grid::CartesianGrid, E::TE, B::TB; species::Species=Proton,
   q::AbstractFloat=1.0, m::AbstractFloat=1.0, order::Int=1, bc::Int=1) where {TE, TB}

   q, m = getchargemass(species, q, m)

   if embeddim(grid) == 3
      gridx, gridy, gridz = makegrid(grid)
      E = TE <: AbstractArray ? getinterp(E, gridx, gridy, gridz, order, bc) : E
      B = TB <: AbstractArray ? getinterp(B, gridx, gridy, gridz, order, bc) : B
   elseif embeddim(grid) == 2
      gridx, gridy = makegrid(grid)
      E = TE <: AbstractArray ? getinterp(E, gridx, gridy, order, bc) : E
      B = TB <: AbstractArray ? getinterp(B, gridx, gridy, order, bc) : B
   end

   q/m, Field(E), Field(B)
end

function prepare(grid::CartesianGrid, E::TE, B::TB, F::TF; species::Species=Proton,
   q::AbstractFloat=1.0, m::AbstractFloat=1.0, order::Int=1, bc::Int=1) where {TE, TB, TF}

   q, m = getchargemass(species, q, m)

   gridx, gridy, gridz = makegrid(grid)

   E = TE <: AbstractArray ? getinterp(E, gridx, gridy, gridz, order, bc) : E
   B = TB <: AbstractArray ? getinterp(B, gridx, gridy, gridz, order, bc) : B
   F = TF <: AbstractArray ? getinterp(F, gridx, gridy, gridz, order, bc) : F

   q, m, Field(E), Field(B), Field(F)
end

function prepare(x::T, y::T, E::TE, B::TB; species::Species=Proton, q::AbstractFloat=1.0,
   m::AbstractFloat=1.0, order::Int=1, bc::Int=1) where {T<:AbstractRange, TE, TB}

   q, m = getchargemass(species, q, m)

   E = TE <: AbstractArray ? getinterp(E, x, y, order, bc) : E
   B = TB <: AbstractArray ? getinterp(B, x, y, order, bc) : B

   q/m, Field(E), Field(B)
end

function prepare(x::T, E::TE, B::TB; species::Species=Proton, q::AbstractFloat=1.0,
   m::AbstractFloat=1.0, order::Int=1, bc::Int=3) where {T<:AbstractRange, TE, TB}

   q, m = getchargemass(species, q, m)

   E = TE <: AbstractArray ? getinterp(E, x, order, bc) : E
   B = TB <: AbstractArray ? getinterp(B, x, order, bc) : B

   q/m, Field(E), Field(B)
end

function prepare(x::T, y::T, z::T, E::TE, B::TB;
   species::Species=Proton, q::AbstractFloat=1.0, m::AbstractFloat=1.0, order::Int=1,
   bc::Int=1) where {T<:AbstractRange, TE, TB}

   q, m = getchargemass(species, q, m)

   E = TE <: AbstractArray ? getinterp(E, x, y, z, order, bc) : E
   B = TB <: AbstractArray ? getinterp(B, x, y, z, order, bc) : B

   q/m, Field(E), Field(B)
end

function prepare(E, B; species::Species=Proton, q::AbstractFloat=1.0, m::AbstractFloat=1.0)
   q, m = getchargemass(species, q, m)

   q/m, Field(E), Field(B)
end

function prepare(E, B, F; species::Species=Proton, q::AbstractFloat=1.0,
   m::AbstractFloat=1.0)
   q, m = getchargemass(species, q, m)

   q, m, Field(E), Field(B), Field(F)
end

function prepare_gc(xv, xrange::T, yrange::T, zrange::T, E::TE, B::TB;
   species::Species=Proton, q::AbstractFloat=1.0, m::AbstractFloat=1.0, order::Int=1,
   bc::Int=1, removeExB=true) where {T<:AbstractRange, TE, TB}

   q, m = getchargemass(species, q, m)
   x, v = @views xv[1:3], xv[4:6]

   E = TE <: AbstractArray ? getinterp(E, xrange, yrange, zrange, order, bc) : E
   B = TB <: AbstractArray ? getinterp(B, xrange, yrange, zrange, order, bc) : B

   bparticle = B(x)
   Bmag_particle = √(bparticle[1]^2 + bparticle[2]^2 + bparticle[3]^2)
   b̂particle = bparticle ./ Bmag_particle
   # vector of Larmor radius
   ρ = (b̂particle × v) ./ (q/m*Bmag_particle)
   # Get the guiding center location
   X = x - ρ
   # Get EM field at guiding center
   b = B(X)
   Bmag = √(b[1]^2 + b[2]^2 + b[3]^2)
   b̂ = b ./ Bmag
   vpar = @views b̂ ⋅ v

   vperp = @. v - vpar * b̂
   if removeExB
      e = E(X)
      vE = e × b̂ / Bmag
      w = vperp - vE
   else
      w = vperp
   end
   μ = m * (w ⋅ w) / (2 * Bmag)

   stateinit_gc = [X..., vpar]

   stateinit_gc, (q, m, μ, Field(E), Field(B))
end

function prepare_gc(xv, E, B; species::Species=Proton, q::AbstractFloat=1.0,
   m::AbstractFloat=1.0, removeExB=true)
   q, m = getchargemass(species, q, m)
   x, v = @views xv[1:3], xv[4:6]

   bparticle = B(x)
   Bmag_particle = √(bparticle[1]^2 + bparticle[2]^2 + bparticle[3]^2)
   b̂particle = bparticle ./ Bmag_particle
   # vector of Larmor radius
   ρ = (b̂particle × v) ./ (q/m*Bmag_particle)
   # Get the guiding center location
   X = x - ρ
   # Get EM field at guiding center
   b = B(X)
   Bmag = √(b[1]^2 + b[2]^2 + b[3]^2)
   b̂ = b ./ Bmag
   vpar = @views b̂ ⋅ v

   vperp = @. v - vpar * b̂
   if removeExB
      e = E(X)
      vE = e × b̂ / Bmag
      w = vperp - vE
   else
      w = vperp
   end
   μ = m * (w ⋅ w) / (2 * Bmag)

   stateinit_gc = [X..., vpar]

   stateinit_gc, (q, m, μ, Field(E), Field(B))
end

"""
    guiding_center(xu, param::Union{TPTuple, FullTPTuple})

Calculate the coordinates of the guiding center according to the phase space coordinates of a particle.
Reference: https://en.wikipedia.org/wiki/Guiding_center

A simple definition:
```math
\\mathbf{X}=\\mathbf{x}-m\\frac{\\mathbf{v}\\times\\mathbf{B}}{qB}
```
"""
function guiding_center(xu, param::TPTuple)
   q2m, _, B_field = param
   t = xu[end]
   v = @view xu[4:6]
   Bv = B_field(xu, t)
   B = sqrt(Bv[1]^2 + Bv[2]^2 + Bv[3]^2)
   # unit vector along B
   b = Bv ./ B
   # vector of Larmor radius
   ρ = (b × v) ./ (q2m*B)

   X = @views xu[1:3] - ρ
end

function guiding_center(xu, param::FullTPTuple)
   q, m, _, B_field = param
   t = xu[end]
   v = @view xu[4:6]
   Bv = B_field(xu, t)
   B = sqrt(Bv[1]^2 + Bv[2]^2 + Bv[3]^2)
   # unit vector along B
   b̂ = Bv ./ B
   # vector of Larmor radius
   ρ = (b̂ × v) ./ (q*B/m)

   X = @views xu[1:3] - ρ
end

function guiding_center(xu, q, m, B_field)
   q2m = q / m
   t = xu[end]
   v = @view xu[4:6]
   Bv = B_field(xu, t)
   B = sqrt(Bv[1]^2 + Bv[2]^2 + Bv[3]^2)
   # unit vector along B
   b̂ = Bv ./ B
   # the vector of Larmor radius
   ρ = (b̂ × v) ./ (q2m*B)

   X = @views xu[1:3] - ρ
end

"""
    get_gc(param::Union{TPTuple, FullTPTuple})

Get three functions for plotting the orbit of guiding center.

For example:
```julia
param = prepare(E, B; species=Proton)
gc = get_gc(param)
# The definitions of stateinit, tspan, E and B are ignored.
prob = ODEProblem(trace!, stateinit, tspan, param)
sol = solve(prob, Vern7(); dt=2e-11)

f = Figure(fontsize=18)
ax = Axis3(f[1, 1], aspect = :data)
gc_plot(x,y,z,vx,vy,vz) = (gc([x,y,z,vx,vy,vz])...,)
lines!(ax, sol, idxs=(gc_plot, 1, 2, 3, 4, 5, 6))
```
"""
function get_gc(param::Union{TPTuple, FullTPTuple})
   gc(xu) = guiding_center(xu, param)
end