# Guiding center.

function prepare_gc(xv, xrange::T, yrange::T, zrange::T, E::TE, B::TB;
	species::Species = Proton, q::AbstractFloat = 1.0, m::AbstractFloat = 1.0,
	order::Int = 1,
	bc::Int = 1, removeExB = true) where {T <: AbstractRange, TE, TB}

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

function prepare_gc(xv, E, B; species::Species = Proton, q::AbstractFloat = 1.0,
	m::AbstractFloat = 1.0, removeExB = true)
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
	 get_gc(xu, param::Union{TPTuple, FullTPTuple})
	 get_gc(x, y, z, vx, vy, vz, bx, by, bz, q2m)

Calculate the coordinates of the guiding center according to the phase space coordinates of a particle.
Reference: [wiki](https://en.wikipedia.org/wiki/Guiding_center)

Nonrelativistic definition:
```math
\\mathbf{X}=\\mathbf{x}-m\\frac{\\mathbf{b}\\times\\mathbf{v}}{qB}
```
"""
function get_gc(xu, param::TPTuple)
	q2m, _, B_field = param
	t = xu[end]
	v = @view xu[4:6]
	B = B_field(xu, t)
	B² = B[1]^2 + B[2]^2 + B[3]^2
	# vector of Larmor radius
	ρ = B × v ./ (q2m*B²)

	X = @views xu[1:3] - ρ
end

function get_gc(xu, param::FullTPTuple)
	q, m, _, B_field = param
	t = xu[end]
	v = @view xu[4:6]
	B = B_field(xu, t)
	B² = B[1]^2 + B[2]^2 + B[3]^2
	# vector of Larmor radius
	ρ = B × v ./ (q*B²/m)

	X = @views xu[1:3] - ρ
end

function get_gc(x, y, z, vx, vy, vz, bx, by, bz, q2m)
	l = SVector{3}(x, y, z)
	B = SVector{3}(bx, by, bz)
	v = SVector{3}(vx, vy, vz)

	B² = bx^2 + by^2 + bz^2
	# vector of Larmor radius
	ρ = B × v ./ (q2m*B²)

	X = l - ρ
end

function get_gc(
	x::T,
	y::T,
	z::T,
	vx::T,
	vy::T,
	vz::T,
	bx::U,
	by::U,
	bz::U,
	q2m,
) where {T <: AbstractVector, U <: AbstractVector}
	X = [zeros(SVector{3, eltype(x)}) for _ in x]
	for i in eachindex(X)
		X[i] =
			get_gc(x[i], y[i], z[i], vx[i], vy[i], vz[i], bx[i], by[i], bz[i], q2m)
	end

	X
end

get_gc(x, y, z, vx, vy, vz, bx, by, bz, q, m) =
	get_gc(x, y, z, vx, vy, vz, bx, by, bz, q/m)

"""
	 get_gc(param::Union{TPTuple, FullTPTuple})

Return the function for plotting the orbit of guiding center.

# Example
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
	gc(xu) = get_gc(xu, param)
end