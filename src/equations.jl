# Tracing equations.
get_x(u) = @inbounds SA[u[1], u[2], u[3]]
get_v(u) = @inbounds SA[u[4], u[5], u[6]]

"""
	trace!(dy, y, p, t)

ODE equations for charged particle moving in static EM field and external force field with in-place form.
"""
function trace!(dy, y, p, t)
	q2m, m, Efunc, Bfunc, Ffunc = p

	v = get_v(y)
	E = Efunc(y, t)
	B = Bfunc(y, t)
	F = Ffunc(y, t)

	dy[1:3] = v
	dy[4:6] = SVector{3}(q2m * (v × B + E) + F / m)

	return
end

"""
	trace(y, p, t) -> SVector{6, Float64}

ODE equations for charged particle moving in static EM field and external force field with out-of-place form.
"""
function trace(y, p, t)
	q2m, m, Efunc, Bfunc, Ffunc = p

	v = get_v(y)
	E = Efunc(y, t)
	B = Bfunc(y, t)
	F = Ffunc(y, t)

	dv = SVector{3}(q2m * (v × B + E) + F / m)

	vcat(v, dv)
end

"""
	 trace_relativistic!(dy, y, p, t)

ODE equations for relativistic charged particle (x, γv) moving in static EM field with in-place form.
"""
function trace_relativistic!(dy, y, p, t)
	q2m, m, Efunc, Bfunc, Ffunc = p
	E = Efunc(y, t)
	B = Bfunc(y, t)
	F = Ffunc(y, t)

	γv = get_v(y)
	γ²v² = γv[1]^2 + γv[2]^2 + γv[3]^2
	if γ²v² > eps(eltype(dy))
		v̂ = normalize(γv)
	else # no velocity
		v̂ = SVector{3, eltype(dy)}(0, 0, 0)
	end
	vmag = √(γ²v² / (1 + γ²v² / c^2))
	v = vmag * v̂

	dy[1:3] = v
	dy[4:6] = SVector{3}(q2m * (v × B + E) + F / m)

	return
end

"""
	 trace_relativistic(y, p, t) -> SVector{6}

ODE equations for relativistic charged particle (x, γv) moving in static EM field with out-of-place form.
"""
function trace_relativistic(y, p, t)
	q2m, m, Efunc, Bfunc, Ffunc = p
	E = Efunc(y, t)
	B = Bfunc(y, t)
	F = Ffunc(y, t)

	γv = get_v(y)
	γ²v² = γv[1]^2 + γv[2]^2 + γv[3]^2
	if γ²v² > eps(eltype(y))
		v̂ = normalize(γv)
	else # no velocity
		v̂ = SVector{3, eltype(y)}(0, 0, 0)
	end
	vmag = √(γ²v² / (1 + γ²v²/c^2))
	v = vmag * v̂
	dv = SVector{3}(q2m * (v × B + E))

	vcat(v, dv)
end

"""
	 trace_normalized!(dy, y, p, t)

Normalized ODE equations for charged particle moving in static EM field with in-place form.
If the field is in 2D X-Y plane, periodic boundary should be applied for the field in z via
the extrapolation function provided by Interpolations.jl.
"""
function trace_normalized!(dy, y, p, t)
	_, m, Efunc, Bfunc = p

	v = get_v(y)
	E = Efunc(y, t)
	B = Bfunc(y, t)

	dy[1:3] = v
	dy[4:6] = SVector{3}(v × B + E)

	return
end

"""
	 trace_relativistic_normalized!(dy, y, p, t)

Normalized ODE equations for relativistic charged particle (x, γv) moving in static EM field with in-place form.
"""
function trace_relativistic_normalized!(dy, y, p, t)
	_, m, Efunc, Bfunc = p
	E = Efunc(y, t)
	B = Bfunc(y, t)
	γv = get_v(y)

	γ²v² = γv[1]^2 + γv[2]^2 + γv[3]^2
	if γ²v² > eps(eltype(dy))
		v̂ = normalize(γv)
	else # no velocity
		v̂ = SVector{3,eltype(dy)}(0, 0, 0)
	end
	vmag = √(γ²v² / (1 + γ²v²))
	v = vmag * v̂

	dy[1:3] = v
	dy[4:6] = SVector{3}(v × B + E)

	return
end

"""
	 trace_relativistic_normalized(y, p, t)

Normalized ODE equations for relativistic charged particle (x, γv) moving in static EM field with out-of-place form.
"""
function trace_relativistic_normalized(y, p, t)
	_, m, Efunc, Bfunc = p
	E = Efunc(y, t)
	B = Bfunc(y, t)
	γv = get_v(y)

	γ²v² = γv[1]^2 + γv[2]^2 + γv[3]^2
	if γ²v² > eps(eltype(y))
		v̂ = normalize(γv)
	else # no velocity
		v̂ = SVector{3,eltype(y)}(0, 0, 0)
	end
	vmag = √(γ²v² / (1 + γ²v²))
	v = vmag * v̂
	dv = SVector{3}(v × B + E)

	vcat(v, dv)
end

"""
	 trace_gc_drifts!(dx, x, p, t)

Equations for tracing the guiding center using analytical drifts, including the grad-B drift, curvature drift, and ExB drift.
Parallel velocity is also added. This expression requires the full particle trajectory `p.sol`.
"""
function trace_gc_drifts!(dx, x, p, t)
	q2m, m, Efunc, Bfunc, sol = p
	xu = sol(t)
	v = get_v(xu)
	E = Efunc(x)
	B = Bfunc(x)

	Bmag(x) = √(Bfunc(x) ⋅ Bfunc(x))
	∇B = ForwardDiff.gradient(Bmag, x)
	b = normalize(B)
	v_par = (v ⋅ b) .* b
	v_perp = v - v_par
	Ω = q2m * norm(B)
	κ = ForwardDiff.jacobian(Bfunc, x) * B  # B⋅∇B
	## v⟂^2*(B×∇|B|)/(2*Ω*B^2) + v∥^2*(B×(B⋅∇B))/(Ω*B^3) + (E×B)/B^2 + v∥
	dx[1:3] =
		norm(v_perp)^2 * (B × ∇B) / (2 * Ω * norm(B)^2) +
		norm(v_par)^2 * (B × κ) / Ω / norm(B)^3 + (E × B) / (B ⋅ B) + v_par
end

"""
	 trace_gc!(dy, y, p::TPTuple, t)

Guiding center equations for nonrelativistic charged particle moving in static EM field with in-place form.
Variable `y = (x, y, z, u)`, where `u` is the velocity along the magnetic field at (x,y,z).
"""
function trace_gc!(dy, y, p::GCTuple, t)
	q, m, μ, Efunc, Bfunc = p
	q2m = q / m
	X = get_x(y)
	E = Efunc(X, t)
	B = Bfunc(X, t)
	b̂ = normalize(B) # unit B field at X

	Bmag(x) = √(Bfunc(x) ⋅ Bfunc(x))
	∇B = SVector{3}(ForwardDiff.gradient(Bmag, X))

	function E2B(x)
		E² = Efunc(x) ⋅ Efunc(x)
		B² = Bfunc(x) ⋅ Bfunc(x)
		vE² = E² / B²
	end

	∇vE² = SVector{3}(ForwardDiff.gradient(E2B, X))

	bfunc(x) = normalize(Bfunc(x, t))
	∇b̂ = ForwardDiff.jacobian(bfunc, X)
	# ∇ × b̂
	curlb = SVector{3}(∇b̂[3, 2] - ∇b̂[2, 3], ∇b̂[1, 3] - ∇b̂[3, 1], ∇b̂[2, 1] - ∇b̂[1, 2])
	# effective EM fields
	Eᵉ = @. E - (μ * ∇B - 0.5 * m * ∇vE²) / q
	Bᵉ = @. B + y[4] / q2m * curlb
	Bparᵉ = b̂ ⋅ Bᵉ # effective B field parallel to B

	dy[1] = (y[4] * Bᵉ[1] + Eᵉ[2] * b̂[3] - Eᵉ[3] * b̂[2]) / Bparᵉ
	dy[2] = (y[4] * Bᵉ[2] + Eᵉ[3] * b̂[1] - Eᵉ[1] * b̂[3]) / Bparᵉ
	dy[3] = (y[4] * Bᵉ[3] + Eᵉ[1] * b̂[2] - Eᵉ[2] * b̂[1]) / Bparᵉ
	dy[4] = q2m / Bparᵉ * Bᵉ ⋅ Eᵉ

	return
end

"1st order approximation of guiding center equations."
function trace_gc_1st!(dy, y, p::GCTuple, t)
	q, m, μ, Efunc, Bfunc = p
	q2m = q / m
	X = get_x(y)
	E = Efunc(X, t)
	B = Bfunc(X, t)
	b̂ = normalize(B) # unit B field at X
	u = y[4]

	Bmag(x) = √(Bfunc(x) ⋅ Bfunc(x))
	∇B = SVector{3}(ForwardDiff.gradient(Bmag, X))
	Ω = q * Bmag(X) / m

	bfunc(x) = normalize(Bfunc(x, t))
	∇b̂ = ForwardDiff.jacobian(bfunc, X)
	# effective EM fields
	Eᵉ = @. E - (μ * ∇B) / q
	κ = ∇b̂ * b̂  # curvature
	vX = u * b̂ + b̂ × (κ * u^2 - q2m * Eᵉ) / Ω

	dy[1] = vX[1]
	dy[2] = vX[2]
	dy[3] = vX[3]
	dy[4] = q2m * (b̂ + u * (b̂ × κ) / Ω) ⋅ Eᵉ

	return
end
