# # Shock Drift Acceleration at a Perpendicular Shock
#
# This example shows Shock Drift Acceleration (SDA): a particle gains energy by drifting
# along the shock front in the convective electric field while it repeatedly hits the shock.
#
# The setup is the textbook *perpendicular* shock (``\theta_{Bn}=90^\circ``):
#
# | Quantity            | Direction      |
# |---------------------|----------------|
# | Shock normal ``\hat{\mathbf{n}}`` | ``+\hat{\mathbf{x}}`` (upstream ``x<0``, downstream ``x>0``) |
# | Upstream flow ``\mathbf{V}_1``    | ``+\hat{\mathbf{x}}`` (plasma flows into the shock) |
# | Magnetic field ``\mathbf{B}_1``   | ``+\hat{\mathbf{z}}`` (purely tangential) |
# | Convection field ``\mathbf{E}=-\mathbf{V}_1\times\mathbf{B}_1`` | ``+\hat{\mathbf{y}}`` |
#
# The field only varies across the shock, so the only drift along the front is the
# gradient-``\mathbf{B}`` drift
# ``\mathbf{v}_{\nabla B}\propto(\mathbf{B}\times\nabla B)/B^3``. Since ``\nabla B`` points
# downstream (``+\hat{\mathbf{x}}``) and ``\mathbf{B}\parallel\hat{\mathbf{z}}``, the drift
# runs along ``+\hat{\mathbf{y}}`` and is *parallel* to ``\mathbf{E}``. Its power input
# ``P=q\,\mathbf{v}_{\nabla B}\!\cdot\!\mathbf{E}`` is therefore positive, so the proton is
# accelerated as it slides along the front — the essence of SDA. (The ``\mathbf{E}\times
# \mathbf{B}`` drift just carries the particle into the shock along ``\hat{\mathbf{x}}``; it
# does not by itself accelerate the particle.)
#
# ## Cross-shock potential
#
# A perpendicular shock also carries an electrostatic potential ``\phi(x)`` (downstream at
# higher potential). Its normal field ``\mathbf{E}_x=-\nabla\phi`` points upstream and
# reflects the proton: after one hit it turns back, the ExB drift brings it in again, and it
# meets the front repeatedly. Each crossing lets the tangential drift add more energy. (Real
# shocks reflect only some ions; here the potential is enhanced so the reflection is easy to
# see.)
#
# We trace a real *proton* (``m = mᵢ``, ``q = qᵢ``). The shock has a finite thickness (a
# smooth tanh ramp of order the proton gyroradius) so the gradient driving the drift is
# resolved.

import DisplayAs #hide
using TestParticle
using TestParticle: mᵢ, qᵢ
using OrdinaryDiffEq
using LinearAlgebra
using CairoMakie
CairoMakie.activate!(type = "png") #hide

# ## Shock Parameters (perpendicular)
#
# Upstream (Region 1) flows toward the shock along ``+x``; downstream (Region 2) is
# compressed and slower by the ratio ``r``. ``\mathbf{B}`` is along ``z`` (tangential), so
# the motional field ``\mathbf{E}=-\mathbf{V}\times\mathbf{B}`` points along ``y`` and stays
# continuous across the shock (de Hoffmann-Teller frame).

V₁ = [545.0, 0.0, 0.0] .* 1.0e3   # upstream flow (into the shock, +x)
r = 4.0                            # field compression ratio
B₁ = [0.0, 0.0, 5.0] .* 1.0e-9      # upstream |B| [T], along z
B₂ = [0.0, 0.0, r * 5.0] .* 1.0e-9  # downstream, same direction
V₂ = V₁ ./ r                        # downstream flow (mass-flux conservation)
E_y = -(V₁ × B₁)[2]                 # motional E along +y, continuous across shock
println("E_y = ", E_y, "  (points along +y, as expected)")

# Cross-shock potential (downstream at higher potential) and its normal field:
#   φ(x) = Δφ · ½(1 + tanh(x/L))      E_x = -dφ/dx = -Δφ/(2L)·sech²(x/L)  (points upstream)
Δφ = 15.0e3                          # cross-shock potential drop [V] (enhanced for visibility)
L = 1500.0e3                        # shock ramp half-width
Ex(xq) = -Δφ / (2L) / cosh(xq / L)^2

# ## Grid and finite-thickness shock
#
# The fields live on a 1D grid in ``x``; the shock sits at ``x = 0`` and is smoothed over a
# width ``L`` so the proton (gyroradius ``\sim 1100`` km) feels a resolved gradient rather
# than a discontinuity.

x_max = 12_000.0e3
x = range(-x_max, x_max; length = 4000)
Bvecs = [B₁ .+ (B₂ .- B₁) .* 0.5 .* (1 .+ tanh(xq / L)) for xq in x]
Bm = hcat(Bvecs...)
Em = hcat([[Ex(xq), E_y, 0.0] for xq in x]...)

param = prepare(x, Em, Bm; bc = ClampExtrap(), q = qᵢ, m = mᵢ)

function field_at(xq)
    i = clamp(searchsortedlast(x, xq) + 1, 1, length(x))
    return Em[:, i], Bm[:, i]
end

function drift_velocity(xq)
    Evec, Bvec = field_at(xq)
    return (Evec × Bvec) / dot(Bvec, Bvec)
end

# ## Proton Injection
#
# The proton starts upstream and gyrates in the ``x``-``y`` plane (velocity along ``-y``),
# so it is not merely convected. Carried into the shock by the normal ExB drift, it feels
# the tangential gradient-B drift and is accelerated along ``+y``; the cross-shock potential
# then reflects it, so it meets the front again and again.

const x₀ = -2500.0e3
const v₀ = [0.0, -545.0e3, 0.0]
const t_max = 35.0

u0 = [x₀, 0.0, 0.0, v₀...]
prob = ODEProblem(trace!, u0, (0.0, t_max), param)
sol = solve(prob, Vern9(); saveat = 0.05)

# ## Energy: the SDA signature
#
# The particle's kinetic energy ``K = \tfrac12 m |\mathbf{v}|^2`` oscillates within each
# gyration. The meaningful SDA quantity is the *work done by the electric field*,
# ``W(t)=\int_0^t q\,\mathbf{v}\cdot\mathbf{E}\,dt``, which isolates the energy fed in by the
# tangential drift (parallel to ``\mathbf{E}``). It rises with every encounter, confirming
# the repeated acceleration. We also count the shock crossings to show how often it hits the
# front.

ts = sol.t
X = [u[1] for u in sol.u]
Y = [u[2] for u in sol.u]
Z = [u[3] for u in sol.u]
vs = [u[4:6] for u in sol.u]

K_lab = [0.5 * mᵢ * dot(v, v) for v in vs]
W_cum = cumsum([qᵢ * dot(v, field_at(u[1])[1]) for (u, v) in zip(sol.u, vs)]) .* sol.t[2]

n_cross = count(i -> (X[i - 1] < 0 <= X[i] || X[i - 1] > 0 >= X[i]), 2:length(X))

K_norm = K_lab ./ K_lab[1]
W_norm = W_cum ./ K_lab[1]

println("Shock encounters (crossings of x=0): ", n_cross)
println(
    "Perpendicular-shock SDA: proton KE ends at ", round(K_norm[end]; digits = 2),
    "x; cumulative work from E ends at ", round(W_norm[end]; digits = 2), "x initial"
) #hide

# ## Visualization

f = Figure(size = (1500, 650), fontsize = 18)

xpad = 500.0e3; ypad = 2000.0e3; zpad = 500.0e3
xlim = (minimum(X) - xpad, maximum(X) + xpad)
ylim = (minimum(Y) - ypad, maximum(Y) + ypad)
zlim = (minimum(Z) - zpad, maximum(Z) + zpad)

ax3d = Axis3(
    f[1, 1];
    title = "Proton SDA trajectory (colored by gained energy)",
    xlabel = "x [km]", ylabel = "y [km]", zlabel = "z [km]",
    aspect = :data,
    limits = (xlim ./ 1.0e3, ylim ./ 1.0e3, zlim ./ 1.0e3),
)

lines!(
    ax3d, X ./ 1.0e3, Y ./ 1.0e3, Z ./ 1.0e3;
    color = K_norm, colormap = :plasma, linewidth = 2.0
)
cb = Colorbar(
    f[1, 2], limits = (1, maximum(K_norm)),
    colormap = :plasma, label = "KE / KE₀"
)
cb.width = 18

# Shock front (semi-transparent plane at x = 0)
mesh!(
    ax3d,
    Rect3f(
        Point3f(-10.0e3, ylim[1], zlim[1]),
        Vec3f(20.0e3, ylim[2] - ylim[1], zlim[2] - zlim[1])
    );
    color = (:red, 0.2),
)

axE = Axis(
    f[1, 3];
    title = "Energy gain",
    xlabel = "Time [s]", ylabel = "Energy / initial",
    limits = (nothing, (0, nothing)),
)
lines!(axE, ts, K_norm; label = "instantaneous KE (lab frame)", color = :gray, linewidth = 1.5)
lines!(axE, ts, W_norm; label = "cumulative work q∫v·E dt (SDA)", color = :black, linewidth = 2.5)
axislegend(axE; position = :lt, framevisible = true, backgroundcolor = (:white, 0.7))

f = DisplayAs.PNG(f) #hide

# The proton is carried into the red shock front by the normal ExB drift, reflected by the
# cross-shock potential, and meets the front repeatedly. Each encounter adds energy from
# ``\mathbf{E}`` (the black curve — the cumulative work — climbs with every crossing), so
# the proton surfs along ``+y`` and leaves with a large, sustained gain: shock drift
# acceleration enabled by the cross-shock potential.
