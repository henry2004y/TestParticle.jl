# # Extra Diagnostics
#
# This example demonstrates tracing one proton in an analytic E field and numerical B field.
#
# The `SavingCallback` from DiffEqCallbacks.jl can be used to save additional outputs for diagnosis. Here we save the magnetic field along the trajectory, together with the parallel velocity.
# Note that `SavingCallback` is currently not compatible with ensemble problems; for multiple particle tracing with customized outputs, see [Demo: ensemble tracing with extra saving](@ref Ensemble-Tracing).
#
# For the fields and the work done by them, `get_fields` and `get_work_rates` compute the values along any saved trajectory without callbacks.

using TestParticle, OrdinaryDiffEq, StaticArrays
using BorisPushers
using TestParticle: WrapExtrap
using TestParticle: qᵢ, mᵢ
using Statistics
using LinearAlgebra: normalize, ×, ⋅
using DiffEqCallbacks


## Number of cells for the field along each dimension
nx, ny, nz = 4, 6, 8
## Spatial coordinates given in customized units
x = range(-0.5, 0.5, length = nx)
y = range(-0.5, 0.5, length = ny)
z = range(-0.5, 0.5, length = nz)
## Numerical magnetic field given in customized units
B = Array{Float32, 4}(undef, 3, nx, ny, nz)

B[1, :, :, :] .= 0.0
B[2, :, :, :] .= 0.0
B[3, :, :, :] .= 2.0

## Reference values for unit conversions between the customized and dimensionless units
const B₀ = sqrt(sum(abs2, B) / (nx * ny * nz)) # RMS magnitude
const U₀ = 1.0
const l₀ = 4 * nx
const t₀ = l₀ / U₀
const E₀ = U₀ * B₀

### Convert from customized to default dimensionless units
## Dimensionless spatial extents [l₀]
x /= l₀
y /= l₀
z /= l₀
## For full EM problems, the normalization of E and B should be done separately.
B ./= B₀
E(x) = SA[0.0 / E₀, 0.0 / E₀, 0.0 / E₀]

## bc=WrapExtrap() uses periodic boundary conditions
param = prepare(;
    E,
    B = build_interpolator(B, (x, y, z); bc = WrapExtrap()),
    m = 1,
    q = 1
)

tspan = (0.0, π) # half averaged gyroperiod based on B₀

## Dummy initial state; positions have units l₀; velocities have units U₀
stateinit = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

prob = ODEProblem(trace!, stateinit, tspan, param)

prob.u0[4:6] = let
    B0 = prob.p.B(prob.u0)
    B0 = normalize(B0)

    Bperp1 = SA[0.0, -B0[3], B0[2]] |> normalize
    Bperp2 = B0 × Bperp1 |> normalize

    ## initial velocity azimuthal angle
    ϕ = 2π * 0
    ## initial velocity pitch angle w.r.t. B
    θ = acos(0.0)

    sinϕ, cosϕ = sincos(ϕ)
    @. (B0 * cos(θ) + Bperp1 * (sin(θ) * cosϕ) + Bperp2 * (sin(θ) * sinϕ)) * U₀
end

saved_values = SavedValues(Float64, Tuple{SVector{3, Float64}, Float64})

function save_B_mu(u, t, integrator)
    b = integrator.p.B(u)
    μ = @views b ⋅ u[4:6] / √(b[1]^2 + b[2]^2 + b[3]^2)

    return b, μ
end

cb = SavingCallback(save_B_mu, saved_values)

sol = solve(prob, Vern9(); callback = cb);

# The extra values are saved in `saved_values`:

saved_values

# ## Fields and Work Rates along a Trajectory
#
# `get_fields(sol)` returns the electric and magnetic fields at every saved state, and `get_work_rates(sol)` the work rates:
#   - $P_{\parallel} = q v_{\parallel} (\mathbf{E} \cdot \hat{b})$: parallel work rate
#   - $P_{\text{Fermi}} = \frac{m v_{\parallel}^2}{B} (\hat{b} \times \boldsymbol{\kappa}) \cdot \mathbf{E}$: Fermi work rate
#   - $P_{\text{Grad}} = \frac{\mu}{B} (\hat{b} \times \nabla B) \cdot \mathbf{E}$: gradient drift work rate
#   - $P_{\text{Betatron}} = \mu \frac{\partial B}{\partial t}$: betatron work rate
#
# Then the total kinetic energy change is $\Delta E_{kin} \approx \int (P_{\parallel} + P_{\text{Fermi}} + P_{\text{Grad}} + P_{\text{Betatron}}) dt$.
#
# Note that $P_{\text{Grad}} = \frac{\mu}{B}(\hat{b} \times \nabla B) \cdot \mathbf{E} = \mu \mathbf{v}_E \cdot \nabla B$,
# where $\mathbf{v}_E = (\mathbf{E} \times \hat{b})/B$ is the ExB drift. It is therefore the
# convective part of the betatron channel rather than a fourth mechanism, and
# [`get_work_channels`](@ref TestParticle.get_work_channels) groups the four terms into
# the three physical channels $(P_\parallel, P_{\text{Fermi}}, P_{\text{Betatron}})$.
#
# Here, we use a constant electric field and design a magnetic field that has gradient, time dependence, and curvature.
# - y dependence creates gradients (gradient drift work)
# - Time dependence creates Betatron work
# - Curvature is implicit in the field structure (Fermi work)

E_boris(x, t) = SA[0.1, 0.1, 0.1]
B_boris(x, t) = SA[x[2], x[1], 1.0] * (1.0 + 0.1 * t)

param_boris = prepare(; E = E_boris, B = B_boris, m = 1, q = 1);

## Create TraceProblem for Boris solver
u0_boris = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0]
tspan_boris = (0.0, 2π)

prob_boris = TraceProblem(u0_boris, tspan_boris, param_boris)

sol_work = solve(prob_boris, Boris(); dt = 0.1);

## Access the magnetic field along the trajectory
_, B_trajectory = get_fields(sol_work);

println("First B field: ", B_trajectory[1])
println("Last B field: ", B_trajectory[end])

# We access the work rates along the trajectory and show at selected time steps:
work = get_work_rates(sol_work)
P_par = [w[1] for w in work]
P_fermi = [w[2] for w in work]
P_grad = [w[3] for w in work]
P_betatron = [w[4] for w in work]

using Printf #hide

println("\nWork Rates Table:") #hide
println("="^70) #hide
@printf("%-10s %-15s %-15s %-15s %-15s\n", "Time", "P_par", "P_fermi", "P_grad", "P_betatron") #hide
println("-"^70) #hide

indices = [1, length(sol_work.u) ÷ 4, length(sol_work.u) ÷ 2, 3 * length(sol_work.u) ÷ 4, length(sol_work.u)] #hide
for i in indices #hide
    t = sol_work.t[i] #hide
    @printf( #hide
        "%-10.2f %-15.3e %-15.3e %-15.3e %-15.3e\n", #hide
        t, P_par[i], P_fermi[i], P_grad[i], P_betatron[i] #hide
    ) #hide
end #hide
println("="^70) #hide
