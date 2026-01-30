# # Extra Diagnostics
#
# This example demonstrates tracing one proton in an analytic E field and numerical B field.
#
# The `SavingCallback` from DiffEqCallbacks.jl can be used to save additional outputs for diagnosis. Here we save the magnetic field along the trajectory, together with the parallel velocity.
# Note that `SavingCallback` is currently not compatible with ensemble problems; for multiple particle tracing with customized outputs, see [Demo: ensemble tracing with extra saving](@ref Ensemble-Tracing).
#
# For the native Boris solvers, we also support `save_fields` and `save_work` keywords to save the fields and work done by the fields. These options are more efficient than using `SavingCallback` and are compatible with ensemble problems.

using TestParticle, OrdinaryDiffEqVerner, StaticArrays
using TestParticle: qᵢ, mᵢ
using TestParticle: get_BField
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
const B₀ = get_mean_magnitude(B)
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

## bc=2 uses periodic boundary conditions
param = prepare(x, y, z, E, B; m = 1, q = 1, bc = 2)

tspan = (0.0, π) # half averaged gyroperiod based on B₀

## Dummy initial state; positions have units l₀; velocities have units U₀
stateinit = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

prob = ODEProblem(trace_normalized!, stateinit, tspan, param)

prob.u0[4:6] = let
    B0 = get_BField(prob)(prob.u0)
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
    b = get_BField(integrator.p)(u)
    μ = @views b ⋅ u[4:6] / √(b[1]^2 + b[2]^2 + b[3]^2)

    return b, μ
end

cb = SavingCallback(save_B_mu, saved_values)

sol = solve(prob, Vern9(); callback = cb);

# The extra values are saved in `saved_values`:

saved_values

# ## Native Boris Solver with Additional Diagnostics
#
# The native Boris solver supports additional diagnostic outputs through the `save_fields` and `save_work` keywords.
# These options allow you to save the electric and magnetic fields along the trajectory, as well as various work rate components without using callbacks.
# Here, we set an E field with constant non-zero to generate parallel work, and a B field that is spatially varying and time-dependent.
#
# When `save_fields=true`, the solution will include 6 additional values per time step:
# - E field components (Ex, Ey, Ez) at indices 7, 8, 9
# - B field components (Bx, By, Bz) at indices 10, 11, 12
# When `save_work=true`, the solution will include 4 additional values per time step representing the work rates:
# - $P_{\parallel} = q v_{\parallel} (\mathbf{E} \cdot \hat{b})$: parallel work rate (index 7 or 13 depending on save_fields)
# - $P_{\text{Fermi}} = \frac{m v_{\parallel}^2}{B} (\hat{b} \times \boldsymbol{\kappa}) \cdot \mathbf{E}$: Fermi work rate (index 8 or 14)
# - $P_{\text{Grad}} = \frac{\mu}{B} (\hat{b} \times \nabla B) \cdot \mathbf{E}$: gradient drift work rate (index 9 or 15)
# - $P_{\text{Betatron}} = \mu \frac{\partial B}{\partial t}$: betatron work rate (index 10 or 16)
# Then the total kinetic energy change is $\Delta E_{kin} \approx \int (P_{\parallel} + P_{\text{Fermi}} + P_{\text{Grad}} + P_{\text{Betatron}}) dt$.
#
# Here, we use a constant electric field and design a magnetic field that has gradient, time dependence, and curvature.
# - y dependence creates gradients (gradient drift work)
# - Time dependence creates Betatron work
# - Curvature is implicit in the field structure (Fermi work)

E_boris(x, t) = SA[0.1, 0.1, 0.1]
B_boris(x, t) = SA[x[2], x[1], 1.0] * (1.0 + 0.1 * t)

param_boris = prepare(E_boris, B_boris; m = 1, q = 1);

## Create TraceProblem for Boris solver
u0_boris = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0]
tspan_boris = (0.0, 2π)

prob_boris = TraceProblem(u0_boris, tspan_boris, param_boris)

## Solve with field saving enabled
sol_fields = TestParticle.solve(prob_boris; dt = 0.1, save_fields = true)[1];

# The solution now contains 12 values per time step (6 for state + 6 for fields)
## Access the magnetic field along the trajectory
B_trajectory = [u[10:12] for u in sol_fields.u];

## Verify that the saved B field matches the expected values
println("First saved B field: ", B_trajectory[1])
println("Last saved B field: ", B_trajectory[end])

## Solve with work saving enabled
sol_work = TestParticle.solve(prob_boris; dt = 0.1, save_work = true)[1];

# The solution now contains 10 values per time step (6 for state + 4 for work rates)
# We access the work rates along the trajectory and show at selected time steps:
P_par = [u[7] for u in sol_work.u]
P_fermi = [u[8] for u in sol_work.u]
P_grad = [u[9] for u in sol_work.u]
P_betatron = [u[10] for u in sol_work.u]

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

# ## Post-processing with `get_fields` and `get_work`
#
# If you didn't save the fields or work during the simulation (e.g., to save memory or if you decided to inspect them later), you can compute them from the trajectory solution using `get_fields` and `get_work`.
# Note: This requires the problem parameters (fields) to be accessible from the solution object.

## Solving without saving extra data
sol = TestParticle.solve(prob_boris; dt = 0.1)[1];

## Compute fields and work post-simulation
E_post, B_post = get_fields(sol);
work_post = get_work(sol);

println("\nPost-processed Data Check:") #hide
println("Computed $(length(E_post)) field points.") #hide
println("Computed $(length(work_post)) work points.") #hide
println("Differences at step 10: P_par: ", abs(work_post[10][1] - P_par[10]), ", P_fermi: ", abs(work_post[10][2] - P_fermi[10]), ", P_grad: ", abs(work_post[10][3] - P_grad[10]), ", P_betatron: ", abs(work_post[10][4] - P_betatron[10])) #hide
