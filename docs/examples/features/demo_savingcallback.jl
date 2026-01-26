# # Additional Diagnostics
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

function getmeanB(B)
    B₀sum = eltype(B)(0)
    for k in axes(B, 4), j in axes(B, 3), i in axes(B, 2)
        B₀sum += B[1, i, j, k]^2 + B[2, i, j, k]^2 + B[3, i, j, k]^2
    end

    return sqrt(B₀sum / prod(size(B)[2:4]))
end

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
const B₀ = getmeanB(B)
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
#
# When `save_fields=true`, the solution will include 6 additional values per time step:
# - E field components (Ex, Ey, Ez) at indices 7, 8, 9
# - B field components (Bx, By, Bz) at indices 10, 11, 12
#
# When `save_work=true`, the solution will include 4 additional values per time step representing the work rates:
# - P_par: parallel work rate (index 7 or 13 depending on save_fields)
# - P_fermi: Fermi work rate (index 8 or 14)
# - P_grad: gradient drift work rate (index 9 or 15)
# - P_betatron: betatron work rate (index 10 or 16)

## Set up a simple test problem
x_boris = range(-0.5, 0.5, length = 4)
y_boris = range(-0.5, 0.5, length = 4)
z_boris = range(-0.5, 0.5, length = 4)

B_boris = zeros(Float32, 3, 4, 4, 4)
B_boris[3, :, :, :] .= 1.0

E_boris(x) = SA[0.0, 0.0, 0.0]

param_boris = prepare(x_boris, y_boris, z_boris, E_boris, B_boris; m = 1, q = 1, bc = 2);

## Create TraceProblem for Boris solver
u0_boris = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0]
tspan_boris = (0.0, 2π)

prob_boris = TraceProblem(u0_boris, tspan_boris, param_boris)

## Solve with field saving enabled
sol_fields = solve(prob_boris; dt = 0.1, save_fields = true);

# The solution now contains 12 values per time step (6 for state + 6 for fields)
## Access the magnetic field along the trajectory
B_trajectory = [u[10:12] for u in sol_fields.u];

## Verify that the saved B field matches the expected values
println("First saved B field: ", B_trajectory[1])
println("Last saved B field: ", B_trajectory[end])

# Solve with work saving enabled
sol_work = solve(prob_boris; dt = 0.1, save_work = true);

# The solution now contains 10 values per time step (6 for state + 4 for work rates)
## Access the work rates along the trajectory
P_par = [u[7] for u in sol_work.u]
P_fermi = [u[8] for u in sol_work.u]
P_grad = [u[9] for u in sol_work.u]
P_betatron = [u[10] for u in sol_work.u];

using Printf #hide

println("\nWork Rates Table:") #hide
println("="^70) #hide
@printf("%-10s %-15s %-15s %-15s %-15s\n", "Time", "P_par", "P_fermi", "P_grad", "P_betatron") #hide
println("-"^70) #hide

# Show work rates at selected time steps
indices = [1, length(sol_work.u) ÷ 4, length(sol_work.u) ÷ 2, 3 * length(sol_work.u) ÷ 4, length(sol_work.u)]
for i in indices #hide
    t = sol_work.t[i] #hide
    @printf( #hide
        "%-10.3f %-15.6e %-15.6e %-15.6e %-15.6e\n", #hide
        t, P_par[i], P_fermi[i], P_grad[i], P_betatron[i] #hide
    ) #hide
end #hide
println("="^70) #hide
