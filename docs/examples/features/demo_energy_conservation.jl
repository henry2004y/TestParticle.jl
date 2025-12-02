# # Energy Conservation
#
# This example demonstrates the energy conservation (or theoretical agreement) of a single particle motion in three cases.
# 1. Constant B field, Zero E field.
# 2. Constant E field, Zero B field.
# 3. Magnetic Mirror.
#
# The tests are performed in dimensionless units with q=1, m=1.
# We compare multiple solvers: Tsit5, Vern7, Vern9, BS3, ImplicitMidpoint, Boris, and Boris Multistep.

import DisplayAs #hide
using TestParticle
using TestParticle: ZeroField
using OrdinaryDiffEq
using StaticArrays
using LinearAlgebra: ×, norm
using CairoMakie
CairoMakie.activate!(type = "png") #hide

const q = 1.0
const m = 1.0
const B₀ = 1.0
const E₀ = 1.0

# Solvers to test
# Dictionary of Name => Algorithm
ode_solvers = [
    ("Tsit5", Tsit5()),
    ("Vern7", Vern7()),
    ("Vern9", Vern9()),
    ("BS3", BS3()),
    ("ImplicitMidpoint", ImplicitMidpoint())
]

# Helper function to run tests
function run_test(case_name, param, x0, v0, tspan, expected_energy_func; dt=0.1)
    u0 = [x0..., v0...]
    prob_ode = ODEProblem(trace_normalized!, u0, tspan, param)
    prob_tp = TraceProblem(u0, tspan, param)

    f = Figure(size = (1000, 600), fontsize = 18)
    ax = Axis(f[1, 1],
        title = "$case_name: Energy Error",
        xlabel = "Time",
        ylabel = "Rel. Energy Error |(E - E_ref)/E_ref|",
        yscale = log10
    )

    # Run ODE solvers
    for (name, alg) in ode_solvers
        # Catch potential errors with solvers (e.g. divergence)
        try
            sol = solve(prob_ode, alg)

            # Calculate energy
            v_mag = [norm(u[4:6]) for u in sol.u]
            E = 0.5 .* m .* v_mag.^2

            # Expected energy
            t = sol.t
            x = [u[1:3] for u in sol.u]
            # Pass velocity to expected_energy_func just in case
            E_ref = [expected_energy_func(ti, xi, u[4:6]) for (ti, xi, u) in zip(t, x, sol.u)]

            # Error
            # Avoid division by zero if E_ref is 0
            error = abs.(E .- E_ref) ./ (abs.(E_ref) .+ 1e-16)

            lines!(ax, t, error, label = name)
        catch e
            println("Solver $name failed: $e")
        end
    end

    # Run Native Solvers
    # Boris
    sol_boris = TestParticle.solve(prob_tp; dt=dt)[1] # returns Vector{TraceSolution}
    v_mag_b = [norm(u[4:6]) for u in sol_boris.u]
    E_b = 0.5 .* m .* v_mag_b.^2
    t_b = sol_boris.t
    x_b = [u[1:3] for u in sol_boris.u]
    E_ref_b = [expected_energy_func(ti, xi, u[4:6]) for (ti, xi, u) in zip(t_b, x_b, sol_boris.u)]
    error_b = abs.(E_b .- E_ref_b) ./ (abs.(E_ref_b) .+ 1e-16)
    lines!(ax, t_b, error_b, label = "Boris (dt=$dt)")

    # Boris Multistep (n=2)
    sol_multi = TestParticle.solve(prob_tp; dt=dt, n=2)[1]
    v_mag_m = [norm(u[4:6]) for u in sol_multi.u]
    E_m = 0.5 .* m .* v_mag_m.^2
    t_m = sol_multi.t
    x_m = [u[1:3] for u in sol_multi.u]
    E_ref_m = [expected_energy_func(ti, xi, u[4:6]) for (ti, xi, u) in zip(t_m, x_m, sol_multi.u)]
    error_m = abs.(E_m .- E_ref_m) ./ (abs.(E_ref_m) .+ 1e-16)
    lines!(ax, t_m, error_m, label = "Boris Multistep (n=2)")

    axislegend(ax, position = :rb) # bottom right to avoid covering start
    return f
end

# ## Case 1: Constant B, Zero E
# Energy should be conserved.

uniform_B(x) = SA[0, 0, B₀]

param1 = prepare(ZeroField(), uniform_B; species = User, q = q, m = m)
x0_1 = [0.0, 0.0, 0.0]
v0_1 = [1.0, 0.0, 0.0]
tspan1 = (0.0, 50.0)
E_func1(t, x, v) = 0.5 * m * norm(v0_1)^2 # Constant energy

f1 = run_test("Constant B", param1, x0_1, v0_1, tspan1, E_func1)
f1 = DisplayAs.PNG(f1) #hide

# ## Case 2: Constant E, Zero B
# Energy increases due to work done by the electric field.
# For a particle starting from rest in constant E field:
# $\mathbf{a} = \frac{q}{m} \mathbf{E}$
# $\mathbf{v}(t) = \mathbf{a} t$ (if $v_0=0$)
# $E_{kin} = \frac{1}{2} m v^2 = \frac{1}{2} m (\frac{q E_0}{m} t)^2$

constant_E(x) = SA[E₀, 0.0, 0.0]

param2 = prepare(constant_E, ZeroField(); species = User, q = q, m = m)
x0_2 = [0.0, 0.0, 0.0]
v0_2 = [0.0, 0.0, 0.0] # Start from rest
tspan2 = (0.0, 5.0)

function E_func2(t, x, v)
    # Analytical energy
    v_theo = (q * E₀ / m) * t
    return 0.5 * m * v_theo^2
end

f2 = run_test("Constant E", param2, x0_2, v0_2, tspan2, E_func2)
f2 = DisplayAs.PNG(f2) #hide

# ## Case 3: Magnetic Mirror
# Energy should be conserved (E=0).
# The particle bounces back and forth between regions of high magnetic field.

α = 0.1
function mirror_B(x)
    Bz = B₀ * (1 + α * x[3]^2)

    # Divergence-free B field in cylindrical symmetry
    # Bx = -0.5 * x * dBz/dz
    # By = -0.5 * y * dBz/dz
    # dBz/dz = 2 * B0 * alpha * z

    Bx = - B₀ * α * x[1] * x[3]
    By = - B₀ * α * x[2] * x[3]
    return SA[Bx, By, Bz]
end

param3 = prepare(ZeroField(), mirror_B; species = User, q = q, m = m)
x0_3 = [0.1, 0.0, 0.0]
v0_3 = [0.5, 0.5, 1.0]
tspan3 = (0.0, 100.0)
E_init_3 = 0.5 * m * norm(v0_3)^2
E_func3(t, x, v) = E_init_3

f3 = run_test("Magnetic Mirror", param3, x0_3, v0_3, tspan3, E_func3)
f3 = DisplayAs.PNG(f3) #hide

println("Demo finished successfully")
