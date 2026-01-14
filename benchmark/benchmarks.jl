# Benchmark Test for TestParticle

using BenchmarkTools
using TestParticle
import TestParticle as TP
using OrdinaryDiffEq, Meshes, StaticArrays

const SUITE = BenchmarkGroup()

SUITE["trace"] = BenchmarkGroup()
SUITE["trace"]["analytic field"] = BenchmarkGroup()
SUITE["trace"]["numerical field"] = BenchmarkGroup()
SUITE["trace"]["time-dependent field"] = BenchmarkGroup()
SUITE["interpolation"] = BenchmarkGroup()

# analytic field
E_analytic(xu) = SA[5.0e-10, 5.0e-10, 0]
B_analytic(xu) = SA[0, 0, 1.0e-8]

function setup_numeric_field()
    x = range(-10, 10, length = 15)
    y = range(-10, 10, length = 20)
    z = range(-10, 10, length = 25)
    B_numeric = fill(0.0, 3, length(x), length(y), length(z))
    E_numeric = fill(0.0, 3, length(x), length(y), length(z))
    B_numeric[3, :, :, :] .= 10.0e-9
    E_numeric[1, :, :, :] .= 5.0e-10
    E_numeric[2, :, :, :] .= 5.0e-10

    mesh = CartesianGrid(
        (first(x), first(y), first(z)), (last(x), last(y), last(z));
        dims = (length(x) - 1, length(y) - 1, length(z) - 1)
    )

    return mesh, E_numeric, B_numeric
end

# numerical field
mesh, E_numeric, B_numeric = setup_numeric_field()

B_td(xu, t) = SA[0, 0, 1.0e-11 * cos(2π * t)]
E_td(xu, t) = SA[5.0e-11 * sin(2π * t), 0, 0]
F_td(xu) = SA[0, 9.10938356e-42, 0]

tspan = (0.0, 100.0)
x0 = [0.0, 0.0, 0.0] # initial position, [m]
u0 = [1.0, 0.0, 0.0] # initial velocity, [m/s]
stateinit = [x0..., u0...]

param_analytic = prepare(E_analytic, B_analytic)
prob_ip = ODEProblem(trace!, stateinit, tspan, param_analytic) # in place
prob_rel_ip = ODEProblem(trace_relativistic!, stateinit, tspan, param_analytic) # in place
prob_oop = ODEProblem(trace, SA[stateinit...], tspan, param_analytic) # out of place
SUITE["trace"]["analytic field"]["in place"] = @benchmarkable solve(
    $prob_ip, Tsit5(); save_idxs = [1, 2, 3]
)
SUITE["trace"]["analytic field"]["out of place"] = @benchmarkable solve(
    $prob_oop, Tsit5(); save_idxs = [1, 2, 3]
)
SUITE["trace"]["analytic field"]["in place relativistic"] = @benchmarkable solve(
    $prob_rel_ip, Tsit5(); save_idxs = [1, 2, 3]
)

tspan = (0.0, 10.0)
param_numeric = prepare(mesh, E_numeric, B_numeric)
prob_ip = ODEProblem(trace!, stateinit, tspan, param_numeric) # in place
prob_oop = ODEProblem(trace, SA[stateinit...], tspan, param_numeric) # out of place
prob_boris = TraceProblem(stateinit, tspan, param_numeric)

SUITE["trace"]["numerical field"]["in place"] = @benchmarkable solve(
    $prob_ip, Tsit5(); save_idxs = [1, 2, 3]
)
SUITE["trace"]["numerical field"]["out of place"] = @benchmarkable solve(
    $prob_oop, Tsit5(); save_idxs = [1, 2, 3]
)
SUITE["trace"]["numerical field"]["Boris"] = @benchmarkable TP.solve(
    $prob_boris; dt = 1 / 7, savestepinterval = 10
)
SUITE["trace"]["numerical field"]["Boris ensemble"] = @benchmarkable TP.solve(
    $prob_boris; dt = 1 / 7, savestepinterval = 10, trajectories = 2
)

function setup_spherical_field()
    r = logrange(0.1, 10.0, length = 15)
    θ = range(0, π, length = 20)
    ϕ = range(0, 2π, length = 25)

    B₀ = 1.0e-8 # [nT]
    B = zeros(3, length(r), length(θ), length(ϕ))

    for (iθ, θ_val) in enumerate(θ)
        sinθ, cosθ = sincos(θ_val)
        B[1, :, iθ, :] .= B₀ * cosθ
        B[2, :, iθ, :] .= -B₀ * sinθ
    end

    B_field = TP.getinterp(TP.StructuredGrid, B, r, θ, ϕ)

    return B_field
end

B_field_car = TP.get_BField(param_numeric)
B_field_sph = setup_spherical_field()
loc = SA[1.0, 1.0, 1.0]
SUITE["interpolation"]["cartesian"] = @benchmarkable $B_field_car(loc)
SUITE["interpolation"]["spherical"] = @benchmarkable $B_field_sph(loc)

param_td = prepare(E_td, B_td, F_td)
prob_ip = ODEProblem(trace!, stateinit, tspan, param_td) # in place
prob_oop = ODEProblem(trace, SA[stateinit...], tspan, param_td) # out of place
SUITE["trace"]["time-dependent field"]["in place"] = @benchmarkable solve(
    $prob_ip, Tsit5(); save_idxs = [1, 2, 3]
)
SUITE["trace"]["time-dependent field"]["out of place"] = @benchmarkable solve(
    $prob_oop, Tsit5(); save_idxs = [1, 2, 3]
)

stateinit_gc,
    param_gc = TP.prepare_gc(
    stateinit, E_analytic, B_analytic,
    species = Proton, removeExB = true
)
prob_gc = ODEProblem(trace_gc!, stateinit_gc, tspan, param_gc)
SUITE["trace"]["GC"]["1st order"] = @benchmarkable solve($prob_gc, Vern9())

# Hybrid solver benchmark
# Define a simple analytic magnetic mirror field
function B_mirror(x, t)
    B0_hybrid = 1.0
    L_hybrid = 2.0
    xx, yy, zz = x[1], x[2], x[3]
    Bx = -B0_hybrid * xx * zz / L_hybrid^2
    By = -B0_hybrid * yy * zz / L_hybrid^2
    Bz = B0_hybrid * (1 + zz^2 / L_hybrid^2)
    return SA[Bx, By, Bz]
end
# 1-arg version for initialization
B_mirror(x) = B_mirror(x, 0.0)

E_zero_hybrid(x, t) = SA[0.0, 0.0, 0.0]
E_zero_hybrid(x) = SA[0.0, 0.0, 0.0]

# Initial condition for hybrid
x0_hybrid = SA[0.5, 0.0, 0.5]
Ek_hybrid = 1.0e6 * TP.qᵢ # 1 MeV
v0_hybrid = TP.energy2velocity(Ek_hybrid; m = TP.Proton.m, q = TP.Proton.q)
v0_vec_hybrid = SA[0.0, v0_hybrid * sind(45), v0_hybrid * cosd(45)]

state_gc_hybrid, params_gc_hybrid = TP.prepare_gc(vcat(x0_hybrid, v0_vec_hybrid), E_zero_hybrid, B_mirror; species = TP.Proton)

prob_gc_hybrid = ODEProblem(trace_gc!, state_gc_hybrid, tspan, params_gc_hybrid)

SUITE["trace"]["GC"]["hybrid"] = @benchmarkable solve_hybrid(
    $prob_gc_hybrid, Tsit5(); epsilon = 0.1, dt = 1.0e-7
)
