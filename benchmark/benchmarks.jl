# Benchmark Test for TestParticle

using BenchmarkTools
using TestParticle
import TestParticle as TP
using OrdinaryDiffEq, Meshes, StaticArrays
using TestParticle: CPU

const SUITE = BenchmarkGroup()

SUITE["trace"] = BenchmarkGroup()
SUITE["trace"]["analytic field"] = BenchmarkGroup()
SUITE["trace"]["numerical field"] = BenchmarkGroup()
SUITE["trace"]["normalized"] = BenchmarkGroup()
SUITE["trace"]["time-dependent field"] = BenchmarkGroup()
SUITE["trace"]["GC"] = BenchmarkGroup()
SUITE["trace"]["Hybrid"] = BenchmarkGroup()
SUITE["interpolation"] = BenchmarkGroup()

# --- Field Definitions ---

# Analytic field
E_analytic(xu) = SA[5.0e-10, 5.0e-10, 0]
B_analytic(xu) = SA[0, 0, 1.0e-8]

# Numerical field setup
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

# Normalized field setup
function setup_normalized_field(x0, u0, tspan)
    # Unit conversion factors (Proton default)
    B_dim = 1.0e-8             # [T]
    U_dim = TP.c               # [m/s]
    E_dim = U_dim * B_dim      # [V/m]
    Ω_dim = abs(TP.qᵢ) * B_dim / TP.mᵢ # [1/s]
    t_dim = 1 / Ω_dim        # [s]
    l_dim = U_dim * t_dim    # [m]

    E_mag_norm = 5.0e-10 / E_dim
    B_mag_norm = 1.0e-8 / B_dim

    # Use closures for fields to capture local variables
    E_norm(xu) = SA[E_mag_norm, E_mag_norm, 0.0]
    B_norm(xu) = SA[0.0, 0.0, B_mag_norm]

    # Normalized Initial Conditions
    x0_norm_val = x0 ./ l_dim
    u0_norm_val = u0 ./ U_dim
    stateinit_norm = [x0_norm_val..., u0_norm_val...]
    tspan_norm = tspan ./ t_dim

    param_norm = prepare(E_norm, B_norm; m = 1.0, q = 1.0)

    return param_norm, stateinit_norm, tspan_norm
end

# Time-dependent field definitions
B_td(xu, t) = SA[0, 0, 1.0e-11 * cos(2π * t)]
E_td(xu, t) = SA[5.0e-11 * sin(2π * t), 0, 0]
F_td(xu) = SA[0, 9.10938356e-42, 0]

# Spherical field setup (for interpolation benchmark)
function setup_spherical_field()
    r = logrange(0.1, 10.0, length = 15)
    θ = range(0, π, length = 20)
    ϕ = range(0, 2π, length = 25)

    B₀ = 1.0e-8 # [T]
    B = zeros(3, length(r), length(θ), length(ϕ))

    for (iθ, θ_val) in enumerate(θ)
        sinθ, cosθ = sincos(θ_val)
        B[1, :, iθ, :] .= B₀ * cosθ
        B[2, :, iθ, :] .= -B₀ * sinθ
    end

    B_field = TP.getinterp(TP.StructuredGrid, B, r, θ, ϕ)

    return B_field
end

# FO mode / Switching B field
function sheared_B_func(x)
    B0 = 0.01
    k = 100.0 # High curvature to force FO mode or switching
    return SA[B0 * cos(k * x[1]), B0 * sin(k * x[1]), 0.0]
end

# --- Global Initial Conditions ---
tspan = (0.0, 100.0)
x0 = [0.0, 0.0, 0.0] # initial position, [m]
u0 = [1.0, 0.0, 0.0] # initial velocity, [m/s]
stateinit = [x0..., u0...]

# --- Benchmarks ---

# Analytic Field
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

# Normalized Field
param_norm, stateinit_norm, tspan_norm = setup_normalized_field(x0, u0, tspan)
prob_norm_oop = ODEProblem(trace_normalized, SA[stateinit_norm...], tspan_norm, param_norm)

SUITE["trace"]["normalized"]["out of place"] = @benchmarkable solve(
    $prob_norm_oop, Tsit5(); save_idxs = [1, 2, 3]
)

# Numerical Field
mesh, E_numeric, B_numeric = setup_numeric_field()
tspan_num = (0.0, 10.0)
param_numeric = prepare(mesh, E_numeric, B_numeric)
prob_ip_num = ODEProblem(trace!, stateinit, tspan_num, param_numeric) # in place
prob_oop_num = ODEProblem(trace, SA[stateinit...], tspan_num, param_numeric) # out of place
prob_boris = TraceProblem(stateinit, tspan_num, param_numeric)

SUITE["trace"]["numerical field"]["in place"] = @benchmarkable solve(
    $prob_ip_num, Tsit5(); save_idxs = [1, 2, 3]
)
SUITE["trace"]["numerical field"]["out of place"] = @benchmarkable solve(
    $prob_oop_num, Tsit5(); save_idxs = [1, 2, 3]
)
SUITE["trace"]["numerical field"]["Boris"] = @benchmarkable TP.solve(
    $prob_boris; dt = 1 / 7, savestepinterval = 10
)
SUITE["trace"]["numerical field"]["Boris with fields"] = @benchmarkable TP.solve(
    $prob_boris; dt = 1 / 7, savestepinterval = 10, save_fields = true
)
SUITE["trace"]["numerical field"]["Boris ensemble"] = @benchmarkable TP.solve(
    $prob_boris; dt = 1 / 7, savestepinterval = 10, trajectories = 2
)
SUITE["trace"]["numerical field"]["Multistep Boris"] = @benchmarkable TP.solve(
    $prob_boris; dt = 1 / 7, savestepinterval = 10, n = 2
)
alg_adaptive = AdaptiveBoris(dtmax = 1.0e-3)
SUITE["trace"]["numerical field"]["Adaptive Boris"] = @benchmarkable TP.solve(
    $prob_boris, $alg_adaptive
)
SUITE["trace"]["numerical field"]["Boris kernel"] = @benchmarkable TP.solve(
    $prob_boris, CPU(); dt = 1 / 7, savestepinterval = 10
)

# Time-Dependent Field
param_td = prepare(E_td, B_td, F_td)
prob_ip_td = ODEProblem(trace!, stateinit, tspan, param_td) # in place
prob_oop_td = ODEProblem(trace, SA[stateinit...], tspan, param_td) # out of place

SUITE["trace"]["time-dependent field"]["in place"] = @benchmarkable solve(
    $prob_ip_td, Tsit5(); save_idxs = [1, 2, 3]
)
SUITE["trace"]["time-dependent field"]["out of place"] = @benchmarkable solve(
    $prob_oop_td, Tsit5(); save_idxs = [1, 2, 3]
)

# Interpolation
B_field_car = TP.get_BField(param_numeric)
B_field_sph = setup_spherical_field()
loc = SA[1.0, 1.0, 1.0]
SUITE["interpolation"]["cartesian"] = @benchmarkable $B_field_car(loc)
SUITE["interpolation"]["spherical"] = @benchmarkable $B_field_sph(loc)

# Lazy Time-Dependent Interpolation
itp_num, t1, t2 = let
    # Define grid
    x_grid = range(0.0, 1.0, length = 4)
    y_grid = range(0.0, 1.0, length = 4)
    z_grid = range(0.0, 1.0, length = 4)

    B_fields = ntuple(_ -> fill(0.0, 3, 4, 4, 4), 4)
    B_fields[1][1, :, :, :] .= 1.0
    B_fields[2][2, :, :, :] .= 2.0
    B_fields[3][3, :, :, :] .= 3.0
    B_fields[4][1, :, :, :] .= 4.0

    times_num = [0.0, 1.0, 2.0, 3.0]

    # Loader for numerical field
    loader_num(i) = TP.getinterp(TP.CartesianGrid, B_fields[i], x_grid, y_grid, z_grid)

    LazyTimeInterpolator(times_num, loader_num), 0.5, 2.5
end

SUITE["interpolation"]["time-dependent"] = @benchmarkable begin
    $itp_num($loc, $t1)
    $itp_num($loc, $t2)
end

# GC
stateinit_gc, param_gc = TP.prepare_gc(
    stateinit, E_analytic, B_analytic,
    species = Proton
)
prob_gc = ODEProblem(trace_gc!, stateinit_gc, tspan, param_gc)
SUITE["trace"]["GC"]["DiffEq Vern6"] = @benchmarkable solve($prob_gc, Vern6())

prob_native_gc = TraceGCProblem(stateinit_gc, tspan, param_gc)
SUITE["trace"]["GC"]["Native RK4"] = @benchmarkable TP.solve(
    $prob_native_gc;
    dt = 1.0e-4, savestepinterval = 100, alg = :rk4
)
SUITE["trace"]["GC"]["Native RK45"] = @benchmarkable TP.solve(
    $prob_native_gc;
    dt = 1.0e-4, savestepinterval = 100, alg = :rk45
)

# Hybrid
alg_hybrid = AdaptiveHybrid(threshold = 0.1, dtmax = 1.0e-4)
param_sheared = TP.prepare(E_analytic, sheared_B_func) # Reuse E_analytic

x0_h = SA[0.0, 0.0, 0.0]
v0_h = SA[1.0e4, 0.0, 1.0e3]
u0_h = vcat(x0_h, v0_h)
tspan_h = (0.0, 1.0e-5)

prob_hybrid_sheared = TraceHybridProblem(u0_h, tspan_h, param_sheared)

SUITE["trace"]["Hybrid"]["Sheared"] = @benchmarkable TP.solve(
    $prob_hybrid_sheared, $alg_hybrid
)
