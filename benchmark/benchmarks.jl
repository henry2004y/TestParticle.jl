# Benchmark Test for TestParticle

using BenchmarkTools
using TestParticle
using BorisPushers
using TestParticle: LazyTimeInterpolator, get_particle_flux
import TestParticle as TP
using OrdinaryDiffEq, OrdinaryDiffEqLowOrderRK, OrdinaryDiffEqVerner, StaticArrays
using BorisPushers: CPU
using Meshes: Vec, Plane, Disk, Point, Sphere, CartesianGrid, StructuredGrid

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

    param_norm = prepare(; E = E_norm, B = B_norm, m = 1.0, q = 1.0)

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

    B_field = build_interpolator(StructuredGrid, B, r, θ, ϕ)

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
param_analytic = prepare(; E = E_analytic, B = B_analytic)
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
prob_norm_oop = ODEProblem(trace, SA[stateinit_norm...], tspan_norm, param_norm)

SUITE["trace"]["normalized"]["out of place"] = @benchmarkable solve(
    $prob_norm_oop, Tsit5(); save_idxs = [1, 2, 3]
)

# Numerical Field
mesh, E_numeric, B_numeric = setup_numeric_field()
tspan_num = (0.0, 10.0)
param_numeric = prepare(;
    E = build_interpolator(
        E_numeric, mesh;
        bc = TestParticle.FillExtrap(NaN), order = 1
    ),
    B = build_interpolator(
        B_numeric, mesh;
        bc = TestParticle.FillExtrap(NaN), order = 1
    )
)
param_numeric_cubic = prepare(;
    E = build_interpolator(
        E_numeric, mesh;
        bc = TestParticle.FillExtrap(NaN), order = 3
    ),
    B = build_interpolator(
        B_numeric, mesh;
        bc = TestParticle.FillExtrap(NaN), order = 3
    )
)
prob_ip_num = ODEProblem(trace!, stateinit, tspan_num, param_numeric) # in place
prob_oop_num = ODEProblem(trace, SA[stateinit...], tspan_num, param_numeric) # out of place
prob_boris = TraceProblem(stateinit, tspan_num, param_numeric)

SUITE["trace"]["numerical field"]["in place"] = @benchmarkable solve(
    $prob_ip_num, Tsit5(); save_idxs = [1, 2, 3]
)
SUITE["trace"]["numerical field"]["out of place"] = @benchmarkable solve(
    $prob_oop_num, Tsit5(); save_idxs = [1, 2, 3]
)
SUITE["trace"]["numerical field"]["Boris"] = @benchmarkable solve(
    $prob_boris, Boris(); dt = 1 / 7, saveat = 10 / 7
)
SUITE["trace"]["numerical field"]["Boris ensemble"] = @benchmarkable solve(
    $prob_boris, Boris(), EnsembleSerial(); dt = 1 / 7, saveat = 10 / 7, trajectories = 2
)
SUITE["trace"]["numerical field"]["Multistep Boris"] = @benchmarkable solve(
    $prob_boris, MultistepBoris2(; n = 2); dt = 1 / 7, saveat = 10 / 7
)
SUITE["trace"]["numerical field"]["Hyper Boris (n=2, N=4)"] = @benchmarkable solve(
    $prob_boris, MultistepBoris4(; n = 2); dt = 1 / 7, saveat = 10 / 7
)
SUITE["trace"]["numerical field"]["Hyper Boris (n=2, N=6)"] = @benchmarkable solve(
    $prob_boris, MultistepBoris6(; n = 2); dt = 1 / 7, saveat = 10 / 7
)
alg_adaptive = AdaptiveBoris(safety = 0.1)
SUITE["trace"]["numerical field"]["Adaptive Boris"] = @benchmarkable solve(
    $prob_boris, $alg_adaptive
)
alg_adaptive_multi2 = AdaptiveMultistepBoris{2}(n = 2, safety = 0.1)
SUITE["trace"]["numerical field"]["Adaptive Multistep Boris (n=2, N=2)"] = @benchmarkable solve(
    $prob_boris, $alg_adaptive_multi2
)
alg_adaptive_multi4 = AdaptiveMultistepBoris{4}(n = 2, safety = 0.1)
SUITE["trace"]["numerical field"]["Adaptive Hyper Boris (n=2, N=4)"] = @benchmarkable solve(
    $prob_boris, $alg_adaptive_multi4
)
alg_adaptive_multi6 = AdaptiveMultistepBoris{6}(n = 2, safety = 0.1)
SUITE["trace"]["numerical field"]["Adaptive Hyper Boris (n=2, N=6)"] = @benchmarkable solve(
    $prob_boris, $alg_adaptive_multi6
)
SUITE["trace"]["numerical field"]["Boris kernel"] = @benchmarkable solve(
    $prob_boris, Boris(), CPU(); dt = 1 / 7, saveat = 10 / 7
)
isoutside(u, p, t) = hypot(u[1], u[2], u[3]) > 2.0
SUITE["trace"]["numerical field"]["Boris kernel with isoutside"] = @benchmarkable solve(
    $prob_boris, Boris(), CPU(); dt = 1 / 7, saveat = 10 / 7, isoutside = $isoutside
)

# Time-Dependent Field
param_td = prepare(; E = E_td, B = B_td, F = F_td)
prob_ip_td = ODEProblem(trace!, stateinit, tspan, param_td) # in place
prob_oop_td = ODEProblem(trace, SA[stateinit...], tspan, param_td) # out of place

SUITE["trace"]["time-dependent field"]["in place"] = @benchmarkable solve(
    $prob_ip_td, Tsit5(); save_idxs = [1, 2, 3]
)
SUITE["trace"]["time-dependent field"]["out of place"] = @benchmarkable solve(
    $prob_oop_td, Tsit5(); save_idxs = [1, 2, 3]
)

# Interpolation
B_field_car_linear = param_numeric.B
B_field_car_cubic = param_numeric_cubic.B
B_field_sph = setup_spherical_field()
loc = SA[1.0, 1.0, 1.0]
SUITE["interpolation"]["cartesian"] = BenchmarkGroup()
SUITE["interpolation"]["cartesian"]["linear"] = @benchmarkable $B_field_car_linear($loc)
SUITE["interpolation"]["cartesian"]["cubic"] = @benchmarkable $B_field_car_cubic($loc)
SUITE["interpolation"]["spherical"] = @benchmarkable $B_field_sph($loc)

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
    loader_num(i) = build_interpolator(CartesianGrid, B_fields[i], x_grid, y_grid, z_grid)

    LazyTimeInterpolator(times_num, loader_num), 0.5, 2.5
end

SUITE["interpolation"]["time-dependent"] = @benchmarkable begin
    $itp_num($loc, $t1)
    $itp_num($loc, $t2)
end

# GC
stateinit_gc, param_gc = prepare_gc(
    stateinit; E = E_analytic, B = B_analytic, species = Proton
)
prob_gc = ODEProblem(trace_gc!, stateinit_gc, tspan, param_gc)
SUITE["trace"]["GC"]["DiffEq Vern6"] = @benchmarkable solve($prob_gc, Vern6())

prob_native_gc = TraceGCProblem(stateinit_gc, tspan, param_gc)
SUITE["trace"]["GC"]["RK4"] = @benchmarkable solve(
    $prob_native_gc, RK4();
    dt = 2.0e-2, adaptive = false, saveat = 2.0, maxiters = 10000
)
SUITE["trace"]["GC"]["Vern7"] = @benchmarkable solve(
    $prob_native_gc, Vern7();
    dt = 2.0e-2, saveat = 2.0, maxiters = 10000
)

# Hybrid
const B0_hybrid = 1.0e-4
const α_hybrid = 1.0e-2
function bottle_B_hybrid(x, t)
    Bz = B0_hybrid * (1 + α_hybrid * x[3]^2)
    Bx = -B0_hybrid * α_hybrid * x[1] * x[3]
    By = -B0_hybrid * α_hybrid * x[2] * x[3]
    return SA[Bx, By, Bz]
end

# Each hybrid entry is paired with the standalone AdaptiveBoris and Vern7 solvers
# on the same trajectory, so the timings can be read as a comparison. The FO
# baseline tracks the local gyroperiod with the same safety factor as the hybrid's
# FO mode and gets the same `dtmin` / `dtmax` bounds, so both follow one step
# resolution; a fixed step tied to the initial gyroperiod is up to 20x longer than
# the local gyroperiod in the bottle and only aliases the gyration. Diagnostics are
# left at their default (on): a normal solve also records the check history.
let
    B_bottle = TP.Field(bottle_B_hybrid)
    E_bottle = TP.Field((x, t) -> SA[0.0, 0.0, 0.0])

    m, q = TP.mᵢ, TP.qᵢ
    q2m = q / m
    T_gyro = 2π / (abs(q2m) * B0_hybrid)
    tspan = (0.0, 10 * T_gyro)
    p = (q2m, m, E_bottle, B_bottle, TP.ZeroField())
    v0 = SA[5.0e4, 0.0, 1.0e5]
    dt_boris = 0.1 * T_gyro
    dtmax_boris = T_gyro
    dtmin_boris = 1.0e-4 * T_gyro
    alg_boris = AdaptiveBoris(safety = 0.1)

    # Pinned thresholds isolate the cost of each mode. On the axis the field is
    # exactly uniform, so a mixed run would be decided by where the particle
    # happens to be at the first check. The field strength grows by more than an
    # order of magnitude along the bounce, so the full-orbit step follows the
    # local gyroperiod (~700 steps here) instead of aliasing the gyration.
    prob = TraceProblem(vcat(SA[0.0, 0.0, 0.0], v0), tspan, p)
    alg_fo = AdaptiveHybrid(;
        threshold = 0.0, dtmax = T_gyro, dtmin = 1.0e-4 * T_gyro,
        maxiters = 10000, check_interval = 100,
    )
    alg_gc = AdaptiveHybrid(;
        threshold = 1.0e9, dtmax = T_gyro, dtmin = 1.0e-4 * T_gyro,
        maxiters = 10000, check_interval = 100,
    )

    # Off the axis the field lines bend and ε crosses the threshold: the solver
    # switches mode three times on this trajectory.
    u0_mixed = vcat(SA[10.0, 0.0, 0.0], v0)
    prob_mixed = TraceProblem(u0_mixed, tspan, p)
    alg_mixed = AdaptiveHybrid(;
        threshold = 0.1, dtmax = T_gyro, dtmin = 1.0e-4 * T_gyro,
        maxiters = 10000, check_interval = 10,
    )
    state_gc, param_gc = TP.prepare_gc(u0_mixed; B = B_bottle, species = Proton)
    prob_gc = TraceGCProblem(state_gc, tspan, param_gc)

    SUITE["trace"]["Hybrid"]["AdaptiveBoris, magnetic bottle"] = @benchmarkable solve(
        $prob_mixed, $alg_boris;
        dt = $dt_boris, dtmax = $dtmax_boris, dtmin = $dtmin_boris
    )
    SUITE["trace"]["Hybrid"]["Vern7, magnetic bottle"] = @benchmarkable solve(
        $prob_gc, Vern7()
    )
    SUITE["trace"]["Hybrid"]["AdaptiveHybrid, full orbit"] = @benchmarkable solve(
        $prob, $alg_fo; verbose = false, seed = 1234
    )
    SUITE["trace"]["Hybrid"]["AdaptiveHybrid, guiding center"] = @benchmarkable solve(
        $prob, $alg_gc; verbose = false, seed = 1234
    )
    SUITE["trace"]["Hybrid"]["AdaptiveHybrid, mixed"] = @benchmarkable solve(
        $prob_mixed, $alg_mixed; verbose = false, seed = 1234
    )
end

# Hybrid in a current sheet: the particle crosses a field reversal, so the solver
# alternates between the adiabatic and non-adiabatic regimes.
let
    B0_sheet = 1.0e-7
    L_sheet = 2.0e4
    Bg_sheet = 2.0e-8
    sheet_B(x, t) = SA[B0_sheet * tanh(x[3] / L_sheet), 0.0, Bg_sheet]

    m, q = TP.mᵢ, TP.qᵢ
    q2m = q / m
    T_gyro = 2π / (abs(q2m) * B0_sheet)
    tspan = (0.0, 200 * T_gyro)

    v = 2.0e5
    u0 = vcat(SA[0.0, 0.0, -2L_sheet], SA[0.0, 0.95v, 0.3v])
    B_sheet = TP.Field(sheet_B)
    E_zero = TP.Field((x, t) -> SA[0.0, 0.0, 0.0])
    p = (q2m, m, E_zero, B_sheet, TP.ZeroField())
    prob_sheet = TraceProblem(u0, tspan, p)
    dt_boris = T_gyro / 20
    dtmax_boris = dt_boris
    dtmin_boris = 1.0e-4 * T_gyro
    alg_boris = AdaptiveBoris(safety = 0.1)

    alg_sheet = AdaptiveHybrid(;
        threshold = 0.1, dtmax = dt_boris, dtmin = 1.0e-4 * T_gyro,
        maxiters = 500_000, check_interval = 10,
    )
    state_gc, param_gc = TP.prepare_gc(u0; B = B_sheet, species = Proton)
    prob_gc = TraceGCProblem(state_gc, tspan, param_gc)

    SUITE["trace"]["Hybrid"]["AdaptiveBoris, current sheet"] = @benchmarkable solve(
        $prob_sheet, $alg_boris;
        dt = $dt_boris, dtmax = $dtmax_boris, dtmin = $dtmin_boris
    )
    SUITE["trace"]["Hybrid"]["Vern7, current sheet"] = @benchmarkable solve(
        $prob_gc, Vern7()
    )
    SUITE["trace"]["Hybrid"]["AdaptiveHybrid, current sheet"] = @benchmarkable solve(
        $prob_sheet, $alg_sheet; verbose = false, seed = 1234
    )
end

# The cases above are deliberately non-adiabatic stress tests (ε ≳ 0.1): there a
# guiding-center step cannot be much longer than a full-orbit step, but each
# stage needs field derivatives, so standalone GC can lose even to adaptive Boris.
# The opposite regime is the reason to trace guiding centers: with ε ≪ 1 one GC
# step can span many gyroperiods while FO still resolves every gyration. Measure
# it on the same footing, over enough gyroperiods for the gain to show up.
let
    B0_adia = 1.0e-4
    α_adia = 1.0e-8
    adia_B(x, t) = SA[
        -B0_adia * α_adia * x[1] * x[3],
        -B0_adia * α_adia * x[2] * x[3],
        B0_adia * (1 + α_adia * x[3]^2),
    ]

    m, q = TP.mᵢ, TP.qᵢ
    q2m = q / m
    T_gyro = 2π / (abs(q2m) * B0_adia)
    tspan = (0.0, 1000 * T_gyro)

    u0 = vcat(SA[1.0, 0.0, 0.0], SA[5.0e4, 0.0, 1.0e5])
    B_adia = TP.Field(adia_B)
    E_zero = TP.Field((x, t) -> SA[0.0, 0.0, 0.0])
    p = (q2m, m, E_zero, B_adia, TP.ZeroField())
    prob_adia = TraceProblem(u0, tspan, p)
    dt_boris = 0.1 * T_gyro
    dtmax_boris = T_gyro
    dtmin_boris = 1.0e-4 * T_gyro
    alg_boris = AdaptiveBoris(safety = 0.1)

    alg_adia = AdaptiveHybrid(;
        threshold = 0.1, dtmax = T_gyro, dtmin = 1.0e-4 * T_gyro,
        maxiters = 10000, check_interval = 10,
    )
    state_gc, param_gc = TP.prepare_gc(u0; B = B_adia, species = Proton)
    prob_gc = TraceGCProblem(state_gc, tspan, param_gc)

    SUITE["trace"]["Hybrid"]["AdaptiveBoris, adiabatic bottle"] = @benchmarkable solve(
        $prob_adia, $alg_boris;
        dt = $dt_boris, dtmax = $dtmax_boris, dtmin = $dtmin_boris
    )
    SUITE["trace"]["Hybrid"]["Vern7, adiabatic bottle"] = @benchmarkable solve(
        $prob_gc, Vern7()
    )
    SUITE["trace"]["Hybrid"]["AdaptiveHybrid, adiabatic bottle"] = @benchmarkable solve(
        $prob_adia, $alg_adia; verbose = false, seed = 1234
    )
end

# --- Flux Benchmarks ---

struct MockSol{T, U}
    t::T
    u::U
end
# Simple linear-like mock for benchmarking intersection logic
(sol::MockSol)(tc) = SA[tc, 0.0, 0.0, 1.0, 0.0, 0.0]

function setup_flux_benchmarks()
    # 10 sols crossing surfaces
    t = collect(0.0:1.0:10.0)
    u = [SA[x, 0.0, 0.0, 1.0, 0.0, 0.0] for x in t]
    sols = [MockSol(t, u) for _ in 1:10]
    weights = collect(1.1:1.0:10.1)

    # Surfaces
    # Disk at x = 5.0
    disk = Disk(Plane(Point(5.0, 0.0, 0.0), Vec(1.0, 0.0, 0.0)), 2.0)
    # Sphere at x = 8.0
    sphere = Sphere(Point(8.0, 0.0, 0.0), 1.0)

    return sols, weights, disk, sphere
end

SUITE["flux"] = BenchmarkGroup()
let (sols, weights, disk, sphere) = setup_flux_benchmarks()
    # Explicitly type the surface vector to match the function signature
    surfaces = [disk, disk]

    SUITE["flux"]["single detector"] = @benchmarkable get_particle_flux($sols, $disk, $weights)
    SUITE["flux"]["multi detector"] = @benchmarkable get_particle_flux($sols, $surfaces, $weights)
end
