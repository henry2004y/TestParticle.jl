# Where does the time of the adaptive hybrid solver go?
#
# The hybrid solver is decomposed into its two modes by forcing each of them for
# the whole trajectory: `threshold = 0` keeps the full-orbit (Boris) mode, a huge
# threshold keeps the guiding-center mode. Comparing those two runs against the
# standalone Boris pusher and the standalone guiding-center solver isolates the
# fixed per-mode overhead of the hybrid loop. The micro-benchmarks price the
# individual operations so the overhead can be attributed, and the last table
# sweeps `check_interval` to show the cost of the adiabaticity check itself.
#
# Run with
#
#     julia --project=test benchmark/hybrid_overhead.jl

using TestParticle, BorisPushers, OrdinaryDiffEq, StaticArrays
using LinearAlgebra, Printf
import TestParticle as TP
import SciMLBase

const E_ZERO = TP.Field((x, t) -> SA[0.0, 0.0, 0.0])

"""
    per_call(f, n = 10_000) -> (ns, sink)

Cost of one `f(delta)` call in nanoseconds, averaged over `n` calls. `f` must
return a number and must depend on its argument, otherwise the compiler hoists
the whole call out of the loop; `delta` sweeps a tiny offset so every iteration
does distinct work. The accumulated result is returned so the loop cannot be
optimised away.
"""
function per_call(f::F, n = 10_000) where {F}
    f(0.0)
    acc = 0.0
    t = @elapsed for i in 1:n
        acc += f(1.0e-6 * Float64(i % 16))
    end
    return 1.0e9 * t / n, acc
end

"Print one row of the per-operation table; `f` must depend on its argument."
function report(f::F, name) where {F}
    t, acc = per_call(f)
    f(0.0)
    bytes = @allocated f(1.0e-6)
    @printf("| `%s` | %.1f | %d |\n", name, t, bytes)
    return acc == -Inf && error("unreachable")
end

"Best wall time [ms] of `f()`, after a warm-up call."
function time_ms(f, n = 5)
    f()
    best = Inf
    for _ in 1:n
        best = min(best, @elapsed f())
    end
    return 1.0e3 * best
end

function harris_B(B0, L, Bg)
    return (x, t) -> SA[B0 * tanh(x[3] / L), 0.0, Bg]
end

function micro_benchmarks(Bfunc, u0, tspan)
    q, m = TP.qᵢ, TP.mᵢ
    q2m = q / m
    E, B = E_ZERO, TP.Field(Bfunc)
    x, v = TP.get_x(u0), TP.get_v(u0)
    t0 = tspan[1]

    X, vpar, _, _, μ, phase = TP._get_gc_parameters(u0, E, B, q, m, t0)
    state_gc = SVector(X[1], X[2], X[3], vpar)
    p_gc = (; q, q2m, μ, E, B)
    p = (; q2m, m, E, B, F = TP.ZeroField())
    fields = BorisPushers._node_fields(p, x, t0)
    dt = 1.0e-2

    shift3 = SVector(1.0, 0.0, 0.0)
    shift4 = SVector(1.0, 0.0, 0.0, 0.0)
    shift6 = SVector(1.0, 0.0, 0.0, 0.0, 0.0, 0.0)

    println("\n### Per-operation cost\n")
    println("| Operation | Time [ns] | Allocated [B] |")
    println("| :--- | ---: | ---: |")

    # Each `report` call is specialised on its own closure type. Iterating a
    # heterogeneous tuple of closures instead would dispatch dynamically and
    # box the arguments, which inflates every measurement.
    report("adiabaticity_components") do δ
        sum(TP.adiabaticity_components(X + δ * shift3, B, q, m, μ, t0))
    end
    report("get_magnetic_properties") do δ
        sum(TP.get_magnetic_properties(X + δ * shift3, t0, B).κ)
    end
    report("_adiabaticity_components (cached mp)") do δ
        sum(TP._adiabaticity_components(
            TP.get_magnetic_properties(X + δ * shift3, t0, B), q, m, μ))
    end
    report("_gc_adia_parameters (X, μ only)") do δ
        Xa, μa = TP._gc_adia_parameters(u0 + δ * shift6, E, B, q, m, t0)
        sum(Xa) + μa
    end
    report("_gc_to_full") do δ
        sum(TP._gc_to_full(state_gc + δ * shift4, E, B, q, m, μ, t0, phase))
    end
    report("_get_gc_parameters") do δ
        sum(TP._get_gc_parameters(u0 + δ * shift6, E, B, q, m, t0)[1])
    end
    report("_gc_state_at_switch (2 it)") do δ
        Xs, vp, μs, ph = TP._gc_state_at_switch(u0 + δ * shift6, E, B, q, m, t0, 2)
        sum(Xs) + vp + μs + ph
    end
    report("update_dp5 (one GC step)") do δ
        dy, dE = TP.update_dp5(state_gc + δ * shift4, p_gc, dt, t0)
        sum(dy) + sum(dE)
    end
    report("advance_boris (one FO step)") do δ
        r, vn = BorisPushers.advance_boris(v, x + δ * shift3, dt, t0, fields, Boris())
        sum(r) + sum(vn)
    end
    report("update_velocity_node") do δ
        sum(BorisPushers.update_velocity_node(v, x + δ * shift3, dt, t0, fields, Boris()))
    end
    return nothing
end

function decompose(name, Bfunc, u0, tspan, dt_boris, hybrid_kw)
    q, m = TP.qᵢ, TP.mᵢ
    p = (q / m, m, E_ZERO, TP.Field(Bfunc), TP.ZeroField())
    prob = TraceProblem(u0, tspan, p)

    fo_alg = AdaptiveHybrid(; threshold = 0.0, dtmax = dt_boris, hybrid_kw...)
    gc_alg = AdaptiveHybrid(; threshold = 1.0e9, dtmax = dt_boris, hybrid_kw...)
    mix_alg = AdaptiveHybrid(; threshold = 0.1, dtmax = dt_boris, hybrid_kw...)

    state_gc, param_gc = TP.prepare_gc(u0; B = Bfunc, species = Proton)
    gc_prob = TraceGCProblem(state_gc, tspan, param_gc)

    runs = (
        ("Boris (same dt)", () -> solve(prob, Boris(); dt = dt_boris)),
        ("GCRK45", () -> solve(gc_prob, GCRK45())),
        ("Hybrid, FO only", () -> solve(prob, fo_alg; verbose = false)),
        ("Hybrid, GC only", () -> solve(prob, gc_alg; verbose = false)),
        ("Hybrid, mixed", () -> solve(prob, mix_alg; verbose = false)),
    )

    println("\n### ", name, "\n")
    println("| Run | Time [ms] | Saved points | Checks | ns/point |")
    println("| :--- | ---: | ---: | ---: | ---: |")
    for (label, f) in runs
        t = time_ms(f)
        sol = f()
        sol = sol isa SciMLBase.EnsembleSolution ? sol.u[1] : sol
        npts = length(sol.t)
        checks = hasproperty(sol.stats, :adiabaticity) ?
            length(sol.stats.adiabaticity.t) : 0
        @printf("| %s | %.2f | %d | %d | %.0f |\n", label, t, npts, checks,
            1.0e6 * t / npts)
    end
    return nothing
end

function check_sweep(name, Bfunc, u0, tspan, dt_boris, threshold, hybrid_kw)
    q, m = TP.qᵢ, TP.mᵢ
    p = (q / m, m, E_ZERO, TP.Field(Bfunc), TP.ZeroField())
    prob = TraceProblem(u0, tspan, p)

    println("\n### ", name, " — cost of the adiabaticity check\n")
    println("| `check_interval` | Time [ms] | Saved points | Checks | Switches |")
    println("| ---: | ---: | ---: | ---: | ---: |")
    for ci in (1, 10, 100, 1000)
        alg = AdaptiveHybrid(; threshold, dtmax = dt_boris, check_interval = ci,
            hybrid_kw...)
        t = time_ms(() -> solve(prob, alg; verbose = false))
        sol = solve(prob, alg; verbose = false).u[1]
        modes = sol.stats.adiabaticity.mode
        switches = count(i -> modes[i] != modes[i - 1], 2:length(modes))
        @printf("| %d | %.2f | %d | %d | %d |\n", ci, t, length(sol.t),
            length(modes), switches)
    end
    return nothing
end

function main()
    # Harris sheet: the particle oscillates across the field reversal, so it
    # keeps alternating between the adiabatic and non-adiabatic regimes.
    B0 = 2.0e-7
    L = 2.0e4
    Bg = 2.0e-8
    v = 4.0e4
    T_gyro = 2π / (abs(TP.qᵢ / TP.mᵢ) * B0)
    u0 = vcat(SA[0.0, 0.0, -2L], SA[0.0, 0.95v, 0.3v])
    dt_boris = T_gyro / 20
    kw = (; dtmin = 1.0e-4 * T_gyro, maxiters = 2_000_000, check_interval = 10)

    micro_benchmarks(harris_B(B0, L, Bg), u0, (0.0, 100 * T_gyro))
    decompose("Harris sheet, 500 gyroperiods", harris_B(B0, L, Bg), u0,
        (0.0, 500 * T_gyro), dt_boris, kw)
    sweep_kw = (; dtmin = 1.0e-4 * T_gyro, maxiters = 2_000_000)

    # Full-orbit only: the trajectory does not depend on `check_interval`, so the
    # time difference is purely the cost of the checks themselves.
    check_sweep("Harris sheet, 500 gyroperiods, full orbit only",
        harris_B(B0, L, Bg), u0, (0.0, 500 * T_gyro), dt_boris, 0.0, sweep_kw)

    check_sweep("Harris sheet, 500 gyroperiods, switching",
        harris_B(B0, L, Bg), u0, (0.0, 500 * T_gyro), dt_boris, 0.1, sweep_kw)
    return nothing
end

main()
