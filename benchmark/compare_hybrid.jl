# Compare the adaptive hybrid solver against a fixed-step Boris pusher and the
# guiding-center solver on two mixed-regime problems.
#
# The reference is a full-orbit Boris trace with a very small step; every solver
# is compared against it in position only, since the guiding-center state carries
# no gyrophase. Run with
#
#     julia --project=test benchmark/compare_hybrid.jl
#
# The output is a Markdown table per scenario: wall time (best of repeated runs),
# full-orbit mode fraction (hybrid only), kinetic-energy drift and the RMS
# position error against the reference. Position errors grow with the number of
# gyroperiods because the bounce phase is exponentially sensitive, so the energy
# drift is the robust accuracy measure.

using TestParticle, BorisPushers, OrdinaryDiffEq, StaticArrays
using LinearAlgebra, Printf
import TestParticle as TP

const E_ZERO = TP.Field((x, t) -> SA[0.0, 0.0, 0.0])

"Best wall time [ms] of `run()` over `n` calls, after one warm-up call."
function timed(run, n = 3)
    run()
    best = Inf
    for _ in 1:n
        best = min(best, @elapsed run())
    end
    return 1.0e3 * best
end

"RMS position error of `sol` against `ref` on `n` uniformly spaced times."
function rms_error(sol, ref, tspan, n)
    err = 0.0
    for t in range(tspan[1], tspan[2]; length = n)
        err += sum(abs2, TP.get_x(sol(t)) - TP.get_x(ref(t)))
    end
    return sqrt(err / n)
end

"Relative drift of the kinetic energy over the trajectory, or `NaN` if unknown."
function energy_drift(sol, m)
    u_first, u_last = sol.u[1], sol.u[end]
    length(u_first) == 6 || return NaN
    ke(u) = 0.5 * m * sum(abs2, u[SA[4, 5, 6]])
    return abs(ke(u_last) - ke(u_first)) / ke(u_first)
end

fo_frac(sol) = sol.stats === nothing ? NaN :
    count(==(:FO), sol.stats.adiabaticity.mode) / length(sol.stats.adiabaticity.mode)

"""
Harris sheet, `B = (B0 tanh(z/L), 0, Bg)`. The guide field `Bg` keeps `|B|` finite
at the reversal, so a guiding center stays defined while crossing.
"""
harris_B(B0, L, Bg) = (x, t) -> SA[B0 * tanh(x[3] / L), 0.0, Bg]

"Bottle field of the hybrid example, where the curvature peaks off the axis."
function bottle_B(B0, α)
    return function (x, t)
        Bz = B0 * (1 + α * x[3]^2)
        Bx = -B0 * α * x[1] * x[3]
        By = -B0 * α * x[2] * x[3]
        return SA[Bx, By, Bz]
    end
end

function compare(name, Bfunc, u0, tspan, dt_reference, dt_boris, hybrid_kw)
    q, m = TP.qᵢ, TP.mᵢ
    p = (q / m, m, E_ZERO, TP.Field(Bfunc), TP.ZeroField())
    prob = TraceProblem(u0, tspan, p)

    ref = solve(prob, Boris(); dt = dt_reference, save_everystep = true)
    boris = solve(prob, Boris(); dt = dt_boris, save_everystep = true)

    alg = AdaptiveHybrid(; dtmax = dt_boris, hybrid_kw...)
    hybrid = solve(prob, alg; verbose = false).u[1]

    state_gc, param_gc = TP.prepare_gc(u0; B = Bfunc, species = Proton)
    gc_prob = TraceGCProblem(state_gc, tspan, param_gc)
    gc = solve(gc_prob, GCRK45()).u[1]

    rows = (
        ("Boris", timed(() -> solve(prob, Boris(); dt = dt_boris)), NaN, boris),
        ("GCRK45", timed(() -> solve(gc_prob, GCRK45())), NaN, gc),
        (
            "AdaptiveHybrid", timed(() -> solve(prob, alg; verbose = false)),
            fo_frac(hybrid), hybrid,
        ),
    )

    println("\n### ", name, "\n")
    println("| Solver | Time [ms] | FO fraction | KE drift | Steps | RMS error [m] |")
    println("| :--- | ---: | ---: | ---: | ---: | ---: |")
    for (label, t, frac, sol) in rows
        frac_str = isnan(frac) ? "-" : @sprintf("%.2f", frac)
        drift = energy_drift(sol, m)
        drift_str = isnan(drift) ? "-" : @sprintf("%.1e", drift)
        @printf(
            "| %s | %.2f | %s | %s | %d | %.3e |\n",
            label, t, frac_str, drift_str, length(sol.t),
            rms_error(sol, ref, tspan, 200)
        )
    end
    return nothing
end

function main()
    # 1. Magnetic bottle: adiabatic near the mirrors, non-adiabatic on the
    #    curved midplane off the axis. Same setup as the hybrid example.
    let B0 = 1.0e-4, α = 1.0e-2, T_gyro = 2π / (abs(TP.qᵢ / TP.mᵢ) * 1.0e-4)
        u0 = vcat(SA[10.0, 0.0, 0.0], SA[5.0e4, 0.0, 1.0e5])
        compare(
            "Magnetic bottle, 100 gyroperiods", bottle_B(B0, α), u0,
            (0.0, 100 * T_gyro), T_gyro / 200, T_gyro / 20,
            (; threshold = 0.1, dtmin = 1.0e-4 * T_gyro, maxiters = 500_000,
                check_interval = 10),
        )
    end

    # 2. Harris sheet: the particle crosses the field reversal, so it is close to
    #    the adiabaticity limit in the sheet and magnetized outside of it.
    let B0 = 2.0e-7, L = 2.0e4, Bg = 2.0e-8, v = 4.0e4,
        T_gyro = 2π / (abs(TP.qᵢ / TP.mᵢ) * 2.0e-7)

        u0 = vcat(SA[0.0, 0.0, -2L], SA[0.0, 0.95v, 0.3v])
        compare(
            "Harris sheet, 100 gyroperiods", harris_B(B0, L, Bg), u0,
            (0.0, 100 * T_gyro), T_gyro / 200, T_gyro / 20,
            (; threshold = 0.1, dtmin = 1.0e-4 * T_gyro, maxiters = 500_000,
                check_interval = 10),
        )
    end
    return nothing
end

main()
