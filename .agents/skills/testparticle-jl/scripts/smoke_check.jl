# Fast smoke check that a TestParticle.jl checkout still traces correctly.
#
#   julia --project=test .agents/skills/testparticle-jl/scripts/smoke_check.jl
#
# Exercises the adaptive full-orbit path, the native Boris pusher and the
# guiding-center solver on the same uniform magnetic field. Exits with a
# non-zero status if any check fails.

using TestParticle, OrdinaryDiffEq, BorisPushers
using StaticArrays, LinearAlgebra
using TestParticle: qᵢ, mᵢ, get_gyrofrequency, get_gyroperiod, get_gyroradius

const B0 = 1.0e-8      # [T], field along +z
const V0 = 1.0e5       # [m/s], perpendicular velocity

uniform_B(x) = SA[0.0, 0.0, B0]
zero_E(x) = SA[0.0, 0.0, 0.0]

const ω = get_gyrofrequency(B0; q = qᵢ, m = mᵢ)
const P = get_gyroperiod(B0; q = qᵢ, m = mᵢ)
const R = get_gyroradius(V0, B0; q = qᵢ, m = mᵢ)

# Start one gyroradius from the origin so the orbit is a circle about the origin.
const xv0 = [-R, 0.0, 0.0, 0.0, V0, 0.0]

const results = NamedTuple{(:name, :ok, :info), Tuple{String, Bool, String}}[]

function check(name, ok, info = "")
    push!(results, (; name, ok, info))
    println(ok ? "  ok   " : "  FAIL ", name, isempty(info) ? "" : "  ($info)")
    return ok
end

println("TestParticle.jl smoke check")
println("  gyroperiod = ", P, " s, gyroradius = ", R, " m")

## Adaptive full orbit: the orbit closes after exactly one gyroperiod ----------
let
    param = prepare(; B = uniform_B, species = Proton)
    prob = TraceProblem(xv0, (0.0, P), param)
    sol = solve(prob, Vern7(); reltol = 1.0e-12, abstol = 1.0e-12)

    u_end = sol.u[end]
    closure = norm(u_end[1:3] - xv0[1:3]) / R
    check("full orbit closes after one gyroperiod", closure < 1.0e-6, "rel err $closure")

    speed = norm(u_end[4:6])
    check("speed conserved", abs(speed - V0) / V0 < 1.0e-9)
end

## Native Boris pusher: energy conservation over ten periods ------------------
let
    param = prepare(; B = uniform_B)
    prob = TraceProblem(xv0, (0.0, 10P), param)
    dt = P / 1000
    sol = solve(prob, Boris(); dt, saveat = P)

    speeds = [norm(u[4:6]) for u in sol.u]
    drift = maximum(abs.(speeds .- V0)) / V0
    check("Boris conserves speed", drift < 1.0e-12, "max rel drift $drift")
    check("Boris produced the expected number of saves", length(sol.u) == 11)
end

## Guiding center: a uniform field leaves the guiding center at rest ----------
let
    state_gc, param_gc = prepare_gc(xv0; B = uniform_B, species = Proton)
    prob = TraceGCProblem(state_gc, (0.0, 5P), param_gc)
    sol = solve(prob, GCRK4(); dt = P / 100)

    X = sol.u[end][1:3]
    check("guiding center stays at rest", norm(X) / R < 1.0e-4, "|X|/ρ = $(norm(X) / R)")
    check("parallel velocity stays zero", abs(sol.u[end][4]) / V0 < 1.0e-6)
end

failed = count(r -> !r.ok, results)
println()
if failed == 0
    println("All $(length(results)) checks passed.")
else
    println("$failed of $(length(results)) checks failed.")
    exit(1)
end
