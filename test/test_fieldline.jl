if !isdefined(Main, :test_common)
    include("test_common.jl")
end

module test_fieldline

using TestParticle
import TestParticle as TP
using OrdinaryDiffEq
using StaticArrays
using Test
using LinearAlgebra: norm
import Magnetostatics as MS
using ..test_common: uniform_grid

const dipole = MS.Dipole(TP.BMoment_Earth)

dipole_field(r) = dipole(r)

"""
    max_L_shell_error(sol, L)

Largest deviation of the L shell from `L` along the traced line, ignoring the
samples too close to the equator or inside the Earth.
"""
function max_L_shell_error(sol, L)
    errs = Float64[]
    for u in sol.u
        r, θ, ϕ = TP.cart2sph(u...)
        if abs(sin(θ)) > 1.0e-2 && r > TP.Rₑ
            push!(errs, abs(r / sin(θ)^2 - L))
        end
    end
    return isempty(errs) ? 0.0 : maximum(errs)
end

@testset "Field line tracing" begin
    @testset "Dipole field" begin
        param = prepare(ZeroField(), dipole_field)
        L = 4.0 * TP.Rₑ
        stateinit = [L, 0.0, 0.0]

        # Use a slightly loose tolerance for cross-platform stability
        # (MacOS M1/M2 vs x86); 1e-2 * Re is about 60 km, under 0.25% of 4 Re.
        tol = 1.0e-2 * TP.Rₑ

        # Trace for a shorter distance to stay well within the valid domain and
        # avoid the singularity at the origin. From L = 4 Re the path to the
        # surface is longer than 4 Re, so 3 Re keeps the trace in space.
        tspan = (0.0, 3.0 * TP.Rₑ)

        # Terminate at Earth's surface to avoid singularity at origin
        isoutside(u, p, t) = norm(u) < TP.Rₑ
        callback = TerminateOutside(isoutside)

        # Forward
        prob = TraceFieldlineProblem(stateinit, param, tspan; mode = :forward)
        sol = solve(prob, Tsit5(); callback)
        @test max_L_shell_error(sol, L) < tol

        # Backward
        prob_back = TraceFieldlineProblem(stateinit, param, tspan; mode = :backward)
        sol_back = solve(prob_back, Tsit5(); callback)
        @test max_L_shell_error(sol_back, L) < tol && prob_back.tspan[2] < 0

        # Both
        probs = TraceFieldlineProblem(stateinit, param, tspan; mode = :both)
        sol1 = solve(probs[1], Tsit5(); callback)
        sol2 = solve(probs[2], Tsit5(); callback)
        @test max_L_shell_error(sol1, L) < tol && max_L_shell_error(sol2, L) < tol

        # Out-of-place equation
        prob_op = ODEProblem(trace_fieldline, SA[stateinit...], tspan, TP.get_BField(param))
        sol_op = solve(prob_op, Tsit5(); callback)
        @test max_L_shell_error(sol_op, L) < tol

        # Function parameter
        prob_func = TraceFieldlineProblem(stateinit, dipole_field, tspan; mode = :forward)
        sol_func = solve(prob_func, Tsit5(); callback)
        @test max_L_shell_error(sol_func, L) < tol

        # AbstractField parameter
        prob_field = TraceFieldlineProblem(
            stateinit, TP.get_BField(param), tspan;
            mode = :forward
        )
        sol_field = solve(prob_field, Tsit5(); callback)
        @test max_L_shell_error(sol_field, L) < tol
    end

    @testset "Interpolated field" begin
        x, y, z = uniform_grid(10, 10, 10)
        B = zeros(3, length(x), length(y), length(z))
        B[1, :, :, :] .= 1.0 # Bx = 1

        param = prepare(x, y, z, zeros(3, length(x), length(y), length(z)), B)
        stateinit = [0.0, 0.0, 0.0]
        tspan = (0.0, 5.0)

        prob = TraceFieldlineProblem(stateinit, param, tspan; mode = :forward)
        sol = solve(prob, Tsit5())

        # Relax tolerance slightly for interpolation artifacts
        @test sol.u[end][1] ≈ 5.0 atol = 5.0e-2
        @test sol.u[end][2] ≈ 0.0 atol = 1.0e-5
        @test sol.u[end][3] ≈ 0.0 atol = 1.0e-5
    end
end

end # module test_fieldline
