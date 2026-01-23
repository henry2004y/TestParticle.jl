using TestParticle
using Test
using StaticArrays
using LinearAlgebra
using ForwardDiff
using SciMLBase: solve, ReturnCode


@testset "Boris Pusher Extras" begin
    # Defines a simple field
    B0 = 1.0
    E0 = 0.1
    L = 10.0

    # B = B0 * (1 + x/L) z_hat
    function Bfunc(x, t)
        return SVector{3}(0.0, 0.0, B0 * (1 + x[1] / L))
    end

    # E = E0 y_hat
    function Efunc(x, t)
        return SVector{3}(0.0, E0, 0.0)
    end

    # Simple particle
    x0 = [1.0, 0.0, 0.0]
    v0 = [0.0, 1.0, 1.0]
    stateinit = [x0..., v0...]
    tspan = (0.0, 2π) # One gyroperiod approx

    # Use prepare to get parameters with mass!
    # Using Proton mass/charge
    param_tuple = prepare(Efunc, Bfunc; species = Proton)
    # param_tuple is (q2m, m, E, B, F)

    prob = TraceProblem(stateinit, tspan, param_tuple)

    # Test 1: Standard solve (backward compatibility)
    sol = solve(prob, dt = 0.1)
    @test length(sol[1].u[1]) == 6
    @test sol[1].retcode == ReturnCode.Default


    # Test 2: Save fields
    sol_fields = solve(prob, dt = 0.1, save_fields = true)
    @test length(sol_fields[1].u[1]) == 12
    # Check E and B
    r1 = sol_fields[1].u[1][1:3]
    E_saved = sol_fields[1].u[1][7:9]
    B_saved = sol_fields[1].u[1][10:12]
    @test E_saved ≈ Efunc(r1, 0.0)
    @test B_saved ≈ Bfunc(r1, 0.0)

    # Test 3: Save work
    sol_work = solve(prob, dt = 0.1, save_work = true)
    @test length(sol_work[1].u[1]) == 10
    # Check work terms are finite
    # P_par, P_gradB, P_curv, P_ind
    # v0 has parallel component (vz=1, B is z), so P_par = q * vz * Ez. Ez=0. So P_par should be 0.
    # Wait, E is y_hat. B is z_hat.
    # v_par is along z. E is along y. v_par . E = 0.

    work_terms = sol_work[1].u[1][7:10]
    @test work_terms[1] ≈ 0.0 atol = 1.0e-10 # P_par

    # P_gradB: gradB is along x. B is along z.
    # ∇B = (B0/L, 0, 0).
    # b x ∇B = z x x = y.
    # v_gradB is along y.
    # E is along y.
    # So P_gradB should be non-zero (positive).
    @test work_terms[2] > 0.0 # P_gradB

    # P_curv: Field lines are straight (only magnitude changes). Curvature is 0.
    # P_curv should be 0.
    @test work_terms[3] ≈ 0.0 atol = 1.0e-10

    # P_ind: Field is static. dB/dt = 0.
    @test work_terms[4] ≈ 0.0 atol = 1.0e-10

    # Test 4: Save both
    sol_both = solve(prob, dt = 0.1, save_fields = true, save_work = true)
    @test length(sol_both[1].u[1]) == 16

    println("Boris Pusher Extras tests passed!")
end
