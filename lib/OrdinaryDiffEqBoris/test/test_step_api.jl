module test_step_api

using Test
using OrdinaryDiffEqBoris
using StaticArrays
using LinearAlgebra: norm

zero_E(x, t) = SA[0.0, 0.0, 0.0]
uniform_B(x, t) = SA[0.0, 0.0, 0.01]

"The Energy of a particle of unit mass, in a potential-free field."
kinetic_energy(v) = 0.5 * sum(abs2, v)

function integrate(alg, p, t0, dt, nsteps)
    r = SA[0.0, 0.0, 0.0]
    v = SA[1.0e5, 0.0, 0.0]
    t = t0
    v_half = boris_half_velocity(v, r, dt, t, p, alg)
    v_sum = zero(eltype(v))

    for _ in 1:nsteps
        r, v_half = boris_advance(v_half, r, dt, t, p, alg)
        t += dt
        v_node = boris_node_velocity(v_half, r, dt, t, p, alg)
        v_sum += v_node[1]
    end

    return v_sum
end

@testset "Step API" begin
    @testset "the step functions reproduce the solve loop" begin
        # Driving the step functions by hand has to walk the same trajectory the
        # SciML loop walks, since that is what any other driver gets.
        let param = (-1.75882001076e11, 9.1093837015e-31, zero_E, uniform_B),
            u0 = SA[0.0, 0.0, 0.0, 1.0e7, 0.0, 0.0], tspan = (0.0, 1.0e-8), dt = 1.0e-11

            for alg in (Boris(), MultistepBoris2(n = 2), MultistepBoris4(n = 2))
                sol = solve(ODEProblem((u, p, t) -> nothing, u0, tspan, param), alg; dt)

                r = SA[u0[1], u0[2], u0[3]]
                v = SA[u0[4], u0[5], u0[6]]
                t = tspan[1]
                v_half = boris_half_velocity(v, r, dt, t, param, alg)

                for n in 1:(length(sol.t) - 1)
                    r, v_half = boris_advance(v_half, r, dt, t, param, alg)
                    t += dt
                    v_node = boris_node_velocity(v_half, r, dt, t, param, alg)
                    @test vcat(r, v_node) ≈ sol.u[n + 1] rtol = 1.0e-10
                end
            end
        end
    end

    @testset "recentering the velocity keeps its magnitude" begin
        # A change of step size goes through the node velocity. Every leg of the
        # round trip is a rotation, so in a magnetic field alone the speed has to
        # survive it, which is what keeps the scheme time-reversible.
        let param = (1.0, 1.0, zero_E, uniform_B),
            r = SA[0.0, 0.0, 0.0], v = SA[1.0e5, 2.0e5, -3.0e5], t = 0.0, dt = 1.0e-9

            for alg in (Boris(), MultistepBoris6(n = 3))
                v_half = boris_half_velocity(v, r, dt, t, param, alg)

                for new_dt in (dt, 2dt, dt / 4)
                    recentered = boris_resync_velocity(v_half, r, dt, new_dt, t, param, alg)
                    v_node = boris_node_velocity(recentered, r, new_dt, t, param, alg)
                    @test kinetic_energy(v_node) ≈ kinetic_energy(v) rtol = 1.0e-12
                end
            end
        end
    end

    @testset "a step is free of allocation" begin
        # Nothing here may allocate if the same functions are to run inside a
        # kernel, where an allocation is not merely slow but impossible.
        let param = (1.0, 1.0, zero_E, uniform_B), dt = 1.0e-9

            for alg in (Boris(), MultistepBoris2(n = 2), MultistepBoris4(n = 2),
                MultistepBoris6(n = 4))
                integrate(alg, param, 0.0, dt, 10)
                @test @allocated(integrate(alg, param, 0.0, dt, 1000)) == 0
            end
        end
    end
end

end
