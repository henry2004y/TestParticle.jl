module test_step_api

using Test
using BorisPushers
import BorisPushers: update_velocity_resync
using StaticArrays
using LinearAlgebra: norm

zero_E(x, t) = SA[0.0, 0.0, 0.0]
uniform_B(x, t) = SA[0.0, 0.0, 0.01]
oscillating_E(x, t) = SA[sin(2π * t), 0.0, 0.0]

const EVALUATIONS = Ref(0)
counting_E(x, t) = (EVALUATIONS[] += 1; SA[0.0, 0.0, 0.0])

"The Energy of a particle of unit mass, in a potential-free field."
kinetic_energy(v) = 0.5 * sum(abs2, v)

function integrate(alg, p, t0, dt, nsteps)
    r = SA[0.0, 0.0, 0.0]
    v = SA[1.0e5, 0.0, 0.0]
    t = t0
    v_half = update_velocity_half(v, r, dt, t, p, alg)
    v_sum = zero(eltype(v))

    for _ in 1:nsteps
        r, v_half = advance_boris(v_half, r, dt, t, p, alg)
        t += dt
        v_node = update_velocity_node(v_half, r, dt, t, p, alg)
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
                v_half = update_velocity_half(v, r, dt, t, param, alg)

                for n in 1:(length(sol.t) - 1)
                    r, v_half = advance_boris(v_half, r, dt, t, param, alg)
                    t += dt
                    v_node = update_velocity_node(v_half, r, dt, t, param, alg)
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
                v_half = update_velocity_half(v, r, dt, t, param, alg)

                for new_dt in (dt, 2dt, dt / 4)
                    recentered = update_velocity_resync(
                        v_half, r, dt, new_dt, t, param, alg
                    )
                    v_node = update_velocity_node(
                        recentered, r, new_dt, t, param, alg
                    )
                    @test kinetic_energy(v_node) ≈ kinetic_energy(v) rtol = 1.0e-12
                end
            end
        end
    end

    @testset "sampling the fields at the node keeps the method second order" begin
        # A step advances the velocity over [t - dt/2, t + dt/2], the midpoint of
        # which is the node t, so the fields have to be evaluated there. Taking
        # them at t + dt/2 instead, the end of that interval, is the rectangle
        # rule in place of the midpoint rule and costs an order, which only shows
        # up once the fields move in time.
        let T = 1.0, param = (1.0, 1.0, oscillating_E, uniform_B),
                u0 = SA[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

            final_state(dt) = solve(
                ODEProblem((u, p, t) -> nothing, u0, (0.0, T), param), Boris();
                dt, save_everystep = false, save_start = false
            ).u[end]

            reference = final_state(T / 8000)
            errors = [norm(final_state(dt) - reference) for dt in (T / 250, T / 500, T / 1000)]
            orders = [log2(errors[i] / errors[i + 1]) for i in eachindex(errors[1:(end - 1)])]

            @test all(>(1.8), orders)
        end
    end

    @testset "a step costs one field evaluation" begin
        # A solver asks for the fields at a node twice, once to synchronise the
        # state at the end of a step and once to advance the next one from that
        # same node, so the pair is carried across the step boundary.
        let nsteps = 500, dt = 1.0e-3, param = (1.0, 1.0, counting_E, uniform_B),
                prob = ODEProblem(
                (u, p, t) -> nothing, SA[0.0, 0.0, 0.0, 1.0e5, 0.0, 0.0],
                (0.0, nsteps * dt), param
            )

            for alg in (Boris(), MultistepBoris4(n = 3))
                EVALUATIONS[] = 0
                solve(prob, alg; dt, save_everystep = false, save_start = false)
                # One to start the staggered velocity, then one per step.
                @test EVALUATIONS[] == nsteps + 1
            end
        end
    end

    @testset "the fields are taken again when the state is moved" begin
        # A callback can move the particle between steps, and fields carried
        # across such a move would belong to where it used to be.
        let nsteps = 500, dt = 1.0e-3, param = (1.0, 1.0, counting_E, uniform_B),
                prob = ODEProblem(
                (u, p, t) -> nothing, SA[0.0, 0.0, 0.0, 1.0e5, 0.0, 0.0],
                (0.0, nsteps * dt), param
            ), moved = Ref(false)

            callback = DiscreteCallback(
                (u, t, integrator) -> t >= 0.25 && !moved[],
                function (integrator)
                    moved[] = true
                    integrator.u = SA[1.0, 0.0, 0.0, 1.0e5, 0.0, 0.0]
                    return
                end
            )

            EVALUATIONS[] = 0
            solve(prob, Boris(); dt, callback, save_everystep = false, save_start = false)
            @test EVALUATIONS[] == nsteps + 2
        end
    end

    @testset "a step is free of allocation" begin
        # Nothing here may allocate if the same functions are to run inside a
        # kernel, where an allocation is not merely slow but impossible.
        let param = (1.0, 1.0, zero_E, uniform_B), dt = 1.0e-9

            for alg in (
                    Boris(), MultistepBoris2(n = 2), MultistepBoris4(n = 2),
                    MultistepBoris6(n = 4),
                )
                integrate(alg, param, 0.0, dt, 10)
                @test @allocated(integrate(alg, param, 0.0, dt, 1000)) == 0
            end
        end
    end
end

end
