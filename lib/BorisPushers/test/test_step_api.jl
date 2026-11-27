module test_step_api

using Test
using BorisPushers
import BorisPushers: update_velocity_resync, update_velocity_boris
using StaticArrays
using LinearAlgebra: norm

function _loop_update_velocity_boris_rel(u, E, B, qdt2m, inv_c2, n)
    for _ in 1:n
        u = update_velocity_boris(u, E, B, qdt2m, inv_c2)
    end
    return u
end

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
        let param = (; q2m = -1.75882001076e11, E = zero_E, B = uniform_B),
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
        let param = (; q2m = 1.0, E = zero_E, B = uniform_B),
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
        let T = 1.0, param = (; q2m = 1.0, E = oscillating_E, B = uniform_B),
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
        let nsteps = 500, dt = 1.0e-3, param = (; q2m = 1.0, E = counting_E, B = uniform_B),
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
        let nsteps = 500, dt = 1.0e-3, param = (; q2m = 1.0, E = counting_E, B = uniform_B),
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
        let param = (; q2m = 1.0, E = zero_E, B = uniform_B), dt = 1.0e-9

            for alg in (
                    Boris(), MultistepBoris2(n = 2), MultistepBoris4(n = 2),
                    MultistepBoris6(n = 4),
                )
                integrate(alg, param, 0.0, dt, 10)
                @test @allocated(integrate(alg, param, 0.0, dt, 1000)) == 0
            end
        end
    end

    @testset "relativistic step API" begin
        c0 = 2.99792458e8
        param_rel = (; q2m = -1.75882001076e11, E = zero_E, B = uniform_B, c = c0)
        u_init = 0.6 * c0 / sqrt(1 - 0.6^2)
        u0_rel = SA[0.0, 0.0, 0.0, u_init, 0.0, 0.0]
        dt = 1.0e-11
        tspan = (0.0, 1.0e-8)

        for alg in (
                Boris(), MultistepBoris2(; n = 2), MultistepBoris4(; n = 2),
                Vay(), HigueraCary(),
            )
            sol = solve(ODEProblem((u, p, t) -> nothing, u0_rel, tspan, param_rel), alg; dt)

            r = SA[u0_rel[1], u0_rel[2], u0_rel[3]]
            v = SA[u0_rel[4], u0_rel[5], u0_rel[6]]
            t = tspan[1]
            v_half = update_velocity_half(v, r, dt, t, param_rel, alg)

            for n in 1:(length(sol.t) - 1)
                r, v_half = advance_boris(v_half, r, dt, t, param_rel, alg)
                t += dt
                v_node = update_velocity_node(v_half, r, dt, t, param_rel, alg)
                @test vcat(r, v_node) ≈ sol.u[n + 1] rtol = 1.0e-10
            end
        end

        for alg in (Boris(), MultistepBoris6(; n = 3), Vay(), HigueraCary())
            v = SA[u_init, 0.2 * c0, -0.3 * c0]
            r = SA[0.0, 0.0, 0.0]
            v_half = update_velocity_half(v, r, dt, 0.0, param_rel, alg)
            for new_dt in (dt, 2dt, dt / 4)
                recentered = update_velocity_resync(
                    v_half, r, dt, new_dt, 0.0, param_rel, alg
                )
                v_node = update_velocity_node(recentered, r, new_dt, 0.0, param_rel, alg)
                @test sum(abs2, v_node) ≈ sum(abs2, v) rtol = 1.0e-12
            end
        end

        param_td = (; q2m = 1.0, E = oscillating_E, B = uniform_B, c = 1.0)
        u0_td = SA[0.0, 0.0, 0.0, 0.8, 0.0, 0.0]
        T_td = 1.0
        for alg in (Boris(), Vay(), HigueraCary())
            final_state_rel(step) = solve(
                ODEProblem((u, p, t) -> nothing, u0_td, (0.0, T_td), param_td), alg;
                dt = step, save_everystep = false, save_start = false
            ).u[end]
            ref = final_state_rel(T_td / 4000)
            errs = [
                norm(final_state_rel(step) - ref)
                    for step in (T_td / 200, T_td / 400, T_td / 800)
            ]
            orders = [log2(errs[i] / errs[i + 1]) for i in 1:(length(errs) - 1)]
            @test all(>(1.8), orders)
        end

        for alg in (Boris(), MultistepBoris4(; n = 2), Vay(), HigueraCary())
            integrate(alg, param_rel, 0.0, dt, 10)
            @test @allocated(integrate(alg, param_rel, 0.0, dt, 1000)) == 0
        end

        @testset "relativistic update_velocity_boris" begin
            let r = SA[0.0, 0.0, 0.0],
                    B_si = uniform_B(r, 0.0),
                    E0 = zero_E(r, 0.0),
                    qdt2m_si = param_rel.q2m * 0.5 * dt,
                    inv_c2_si = inv(c0^2),
                    u_si = SA[u_init, 0.2 * c0, -0.3 * c0]
                u_step = update_velocity_boris(u_si, E0, B_si, qdt2m_si, inv_c2_si)
                u_alg = update_velocity(u_si, r, dt, 0.0, param_rel, Boris())
                @test u_step ≈ u_alg rtol = 1.0e-14
                @test sum(abs2, u_step) ≈ sum(abs2, u_si) rtol = 1.0e-14

                E_n = SA[0.05, -0.02, 0.01]
                B_n = SA[0.0, 0.0, 1.0]
                for FT in (Float64, Float32)
                    u_n = SVector{3, FT}(0.4, -0.3, 0.2)
                    E_T = SVector{3, FT}(E_n)
                    B_T = SVector{3, FT}(B_n)
                    qdt2m_T = FT(0.05)
                    u_nonrel = update_velocity_boris(u_n, E_T, B_T, qdt2m_T)
                    u_c_inf = update_velocity_boris(u_n, E_T, B_T, qdt2m_T, zero(FT))
                    @test u_c_inf == u_nonrel

                    u_B = update_velocity_boris(u_n, zero(E_T), B_T, qdt2m_T, one(FT))
                    @test sum(abs2, u_B) ≈ sum(abs2, u_n) rtol = 10 * eps(FT)
                end

                _loop_update_velocity_boris_rel(u_si, E0, B_si, qdt2m_si, inv_c2_si, 10)
                @test @allocated(
                    _loop_update_velocity_boris_rel(u_si, E0, B_si, qdt2m_si, inv_c2_si, 100)
                ) == 0
            end
        end
    end
end

end
