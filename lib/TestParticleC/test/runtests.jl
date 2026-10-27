module test_testparticlec

using Test
using TestParticleC
using StaticArrays
using LinearAlgebra

@testset "TestParticleC C-ABI" begin
    @testset "Single Particle Push Float64" begin
        let
            x = [0.0, 0.0, 0.0]
            v = [1.0e5, 0.0, 0.0]
            E = [0.0, 0.0, 0.0]
            B = [0.0, 0.0, 1.0] # 1 Tesla along z
            q2m = 1.0e8
            dt = 1.0e-9

            TestParticleC.testparticle_boris_push_f64(
                pointer(x), pointer(v), q2m, dt, pointer(E), pointer(B)
            )

            # In pure magnetic field, kinetic energy is strictly conserved
            @test norm(v) ≈ 1.0e5 atol = 1.0e-8
            @test norm(x) > 0.0
        end
    end

    @testset "Single Particle Push Float32" begin
        let
            x = Float32[0.0, 0.0, 0.0]
            v = Float32[1.0e5, 0.0, 0.0]
            E = Float32[0.0, 0.0, 0.0]
            B = Float32[0.0, 0.0, 1.0]
            q2m = 1.0f8
            dt = 1.0f-9

            TestParticleC.testparticle_boris_push_f32(
                pointer(x), pointer(v), q2m, dt, pointer(E), pointer(B)
            )

            @test norm(v) ≈ 1.0f5 atol = 0.1f0
            @test norm(x) > 0.0f0
        end
    end

    @testset "Staggered Push" begin
        let
            x = [0.0, 0.0, 0.0]
            v_half = [1.0e5, 0.0, 0.0]
            E = [0.0, 0.0, 0.0]
            B = [0.0, 0.0, 1.0]
            q2m = 1.0e8
            dt = 1.0e-9

            TestParticleC.testparticle_boris_step_staggered_f64(
                pointer(x), pointer(v_half), q2m, dt, pointer(E), pointer(B)
            )

            @test norm(v_half) ≈ 1.0e5 atol = 1.0e-8
            @test x ≈ v_half * dt
        end
    end

    @testset "Batch Push Float64" begin
        let
            n = 10
            x = zeros(3 * n)
            v = zeros(3 * n)
            E = zeros(3 * n)
            B = zeros(3 * n)

            for i in 0:(n - 1)
                v[3 * i + 1] = 1.0e5 # vx
                B[3 * i + 3] = 1.0   # Bz
            end

            TestParticleC.testparticle_boris_push_batch_f64(
                n, pointer(x), pointer(v), 1.0e8, 1.0e-9, pointer(E), pointer(B)
            )

            for i in 0:(n - 1)
                vi = [v[3 * i + 1], v[3 * i + 2], v[3 * i + 3]]
                @test norm(vi) ≈ 1.0e5 atol = 1.0e-8
            end
        end
    end

    @testset "Batch Push Float32" begin
        let
            n = 5
            x = zeros(Float32, 3 * n)
            v = zeros(Float32, 3 * n)
            E = zeros(Float32, 3 * n)
            B = zeros(Float32, 3 * n)

            for i in 0:(n - 1)
                v[3 * i + 1] = 1.0f5
                B[3 * i + 3] = 1.0f0
            end

            TestParticleC.testparticle_boris_push_batch_f32(
                n, pointer(x), pointer(v), 1.0f8, 1.0f-9, pointer(E), pointer(B)
            )

            for i in 0:(n - 1)
                vi = Float32[v[3 * i + 1], v[3 * i + 2], v[3 * i + 3]]
                @test norm(vi) ≈ 1.0f5 atol = 0.1f0
            end
        end
    end

    @testset "Update Velocity Boris" begin
        let
            v = [1.0e5, 0.0, 0.0]
            E = [0.0, 1.0e2, 0.0]
            B = [0.0, 0.0, 1.0e-2]
            qdt2m = 1.0e8 * 1.0e-9 / 2

            TestParticleC.testparticle_update_velocity_boris_f64(
                pointer(v), pointer(E), pointer(B), qdt2m
            )
            @test norm(v) > 0.0
        end
    end

    @testset "Multistep Boris Push" begin
        let
            for order in (2, 4, 6)
                x = [0.0, 0.0, 0.0]
                v = [1.0e5, 0.0, 0.0]
                E = [0.0, 0.0, 0.0]
                B = [0.0, 0.0, 1.0]
                q2m = 1.0e8
                dt = 1.0e-9

                TestParticleC.testparticle_multistep_boris_push_f64(
                    pointer(x), pointer(v), q2m, dt, pointer(E), pointer(B),
                    Cint(2), Cint(order)
                )
                @test norm(v) ≈ 1.0e5 atol = 1.0e-8
                @test norm(x) > 0.0
            end
        end
    end

    @testset "Hyper-Boris 4 and 6 Specialized" begin
        let
            x4 = [0.0, 0.0, 0.0]
            v4 = [1.0e5, 0.0, 0.0]
            x6 = [0.0, 0.0, 0.0]
            v6 = [1.0e5, 0.0, 0.0]
            E = [0.0, 0.0, 0.0]
            B = [0.0, 0.0, 1.0]
            q2m = 1.0e8
            dt = 1.0e-9

            TestParticleC.testparticle_hyper_boris4_push_f64(
                pointer(x4), pointer(v4), q2m, dt, pointer(E), pointer(B), Cint(2)
            )
            TestParticleC.testparticle_hyper_boris6_push_f64(
                pointer(x6), pointer(v6), q2m, dt, pointer(E), pointer(B), Cint(2)
            )

            @test norm(v4) ≈ 1.0e5 atol = 1.0e-8
            @test norm(v6) ≈ 1.0e5 atol = 1.0e-8
            @test norm(x4) > 0.0
            @test norm(x6) > 0.0
        end
    end

    @testset "Multistep Staggered and Batch" begin
        let
            x = [0.0, 0.0, 0.0]
            v_half = [1.0e5, 0.0, 0.0]
            E = [0.0, 0.0, 0.0]
            B = [0.0, 0.0, 1.0]
            q2m = 1.0e8
            dt = 1.0e-9

            TestParticleC.testparticle_multistep_boris_step_staggered_f64(
                pointer(x), pointer(v_half), q2m, dt, pointer(E), pointer(B),
                Cint(2), Cint(4)
            )
            @test norm(v_half) ≈ 1.0e5 atol = 1.0e-8
            @test x ≈ v_half * dt

            n = 4
            xb = zeros(3 * n)
            vb = zeros(3 * n)
            Eb = zeros(3 * n)
            Bb = zeros(3 * n)
            for i in 0:(n - 1)
                vb[3 * i + 1] = 1.0e5
                Bb[3 * i + 3] = 1.0
            end

            TestParticleC.testparticle_multistep_boris_push_batch_f64(
                n, pointer(xb), pointer(vb), q2m, dt, pointer(Eb), pointer(Bb),
                2, 4
            )
            for i in 0:(n - 1)
                vi = [vb[3 * i + 1], vb[3 * i + 2], vb[3 * i + 3]]
                @test norm(vi) ≈ 1.0e5 atol = 1.0e-8
            end
        end
    end

    @testset "Adaptive Boris and Resync" begin
        let
            B = [0.0, 0.0, 2.0]
            q2m = 1.0e8
            safety = 0.05
            dt_prop = TestParticleC.testparticle_boris_propose_dt_f64(
                q2m, pointer(B), safety
            )
            @test dt_prop ≈ (2π * 0.05) / (1.0e8 * 2.0)

            # Test Resync
            v_half = [1.0e5, 0.0, 0.0]
            E = [0.0, 0.0, 0.0]
            dt_prev = 1.0e-9
            dt_next = 2.0e-9
            TestParticleC.testparticle_boris_resync_velocity_half_f64(
                pointer(v_half), q2m, dt_prev, dt_next, pointer(E), pointer(B)
            )
            @test norm(v_half) ≈ 1.0e5 atol = 1.0e-8

            # Test Adaptive Step
            x = [0.0, 0.0, 0.0]
            vh = [1.0e5, 0.0, 0.0]
            dt_out = [0.0]
            TestParticleC.testparticle_adaptive_boris_step_f64(
                pointer(x), pointer(vh), dt_prev, pointer(dt_out),
                q2m, pointer(E), pointer(B), safety
            )
            @test norm(vh) ≈ 1.0e5 atol = 1.0e-8
            @test norm(x) > 0.0
            @test dt_out[1] ≈ dt_prop
        end
    end

    @testset "SoA and Multi-Species Batch Push" begin
        let
            n = 5
            x = zeros(n)
            y = zeros(n)
            z = zeros(n)
            vx = fill(1.0e5, n)
            vy = zeros(n)
            vz = zeros(n)
            Ex = zeros(n)
            Ey = zeros(n)
            Ez = zeros(n)
            Bx = zeros(n)
            By = zeros(n)
            Bz = fill(1.0, n)
            q2m_arr = [1.0e8, 2.0e8, 0.5e8, 1.5e8, 1.0e8]
            dt = 1.0e-9

            # Test SoA with array of q2m
            TestParticleC.testparticle_boris_push_batch_soa_f64(
                Csize_t(n),
                pointer(x), pointer(y), pointer(z),
                pointer(vx), pointer(vy), pointer(vz),
                pointer(q2m_arr), 0.0, dt,
                pointer(Ex), pointer(Ey), pointer(Ez),
                pointer(Bx), pointer(By), pointer(Bz)
            )

            for i in 1:n
                vi = [vx[i], vy[i], vz[i]]
                @test norm(vi) ≈ 1.0e5 atol = 1.0e-8
            end
            @test any(x .!= 0.0)

            # Test Staggered SoA with scalar q2m
            TestParticleC.testparticle_boris_step_staggered_batch_soa_f64(
                Csize_t(n),
                pointer(x), pointer(y), pointer(z),
                pointer(vx), pointer(vy), pointer(vz),
                Ptr{Cdouble}(C_NULL), 1.0e8, dt,
                pointer(Ex), pointer(Ey), pointer(Ez),
                pointer(Bx), pointer(By), pointer(Bz)
            )
            for i in 1:n
                vi = [vx[i], vy[i], vz[i]]
                @test norm(vi) ≈ 1.0e5 atol = 1.0e-8
            end

            # Test Multi-species AoS
            x_aos = zeros(3 * n)
            v_aos = zeros(3 * n)
            E_aos = zeros(3 * n)
            B_aos = zeros(3 * n)
            for i in 0:(n - 1)
                v_aos[3 * i + 1] = 1.0e5
                B_aos[3 * i + 3] = 1.0
            end
            TestParticleC.testparticle_boris_push_batch_multispecies_f64(
                Csize_t(n), pointer(x_aos), pointer(v_aos), pointer(q2m_arr),
                dt, pointer(E_aos), pointer(B_aos)
            )
            for i in 0:(n - 1)
                vi = [v_aos[3 * i + 1], v_aos[3 * i + 2], v_aos[3 * i + 3]]
                @test norm(vi) ≈ 1.0e5 atol = 1.0e-8
            end
        end
    end

    @testset "Trajectory Field Callback Tracing" begin
        let
            function dummy_field(t, x_ptr, E_ptr, B_ptr, user_data)
                unsafe_store!(E_ptr, 0.0, 1)
                unsafe_store!(E_ptr, 0.0, 2)
                unsafe_store!(E_ptr, 0.0, 3)
                unsafe_store!(B_ptr, 0.0, 1)
                unsafe_store!(B_ptr, 0.0, 2)
                unsafe_store!(B_ptr, 1.0, 3)
                return nothing
            end

            c_field = @cfunction(
                $dummy_field, Cvoid,
                (Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cvoid})
            )

            x = [0.0, 0.0, 0.0]
            v = [1.0e5, 0.0, 0.0]
            t_start = 0.0
            t_end = 1.0e-8
            dt = 1.0e-9
            q2m = 1.0e8
            max_steps = Csize_t(100)
            n_taken = [Csize_t(0)]
            traj = zeros(7 * 100)

            ret = TestParticleC.testparticle_trace_boris_f64(
                pointer(x), pointer(v), t_start, t_end, dt, q2m,
                c_field, C_NULL, max_steps, pointer(n_taken), pointer(traj)
            )

            @test ret == 0
            @test n_taken[1] == 10
            @test norm(v) ≈ 1.0e5 atol = 1.0e-8
            @test norm(x) > 0.0
            @test traj[1] == 0.0 # initial time
        end
    end

    @testset "Empty Batch Push (n == 0)" begin
        let
            null64 = Ptr{Cdouble}(C_NULL)
            null32 = Ptr{Cfloat}(C_NULL)

            @test isnothing(
                TestParticleC.testparticle_boris_push_batch_f64(
                    Csize_t(0), null64, null64, 1.0e8, 1.0e-9, null64, null64
                )
            )
            @test isnothing(
                TestParticleC.testparticle_boris_push_batch_f32(
                    Csize_t(0), null32, null32, 1.0f8, 1.0f-9, null32, null32
                )
            )
            @test isnothing(
                TestParticleC.testparticle_multistep_boris_push_batch_f64(
                    Csize_t(0), null64, null64, 1.0e8, 1.0e-9, null64, null64,
                    Cint(2), Cint(4)
                )
            )
            @test isnothing(
                TestParticleC.testparticle_multistep_boris_push_batch_f32(
                    Csize_t(0), null32, null32, 1.0f8, 1.0f-9, null32, null32,
                    Cint(2), Cint(4)
                )
            )
            @test isnothing(
                TestParticleC.testparticle_boris_push_batch_multispecies_f64(
                    Csize_t(0), null64, null64, null64, 1.0e-9, null64, null64
                )
            )
        end
    end
end

end # module test_testparticlec
