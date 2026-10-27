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
end

end # module test_testparticlec
