using Test
using TestParticleCLI
using LinearAlgebra: dot, norm
using JLD2: jldopen

@testset "TestParticleCLI" begin
    example_uniform = joinpath(@__DIR__, "..", "examples", "uniform_electron.toml")
    example_dipole = joinpath(@__DIR__, "..", "examples", "dipole_proton.toml")
    out_csv = joinpath(@__DIR__, "test_out.csv")
    out_h5 = joinpath(@__DIR__, "test_out.h5")
    out_jld2 = joinpath(@__DIR__, "test_out.jld2")

    @test TestParticleCLI.main(String[]) == 0
    @test TestParticleCLI.main(["-h"]) == 0
    @test TestParticleCLI.main(["--help"]) == 0
    @test TestParticleCLI.main(["nonexistent.toml"]) == 1

    @testset "Pitch-angle initialization" begin
        field_cfg = TestParticleCLI.parse_field_config(
            Dict{String, Any}("type" => "uniform", "B" => [0.01, 0.0, 0.0])
        )
        _, B_field = TestParticleCLI.make_field_functions(field_cfg)
        particle_cfg = TestParticleCLI.parse_particle_config(
            Dict{String, Any}("energy_eV" => 1000.0, "pitch_angle_deg" => 45.0),
            B_field,
            0.0,
        )
        velocity = particle_cfg.stateinit[4:6]
        magnetic_field = B_field(particle_cfg.stateinit[1:3], 0.0)

        @test dot(velocity, magnetic_field) / (norm(velocity) * norm(magnetic_field)) ≈
            cosd(45.0)
        @test norm(velocity) ≈ sqrt(2.0 * 1000.0 * TestParticleCLI.TP.eV / particle_cfg.m)

        zero_field_cfg = TestParticleCLI.parse_field_config(
            Dict{String, Any}("type" => "uniform", "B" => [0.0, 0.0, 0.0])
        )
        _, zero_B_field = TestParticleCLI.make_field_functions(zero_field_cfg)
        @test_throws ArgumentError TestParticleCLI.parse_particle_config(
            Dict{String, Any}("energy_eV" => 1000.0),
            zero_B_field,
            0.0,
        )
    end

    let ret = TestParticleCLI.main(
            [
                example_uniform,
                "--tspan1", "3.0e-9",
                "--dt", "3.0e-11",
                "--saveat", "3.0e-10",
                "--output", out_csv,
            ]
        )
        @test ret == 0
        @test isfile(out_csv)
        lines = readlines(out_csv)
        @test length(lines) == 12
        @test lines[1] == "t,x,y,z,vx,vy,vz"
        rm(out_csv; force = true)
    end

    let ret = TestParticleCLI.main(
            [
                example_dipole,
                "--tspan1", "0.05",
                "--dt", "0.001",
                "--saveat", "0.01",
                "--output", out_csv,
            ]
        )
        @test ret == 0
        @test isfile(out_csv)
        lines = readlines(out_csv)
        @test length(lines) == 7
        rm(out_csv; force = true)
    end

    let ret = TestParticleCLI.main(
            [
                example_uniform,
                "--tspan1", "3.0e-9",
                "--dt", "3.0e-11",
                "--saveat", "3.0e-10",
                "--output", out_h5,
            ]
        )
        @test ret == 0
        @test isfile(out_h5)
        jldopen(out_h5) do f
            @test haskey(f, "t")
            @test haskey(f, "x")
            @test haskey(f, "states")
            @test length(f["t"]) == 11
            @test size(f["states"]) == (6, 11)
            @test f["x"] == f["states"][1, :]
        end
        rm(out_h5; force = true)
    end

    let ret = TestParticleCLI.main(
            [
                example_uniform,
                "--tspan1", "3.0e-9",
                "--dt", "3.0e-11",
                "--saveat", "3.0e-10",
                "--output", out_jld2,
            ]
        )
        @test ret == 0
        @test isfile(out_jld2)
        jldopen(out_jld2) do f
            @test haskey(f, "t")
            @test haskey(f, "states")
            @test length(f["t"]) == 11
        end
        rm(out_jld2; force = true)
    end
end
