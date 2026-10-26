using Test
using TestParticleCLI
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

    let ret = TestParticleCLI.main([
            example_uniform,
            "--tspan1", "3.0e-9",
            "--dt", "3.0e-11",
            "--saveat", "3.0e-10",
            "--output", out_csv,
        ])
        @test ret == 0
        @test isfile(out_csv)
        lines = readlines(out_csv)
        @test length(lines) == 12
        @test lines[1] == "t,x,y,z,vx,vy,vz"
        rm(out_csv; force = true)
    end

    let ret = TestParticleCLI.main([
            example_dipole,
            "--tspan1", "0.05",
            "--dt", "0.001",
            "--saveat", "0.01",
            "--output", out_csv,
        ])
        @test ret == 0
        @test isfile(out_csv)
        lines = readlines(out_csv)
        @test length(lines) == 7
        rm(out_csv; force = true)
    end

    let ret = TestParticleCLI.main([
            example_uniform,
            "--tspan1", "3.0e-9",
            "--dt", "3.0e-11",
            "--saveat", "3.0e-10",
            "--output", out_h5,
        ])
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

    let ret = TestParticleCLI.main([
            example_uniform,
            "--tspan1", "3.0e-9",
            "--dt", "3.0e-11",
            "--saveat", "3.0e-10",
            "--output", out_jld2,
        ])
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
