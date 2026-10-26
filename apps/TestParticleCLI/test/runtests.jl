using Test
using TestParticleCLI

@testset "TestParticleCLI" begin
    example_uniform = joinpath(@__DIR__, "..", "examples", "uniform_electron.toml")
    example_dipole = joinpath(@__DIR__, "..", "examples", "dipole_proton.toml")
    out_csv = joinpath(@__DIR__, "test_out.csv")

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
end
