@testset "Output saving flags" begin
    # Setup
    x0 = [0.0, 0.0, 0.0]
    v0 = [0.0, 1e5, 0.0]
    stateinit = [x0..., v0...]
    tspan = (0.0, 3e-8)
    dt = 3e-11
    zero_E = TP.ZeroField()
    uniform_B2(x) = SA[0.0, 0.0, 0.01]
    param = prepare(zero_E, uniform_B2, species = Electron)
    prob = TraceProblem(stateinit, tspan, param)

    # Baseline: save_everystep=true (default), save_start=true (default implicit), save_end=true (default implicit)
    # savestepinterval = 10
    # nt = 1000. steps = 1000/10 = 100.
    # nout = 101 (0, 10, ..., 1000)
    sol = TP.solve(prob; dt=dt, savestepinterval=10)[1]
    @test length(sol) == 101
    @test sol.t[1] == tspan[1]
    @test sol.t[end] == tspan[2]

    # Scenario 2: Only final state
    sol = TP.solve(prob; dt=dt, savestepinterval=10, save_everystep=false, save_start=false, save_end=true)[1]
    @test length(sol) == 1
    @test sol.t[1] == tspan[2]

    # Scenario 3: Start and End
    sol = TP.solve(prob; dt=dt, savestepinterval=10, save_everystep=false, save_start=true, save_end=true)[1]
    @test length(sol) == 2
    @test sol.t[1] == tspan[1]
    @test sol.t[end] == tspan[2]

    # Scenario 4: Only start
    sol = TP.solve(prob; dt=dt, savestepinterval=10, save_everystep=false, save_start=true, save_end=false)[1]
    @test length(sol) == 1
    @test sol.t[1] == tspan[1]

    # Scenario 5: Every step but no start/end
    sol = TP.solve(prob; dt=dt, savestepinterval=10, save_everystep=true, save_start=false, save_end=false)[1]
    # Steps: 10, 20, ..., 990.
    @test length(sol) == 99
    @test sol.t[1] ≈ tspan[1] + 10*dt
    @test sol.t[end] ≈ tspan[1] + 990*dt

    # Scenario 6: Irregular interval
    sol = TP.solve(prob; dt=dt, savestepinterval=3, save_everystep=true, save_start=true, save_end=true)[1]
    @test length(sol) == 335
    @test sol.t[1] == tspan[1]
    @test sol.t[2] ≈ tspan[1] + 3*dt
    @test sol.t[end-1] ≈ tspan[1] + 999*dt
    @test sol.t[end] == tspan[2]

    # Multistep Boris test
    sol_ms = TP.solve(prob; dt=dt, savestepinterval=10, n=2, save_everystep=false, save_start=true, save_end=true)[1]
    @test length(sol_ms) == 2
    @test sol_ms.t[1] == tspan[1]
    @test sol_ms.t[end] == tspan[2]
end
