# Serial Performance Benchmark for TestParticle

using BenchmarkTools
using TestParticle
using OrdinaryDiffEq
using StaticArrays
using Printf

# Setup uniform field
uniform_B(x) = SA[0.0, 0.0, 1.0e-8]
uniform_E(x) = SA[0.0, 0.0, 0.0]

# Proton parameters by default
param = prepare(uniform_E, uniform_B; species = Proton)

# Initial conditions
x0 = [0.0, 0.0, 0.0]
v0 = [1.0e5, 0.0, 0.0]
stateinit = [x0..., v0...]

# Simulate for 0.1 second with a 1 nanosecond time step
# This represents 100,000 steps.
tspan = (0.0, 0.1)
dt = 1.0e-9

println("="^70)
println("Serial Performance Benchmark: Boris vs. OrdinaryDiffEq")
println("Simulating 1 particle for 0.1 second with dt = 1 ns (10^5 steps).")
println("="^70)

# 1. TestParticle's Custom Boris Solver
println("\nSetting up TestParticle TraceProblem (Boris)...")
prob_boris = TraceProblem(stateinit, tspan, param)

# Precompile
TestParticle.solve(prob_boris; dt = dt, savestepinterval = 10000, maxiters = 200_000_000)

println("Running Boris Benchmark...")
bench_boris = @benchmark TestParticle.solve($prob_boris; dt = $dt, savestepinterval = 100000, maxiters = 200_000_000)
display(bench_boris)

# 2. OrdinaryDiffEq's Tsit5 (Runge-Kutta 5/4)
println("\nSetting up ODEProblem (OrdinaryDiffEq Tsit5)...")
prob_tsit5 = ODEProblem(trace!, stateinit, tspan, param)

# Precompile (limit maxiters for tests to prevent hanging just in case, though it should be fine)
OrdinaryDiffEq.solve(prob_tsit5, Tsit5(); dt = dt, adaptive = false, save_everystep = false)

println("Running Tsit5 Benchmark (fixed dt)...")
bench_tsit5 = @benchmark OrdinaryDiffEq.solve($prob_tsit5, Tsit5(); dt = $dt, adaptive = false, save_everystep = false)
display(bench_tsit5)

# 3. OrdinaryDiffEq's Vern9 (Runge-Kutta 9/8)
println("\nSetting up ODEProblem (OrdinaryDiffEq Vern9)...")
prob_vern9 = ODEProblem(trace!, stateinit, tspan, param)

# Precompile
OrdinaryDiffEq.solve(prob_vern9, Vern9(); dt = dt, adaptive = false, save_everystep = false)

println("Running Vern9 Benchmark (fixed dt)...")
bench_vern9 = @benchmark OrdinaryDiffEq.solve($prob_vern9, Vern9(); dt = $dt, adaptive = false, save_everystep = false)
display(bench_vern9)


# --- Summary Comparison ---
println("\n" * "="^70)
println("Summary Results (Median Times):")

time_boris = median(bench_boris).time / 1.0e6 # in ms
time_tsit5 = median(bench_tsit5).time / 1.0e6 # in ms
time_vern9 = median(bench_vern9).time / 1.0e6 # in ms

@printf("TestParticle Boris: %8.2f ms\n", time_boris)
@printf("ODE Tsit5 (fixed):  %8.2f ms  [ %.2fx slower ]\n", time_tsit5, time_tsit5 / time_boris)
@printf("ODE Vern9 (fixed):  %8.2f ms  [ %.2fx slower ]\n", time_vern9, time_vern9 / time_boris)
println("="^70)
