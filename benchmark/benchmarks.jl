# Benchmark Test for TestParticle

using BenchmarkTools
using TestParticle
using OrdinaryDiffEq, Meshes, StaticArrays

const SUITE = BenchmarkGroup()

SUITE["trace"] = BenchmarkGroup()
SUITE["trace"]["analytic field"] = BenchmarkGroup()
SUITE["trace"]["numerical field"] = BenchmarkGroup()
SUITE["trace"]["time-dependent field"] = BenchmarkGroup()

# analytic field
E_analytic(xu) = SA[5e-10, 5e-10, 0]
B_analytic(xu) = SA[0, 0, 1e-8]

# numerical field
x = range(-10, 10, length=15)
y = range(-10, 10, length=20)
z = range(-10, 10, length=25)
B_numeric = fill(0.0, 3, length(x), length(y), length(z)) # [T]
E_numeric = fill(0.0, 3, length(x), length(y), length(z)) # [V/m]
B_numeric[3,:,:,:] .= 10e-9
E_numeric[1,:,:,:] .= 5e-10
E_numeric[2,:,:,:] .= 5e-10
B_td(xu, t) = SA[0, 0, 1e-11*cos(2π*t)]
E_td(xu, t) = SA[5e-11*sin(2π*t), 0, 0]
F_td(xu) = SA[0, 9.10938356e-42, 0]

Δx = x[2] - x[1]
Δy = y[2] - y[1]
Δz = z[2] - z[1]

mesh = CartesianGrid((length(x)-1, length(y)-1, length(z)-1),
    (x[1], y[1], z[1]),
    (Δx, Δy, Δz))

tspan = (0.0, 1.0)
x0 = [0.0, 0.0, 0.0] # initial position, [m]
u0 = [1.0, 0.0, 0.0] # initial velocity, [m/s]
stateinit = [x0..., u0...]

param_analytic = prepare(E_analytic, B_analytic)
prob_ip = ODEProblem(trace!, stateinit, tspan, param_analytic) # in place
prob_oop = ODEProblem(trace, SA[stateinit...], tspan, param_analytic) # out of place
SUITE["trace"]["analytic field"]["in place"] = @benchmarkable solve($prob_ip, Tsit5(); save_idxs=[1,2,3])
SUITE["trace"]["analytic field"]["out of place"] = @benchmarkable solve($prob_oop, Tsit5(); save_idxs=[1,2,3])

param_numeric = prepare(mesh, E_numeric, B_numeric)
prob_ip = ODEProblem(trace!, stateinit, tspan, param_numeric) # in place
prob_oop = ODEProblem(trace, SA[stateinit...], tspan, param_numeric) # out of place
prob_boris = let dt = 1/7, paramBoris = BorisMethod(param_numeric)
    TraceProblem(stateinit, tspan, dt, paramBoris)
end
SUITE["trace"]["numerical field"]["in place"] = @benchmarkable solve($prob_ip, Tsit5(); save_idxs=[1,2,3])
SUITE["trace"]["numerical field"]["out of place"] = @benchmarkable solve($prob_oop, Tsit5(); save_idxs=[1,2,3])
SUITE["trace"]["numerical field"]["Boris"] = @benchmarkable trace_trajectory($prob_boris; savestepinterval=10)

param_td = prepare(E_td, B_td, F_td)
prob_ip = ODEProblem(trace!, stateinit, tspan, param_td) # in place
prob_oop = ODEProblem(trace, SA[stateinit...], tspan, param_td) # out of place
SUITE["trace"]["time-dependent field"]["in place"] = @benchmarkable solve($prob_ip, Tsit5(); save_idxs=[1,2,3])
SUITE["trace"]["time-dependent field"]["out of place"] = @benchmarkable solve($prob_oop, Tsit5(); save_idxs=[1,2,3])