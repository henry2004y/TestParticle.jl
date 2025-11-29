
# GPU Ensemble tracing

This example demonstrates the usage of GPU for [ensemble tracing](@ref Ensemble-Tracing).
Since GitHub Actions do not have GPU runners for now, we do not show the results on page.

```julia
using TestParticle
using DiffEqGPU, OrdinaryDiffEq, CUDA, StaticArrays
using CairoMakie

"Set initial state for EnsembleProblem."
function prob_func(prob, i, repeat)
   prob = @views remake(prob, u0=[prob.u0[1:3]..., i/3, 0.0, 0.0])
end

## Initialization

B(x) = SA[0, 0, 1e-11]
E(x) = SA[0, 0, 1e-13]

x0 = [0.0, 0.0, 0.0] # initial position, [m]
u0 = [1.0, 0.0, 0.0] # initial velocity, [m/s]
stateinit = [x0..., u0...]

param = prepare(E, B, species=Electron)
tspan = (0.0, 10.0)

trajectories = 3

## Solve for the trajectories

prob = ODEProblem(trace!, stateinit, tspan, param)
ensemble_prob = EnsembleProblem(prob; prob_func, safetycopy=false)
sols = solve(ensemble_prob, Tsit5(), EnsembleGPUArray(CUDA.CUDABackend()); trajectories)

## Visualization

f = Figure(fontsize = 18)
ax = Axis3(f[1, 1],
   title = "Electron trajectories",
   xlabel = "X",
   ylabel = "Y",
   zlabel = "Z",
   aspect = :data,
)

for i in eachindex(sols)
   lines!(ax, sols[i], idxs=(1,2,3), label="$i", color=Makie.wong_colors()[i])
end

f
```

While `EnsembleGPUArray` has a bit of overhead due to its form of GPU code construction, `EnsembleGPUKernel` is a more restrictive GPU-itizing algorithm that achieves a much lower overhead in kernel launching costs. However, it requires this problem to be _written in out-of-place form_ and _use special solvers_. Additionally, a timestep `dt` or `saveat` keyword is required for dense outputs.

```julia
using TestParticle
using DiffEqGPU, OrdinaryDiffEq, CUDA, StaticArrays
using CairoMakie

"Set initial state for EnsembleProblem."
function prob_func(prob, i, repeat)
   prob = @views remake(prob, u0=SA[prob.u0[1:3]..., i/3, 0.0, 0.0])
end

## Initialization

B(x) = SA[0, 0, 1e-11]
E(x) = SA[0, 0, 1e-13]

x0 = SA[0.0, 0.0, 0.0] # initial position, [m]
u0 = SA[1.0, 0.0, 0.0] # initial velocity, [m/s]
stateinit = SA[x0..., u0...]

param = prepare(E, B, species=Electron)
tspan = (0.0, 10.0)

trajectories = 3

## Solve for the trajectories

prob = ODEProblem(trace, stateinit, tspan, param)
ensemble_prob = EnsembleProblem(prob; prob_func, safetycopy=false)
## saving time interval is required for dense output! 
sols = solve(ensemble_prob, GPUTsit5(), EnsembleGPUKernel(CUDA.CUDABackend());
   trajectories, saveat=0.4)

## Visualization

f = Figure(fontsize = 18)
ax = Axis3(f[1, 1],
   title = "Electron trajectories",
   xlabel = "X",
   ylabel = "Y",
   zlabel = "Z",
   aspect = :data,
)

for i in eachindex(sols)
   lines!(ax, sols[i], idxs=(1,2,3), label="$i", color=Makie.wong_colors()[i])
end

f
```
