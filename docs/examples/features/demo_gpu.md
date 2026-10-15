
# GPU Ensemble Tracing

This example demonstrates the usage of GPU for [ensemble tracing](@ref Ensemble-Tracing).
Since GitHub Actions do not have GPU runners for now, we do not show the results on page.

```julia
using TestParticle
using DiffEqGPU, OrdinaryDiffEq, CUDA, StaticArrays
using CairoMakie

"Set initial state for EnsembleProblem."
function prob_func(prob, ctx)
    return prob = @views remake(prob, u0 = [prob.u0[1:3]..., ctx.sim_id / 3, 0.0, 0.0])
end

## Initialization

B(x) = SA[0, 0, 1.0e-11]
E(x) = SA[0, 0, 1.0e-13]

x0 = [0.0, 0.0, 0.0] # initial position, [m]
u0 = [1.0, 0.0, 0.0] # initial velocity, [m/s]
stateinit = [x0..., u0...]

param = prepare(E, B; species = Electron)
tspan = (0.0, 10.0)

trajectories = 3

## Solve for the trajectories

prob = ODEProblem(trace!, stateinit, tspan, param)
ensemble_prob = EnsembleProblem(prob; prob_func, safetycopy = false)
sols = solve(ensemble_prob, Tsit5(), EnsembleGPUArray(CUDA.CUDABackend()); trajectories)

## Visualization

f = Figure(fontsize = 18)
ax = Axis3(
    f[1, 1],
    title = "Electron trajectories",
    xlabel = "X",
    ylabel = "Y",
    zlabel = "Z",
    aspect = :data,
)

for i in 1:trajectories
    lines!(ax, sols.u[i], idxs = (1, 2, 3), label = "$i", color = Makie.wong_colors()[i])
end

f
```

While `EnsembleGPUArray` has a bit of overhead due to its form of GPU code construction, `EnsembleGPUKernel` is a more restrictive GPU-itizing algorithm that achieves a much lower overhead in kernel launching costs. However, it requires this problem to be _written in out-of-place form_ and _use special solvers_. Additionally, a timestep `dt` or `saveat` keyword is required for dense outputs.

```julia
using TestParticle
using DiffEqGPU, OrdinaryDiffEq, CUDA, StaticArrays
using CairoMakie

"Set initial state for EnsembleProblem."
function prob_func(prob, ctx)
    return prob = @views remake(prob, u0 = SA[prob.u0[1:3]..., ctx.sim_id / 3, 0.0, 0.0])
end

## Initialization

B(x) = SA[0, 0, 1.0e-11]
E(x) = SA[0, 0, 1.0e-13]

x0 = SA[0.0, 0.0, 0.0] # initial position, [m]
u0 = SA[1.0, 0.0, 0.0] # initial velocity, [m/s]
stateinit = SA[x0..., u0...]

param = prepare(E, B; species = Electron)
tspan = (0.0, 10.0)

trajectories = 3

## Solve for the trajectories

prob = ODEProblem(trace, stateinit, tspan, param)
ensemble_prob = EnsembleProblem(prob; prob_func, safetycopy = false)
## saving time interval is required for dense output!
sols = solve(
    ensemble_prob, GPUTsit5(), EnsembleGPUKernel(CUDA.CUDABackend());
    trajectories, saveat = 0.4
)

## Visualization

f = Figure(fontsize = 18)
ax = Axis3(
    f[1, 1],
    title = "Electron trajectories",
    xlabel = "X",
    ylabel = "Y",
    zlabel = "Z",
    aspect = :data,
)

for i in 1:trajectories
    lines!(ax, sols.u[i], idxs = (1, 2, 3), label = "$i", color = Makie.wong_colors()[i])
end

f
```

## Native GPU Boris Solver

TestParticle provides a native GPU Boris solver implemented with
[KernelAbstractions.jl](https://github.com/JuliaGPU/KernelAbstractions.jl), which enables
backend-agnostic GPU execution. Each `EnsembleKernel` solve runs on one selected
backend, including one GPU, using method dispatch on the `KA.Backend` type.

Each trajectory is given its own initial state by `prob_func`, which is called once per
trajectory on the host before the kernel is launched. Being host code, it may run anything
Julia can; only the states it produces reach the device. Its `ctx` carries `sim_id`, the
index of the trajectory within the ensemble, which is all that is needed to place the
particles deterministically. Here 1000 protons start from the origin with perpendicular
speeds spanning a factor of ten, hence energies spanning a factor of a hundred:

```julia
using TestParticle, KernelAbstractions, StaticArrays

# Define fields
B(x) = SA[0.0, 0.0, 1.0e-8]  # Uniform B field
E(x) = SA[0.0, 0.0, 0.0]     # No E field

x0 = [0.0, 0.0, 0.0]         # initial position, [m]
v0 = [1.0e5, 0.0, 0.0]       # initial velocity, [m/s]
stateinit = [x0..., v0...]   # template state, overwritten per trajectory
tspan = (0.0, 1.0e-6)

param = prepare(E, B; species = Proton)

# Give every particle a perpendicular speed of its own
trajectories = 1000
speeds = range(1.0e5, 1.0e6; length = trajectories)

function prob_func(prob, ctx)
    return remake(prob; u0 = SA[x0..., speeds[ctx.sim_id], 0.0, 0.0])
end

prob = TraceProblem(stateinit, tspan, param; prob_func)

# Solve on CPU backend
backend = CPU()
sols = solve(prob, Boris(), backend; dt = 1.0e-9, trajectories, saveat = 1.0e-8)
```

`sols.u[i]` is the orbit of the particle started at `speeds[i]`, so the ensemble covers a
band of gyroradii `r_L = m v_perp / |q| B` rather than retracing one orbit a thousand times.

For very large ensembles, where one `remake` per trajectory starts to show, the same states
can be handed over as a `trajectories × 6` matrix: `solve(prob, Boris(), backend; dt,
trajectories, u0 = states, saveat)`. That skips `prob_func` entirely, copying the states to
the device in one go.

The native GPU Boris solver:
- Uses `@kernel` macro from KernelAbstractions.jl for backend-agnostic execution
- Supports multiple GPU backends (CUDA, ROCm, Metal, oneAPI) and CPU fallback
- Dispatches on `KA.Backend` type for GPU execution
- Processes particles in parallel on the GPU
- Returns solutions in the same format as the CPU solver
- Runs every fixed step solver, `Boris()` as well as `MultistepBoris{N}`, sharing
  the same stateless step implementation as the CPU solver, see
  [Boris Pusher](@ref Boris-Pusher)

The adaptive solvers are not available here: choosing a time step happens on the
host, once per step for the whole ensemble, so trace those with `Boris()` on the
CPU instead.

To use actual GPU acceleration, install the appropriate backend package and create the corresponding backend:
```julia
# For NVIDIA GPUs
using CUDA
backend = CUDABackend()
sols = solve(prob, Boris(), backend; dt = 1.0e-9, trajectories, saveat = 1.0e-8)

# For AMD GPUs
using AMDGPU
backend = ROCBackend()
sols = solve(prob, Boris(), backend; dt = 1.0e-9, trajectories, saveat = 1.0e-8)
```

> **Note**: The native GPU solver supports both analytic and numerical (interpolated) fields. Numerical fields see a particularly large performance benefit from GPU acceleration.

### Using multiple GPUs today

The native solver does not yet distribute one ensemble across multiple GPUs
automatically. For independent particle trajectories, launch one Julia worker per
GPU, select a device in each worker, split the initial-state matrix into chunks,
and collect the chunk solutions in order. This example is for multiple GPUs on a
single machine; each worker must use the same Julia project and have CUDA available.

```julia
using Distributed
using TestParticle, CUDA, StaticArrays

trajectories = 1000
gpu_count = length(CUDA.devices())
gpu_count > 0 || error("No CUDA GPUs are available")
trajectories >= gpu_count || error("Use at least one trajectory per GPU")

project = Base.active_project()
isnothing(project) && error("Activate the project used to run this example")
pids = addprocs(gpu_count; exeflags = "--project=$(dirname(project))")

@everywhere using TestParticle, CUDA, StaticArrays

@everywhere begin
    B(x) = SA[0.0, 0.0, 1.0e-8]
    E(x) = SA[0.0, 0.0, 0.0]

    function solve_gpu_chunk(prob, alg, u0_chunk, device_index, dt)
        CUDA.device!(CUDA.devices()[device_index])
        return solve(
            prob, alg, CUDA.CUDABackend();
            dt, trajectories = size(u0_chunk, 1), u0 = u0_chunk
        )
    end
end

stateinit = [0.0, 0.0, 0.0, 1.0e5, 0.0, 0.0]
tspan = (0.0, 1.0e-6)
dt = 1.0e-9

param = prepare(E, B; species = Proton)
prob = TraceProblem(stateinit, tspan, param)
u0 = repeat(reshape(stateinit, 1, 6), trajectories, 1)
for i in 1:trajectories
    u0[i, 4] = 1.0e5 * i / trajectories
end

# Use contiguous, balanced ranges so concatenating results restores input order.
ranges = [
    ((i - 1) * trajectories ÷ gpu_count + 1):(i * trajectories ÷ gpu_count)
    for i in 1:gpu_count
]

try
    futures = [
        remotecall(
            solve_gpu_chunk, pids[i], prob, Boris(), u0[r, :], i, dt
        )
        for (i, r) in enumerate(ranges)
    ]
    solutions = reduce(vcat, (fetch(f).u for f in futures))
finally
    rmprocs(pids)
end
```

Each worker runs the existing single-backend `EnsembleKernel` solve for its
chunk, and the caller concatenates the returned solutions in trajectory order.
This is explicit orchestration rather than a built-in multi-GPU scheduler; it
assumes all workers can see the same GPU list and that there are at least as many
trajectories as GPUs. A future multi-GPU API should manage device assignment and
chunking, retain stable trajectory IDs and seeded results regardless of scheduling,
and gather results consistently for both solution and raw-output modes.

## Boundary Checks with isoutside

The Boris kernel solver supports early termination through the `isoutside` keyword
argument on both CPU and device backends. The boundary predicate `isoutside(u, p, t)`
accepts the current state vector `u`, parameters `p`, and time `t`, returning `true`
when the particle has left the computational domain.

When `isoutside` returns `true`, the solver rolls back the step that exited the domain,
terminates the trace at the last state inside it, and marks the solution with
`ReturnCode.Terminated`.

Because the check executes inside the device kernel, the predicate must be an `isbits`
function with no heap allocations:

```julia
using TestParticle, KernelAbstractions, StaticArrays
using CUDA
using SciMLBase: ReturnCode

# Uniform magnetic field
B(x) = SA[0.0, 0.0, 1.0e-8]
E(x) = SA[0.0, 0.0, 0.0]

x0 = SA[0.0, 0.0, 0.0]
v0 = SA[1.0e5, 0.0, 0.0]
stateinit = SA[x0..., v0...]
tspan = (0.0, 1.0e-5)
dt = 1.0e-7

param = prepare(E, B; species = Proton)

# Ensemble with varying initial velocities along x
trajectories = 8
prob_func(prob, ctx) = remake(
    prob; u0 = SA[0.0, 0.0, 0.0, ctx.sim_id * 2.0e4, 0.0, 0.0]
)
prob = TraceProblem(stateinit, tspan, param; prob_func)

# Terminate when the particle passes x = 0.5 m
isoutside(u, p, t) = u[1] > 0.5

backend = CUDABackend() # or CPU()
sols = solve(prob, Boris(), backend; dt, trajectories, isoutside)

for (i, sol) in enumerate(sols.u)
    println(
        "Trajectory ", i, ": retcode = ", sol.retcode,
        ", final x = ", sol.u[end][1]
    )
end
```

The check is compatible with all output modes, including `saveat`,
`save_everystep = false`, and `raw_output = true`. When `raw_output = true` is used,
early-terminated trajectories repeat their final inside state across any remaining
output slots, and the `terminated` boolean vector indicates which trajectories
stopped early:

```julia
raw = solve(prob, Boris(), backend; dt, trajectories, isoutside, raw_output = true)
# raw.terminated is a Vector{Bool} indicating which trajectories left the domain
```

## Spherical Grids on the Device

Fields stored on a spherical grid are traced on a device without any extra step. Preparing
them with `gridtype = StructuredGrid` builds a spherical interpolator, and handing the
parameters to a device backend is what converts it into a [`GPUSphericalGrid`](@ref), a
trilinear interpolator over `(r, θ, ϕ)` that carries its axes along in device memory. Uniform
grid vectors and non-uniform ones, a logarithmic `r` for instance, are both supported; the
interpolator itself is described in [Field Interpolation](@ref).

```julia
using TestParticle, KernelAbstractions, StaticArrays
using CUDA
import TestParticle as TP

# A uniform 10 nT field along z, stored in spherical components
r = logrange(1.0, 10.0, length = 32)   # non-uniform in r
θ = range(0, π, length = 32)
ϕ = range(0, 2π, length = 32)

B = zeros(3, length(r), length(θ), length(ϕ))
for (iθ, θv) in enumerate(θ)
    sinθ, cosθ = sincos(θv)
    B[1, :, iθ, :] .= 1.0e-8 * cosθ   # Br
    B[2, :, iθ, :] .= -1.0e-8 * sinθ  # Bθ
end

stateinit = [2.0, 2.0, 2.0, 1.0, 0.0, 0.0]  # [m], [m/s]
tspan = (0.0, 1.0)

param = prepare(r, θ, ϕ, ZeroField(), B; species = Proton, gridtype = TP.StructuredGrid)
prob = TraceProblem(stateinit, tspan, param)

# Moving the parameters to the device converts the field into a GPUSphericalGrid
backend = CUDABackend()
sols = solve(prob, Boris(), backend; dt = 1.0e-4, trajectories = 1000, saveat = 1.0e-2)
```

Locations stay Cartesian: the grid converts them to `(r, θ, ϕ)`, interpolates the spherical
components, and rotates a vector result back into the Cartesian basis, so a field stored as
`(Br, Bθ, Bϕ)` comes back as `(Bx, By, Bz)`. Outside the grid the spherical interpolator fills
with `NaN` in `r` and in `θ`, and wraps periodically in `ϕ`.

One thing is worth keeping in mind when comparing a device run against the host: the device
grid always interpolates linearly, whatever `order` the host interpolator was built with, so
keep `order = 1`, the default, when the two are meant to agree.

On the `CPU()` backend the host interpolator is used unchanged and no conversion happens; a
`GPUSphericalGrid` therefore only appears for an actual device backend, or when one is built
explicitly as shown in [Field Interpolation](@ref).
