# App

`TestParticleCLI` is a standalone, statically compilable command-line interface for
tracing charged particles in electromagnetic fields using [TestParticle.jl](index.md).

It is designed for high-performance computing (HPC) batch jobs, parameter sweeps,
and automated pipelines where you want to execute particle simulations directly from
the terminal or as a compiled native binary without interactive Julia sessions.

The source code is located in [`apps/TestParticleCLI`](https://github.com/henry2004y/TestParticle.jl/tree/master/apps/TestParticleCLI).

---

## Features

- **TOML-Driven Configuration**: Configure particles, fields, and simulation
  parameters declaratively.
- **Terminal Overrides**: Override runtime parameters (such as `dt`, `saveat`, or
  `output`) directly via CLI arguments.
- **Analytic Field Models**: Built-in support for uniform fields, planetary dipole
  fields, and Harris current sheets.
- **Dual Output Formats**: Stream trajectories to standard CSV (`.csv`) or
  HDF5-compatible binary containers (`.h5`, `.hdf5`, `.jld2`) via JLD2.
- **Static Compilation**: Compiles into a standalone native binary with `JuliaC`.

---

## Quick Start

### 1. Direct Execution via Julia

You can run `TestParticleCLI` directly using Julia without compiling:

```bash
# Display help and options
julia --project=apps/TestParticleCLI -m TestParticleCLI --help

# Run simulation from a TOML configuration
julia --project=apps/TestParticleCLI -m TestParticleCLI apps/TestParticleCLI/examples/uniform_electron.toml
```

### 2. Static Compilation with JuliaC

With Julia 1.12+, you can compile `TestParticleCLI` into a native standalone executable
using [JuliaC](https://github.com/JuliaLang/JuliaC.jl):

```bash
juliac \
  --output-exe tp_cli \
  --bundle apps/TestParticleCLI/build \
  --trim=no \
  apps/TestParticleCLI
```

The compiled binary will be placed at `apps/TestParticleCLI/build/bin/tp_cli`
(or `tp_cli.exe` on Windows). It is self-contained and runs without needing Julia
installed:

```bash
./apps/TestParticleCLI/build/bin/tp_cli config.toml --output trajectory.h5
```

---

## Configuration Schema

A simulation is defined in a TOML file containing three main sections:

```toml
[simulation]
tspan = [0.0, 3.0e-8]         # [start_time, end_time] in seconds
dt = 3.0e-11                  # Time step
saveat = 3.0e-10              # Output recording interval (optional)
output = "trajectory.h5"      # Output destination (optional)

[particle]
species = "electron"          # "electron", "proton", or custom
position = [0.0, 0.0, 0.0]    # Initial position [x, y, z] in meters
velocity = [0.0, 1.0e5, 0.0]  # Initial velocity [vx, vy, vz] in m/s

# Alternatively, specify energy and pitch angle instead of velocity:
# energy_eV = 1000.0          # Kinetic energy in eV
# pitch_angle_deg = 45.0      # Pitch angle relative to B in degrees

[field]
type = "uniform"              # "uniform", "dipole", or "harris"
B = [0.0, 0.0, 0.01]          # Magnetic field [Bx, By, Bz] in Tesla
E = [0.0, 0.0, 0.0]           # Electric field [Ex, Ey, Ez] in V/m
```

---

## Supported Field Types

### Uniform Field (`type = "uniform"`)
A constant background electric and magnetic field:
```toml
[field]
type = "uniform"
B = [0.0, 0.0, 0.01]
E = [0.0, 0.0, 0.0]
```

### Dipole Field (`type = "dipole"`)
Earth-centered magnetic dipole field $\mathbf{B}(\mathbf{r})$ with magnetic moment $M$:
```toml
[field]
type = "dipole"
M = 8.0e22                    # Magnetic moment in A·m² (default: Earth)
E = [0.0, 0.0, 0.0]
```

### Harris Current Sheet (`type = "harris"`)
Magnetic field with reversal across a current sheet: $B_x(z) = B_0 \tanh(z/L)$,
guide field $B_y = B_g$:
```toml
[field]
type = "harris"
B0 = 0.01                     # Asymptotic field in Tesla
L = 1.0e5                     # Current sheet half-thickness in meters
Bg = 0.0                      # Guide field in Tesla
E = [0.0, 0.0, 0.0]
```

---

## CLI Options and Overrides

You can override values specified in the TOML configuration file on the command line:

```text
TestParticleCLI - Command Line Particle Tracer
Usage: tp_cli <config.toml> [options]

Options:
  --dt <val>       Override time step dt
  --saveat <val>   Override output save interval
  --output <path>  Write trajectory to file (.csv, .h5, .jld2)
  --tspan1 <val>   Override simulation end time
  -h, --help       Show help message
```

Example overriding time step and output filename:
```bash
tp_cli dipole.toml --dt 1e-4 --output run1.csv
```

---

## Output Formats

`TestParticleCLI` detects the output format automatically based on the file extension:

### Plain-Text CSV (`.csv`)
Produces a standard comma-separated table:
```csv
t,x,y,z,vx,vy,vz
0.0,0.0,0.0,0.0,0.0,100000.0,0.0
3.0e-10,-7.73e-6,2.86e-5,0.0,-50339.5,86405.6,0.0
...
```

### HDF5 Binary (`.h5`, `.hdf5`, `.jld2`)
Uses `JLD2` to produce an HDF5-compatible binary file containing the following datasets:
* **`t`**: 1D array of timestamps $[t_0, t_1, \dots, t_N]$.
* **`x`**, **`y`**, **`z`**: 1D arrays of particle coordinates in meters.
* **`vx`**, **`vy`**, **`vz`**: 1D arrays of particle velocities in m/s.
* **`states`**: Full state matrix of shape `(6, N)`.

You can inspect and read this file directly in Python (`h5py`) or Julia:

```julia
using JLD2

jldopen("trajectory.h5") do f
    t = f["t"]
    states = f["states"] # 6 x N matrix
end
```
