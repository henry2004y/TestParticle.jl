# TestParticleCLI

`TestParticleCLI` is a standalone, statically compilable command-line interface for particle tracing using [TestParticle.jl](../../).

## Features

- **Config-Driven Runs**: Configure initial conditions, simulation parameters, and electromagnetic fields via TOML.
- **CLI Overrides**: Override parameters like `--dt`, `--saveat`, `--output`, and `--tspan1` directly from the terminal.
- **Field Catalog**: Supports `uniform`, `dipole`, and `harris` current sheet fields.
- **Fast Trajectory I/O**: Streams trajectory data to plain-text CSV (`.csv`) or HDF5-compatible binary containers (`.h5`, `.hdf5`, `.jld2`) using JLD2 (stores `t`, `x`, `y`, `z`, `vx`, `vy`, `vz`, and `states`).
- **Statically Compilable**: Compatible with `JuliaC` to compile into a native standalone binary.

## Direct Execution via Julia

```bash
julia --project=apps/TestParticleCLI -m TestParticleCLI config.toml [options]
```

Example:
```bash
julia --project=apps/TestParticleCLI -m TestParticleCLI apps/TestParticleCLI/examples/uniform_electron.toml --output trajectory.csv
```

## Static Compilation with JuliaC

To build a standalone executable:

```bash
juliac \
  --output-exe tp_cli \
  --bundle apps/TestParticleCLI/build \
  --trim=no \
  apps/TestParticleCLI
```

The compiled binary will be placed at `apps/TestParticleCLI/build/bin/tp_cli` (or `.exe` on Windows) and can be executed without a Julia installation:

```bash
./apps/TestParticleCLI/build/bin/tp_cli config.toml --output trajectory.csv
```

## Configuration Schema

```toml
[simulation]
tspan = [0.0, 1.0e-7]       # [t_start, t_end]
dt = 1.0e-10                # Simulation time step
saveat = 1.0e-9             # Output interval (optional)
output = "trajectory.csv"   # Destination CSV path (optional)

[particle]
species = "electron"        # "electron", "proton", or custom
position = [0.0, 0.0, 0.0]  # Initial [x, y, z]
velocity = [0.0, 1.0e5, 0.0]# Initial [vx, vy, vz]
# Alternatively, specify energy and pitch angle:
# energy_eV = 1000.0
# pitch_angle_deg = 45.0

[field]
type = "uniform"            # "uniform", "dipole", or "harris"
B = [0.0, 0.0, 0.01]        # Magnetic field [Bx, By, Bz]
E = [0.0, 0.0, 0.0]         # Electric field [Ex, Ey, Ez]
```
