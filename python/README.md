# testparticle-jl

Python wrapper for [TestParticle.jl](https://github.com/JuliaSpacePhysics/TestParticle.jl) - trace charged particles in electromagnetic fields via Julia.

## Installation

```bash
pip install testparticle-jl
```

The Julia package installs automatically on first use.

## Usage

```python
import testparticle as tp
from juliacall import Main as jl
import numpy as np

# Initialize a proton
proton = tp.Proton

# Define E and B fields (simple uniform fields for demonstration)
def get_B(x):
    return [0.0, 0.0, 1e-8] # 10 nT in z-direction

def get_E(x):
    return [0.0, 0.0, 0.0]

# Prepare the trace problem
# Note: Complex setup might require using Julia types directly via `jl` or helper functions
# This is a placeholder for actual usage logic dependent on TestParticle.jl API
```

See the [TestParticle.jl documentation](https://juliaspacephysics.github.io/TestParticle.jl/dev/) for full details on the underlying Julia functions.

## API

Available functions mirror the Julia package.
Functions ending with `!` in Julia are available with `_b` suffix (e.g., `trace!` -> `trace_b`).
