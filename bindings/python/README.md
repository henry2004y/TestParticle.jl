# testparticle

Python wrapper for [TestParticle.jl](https://github.com/JuliaSpacePhysics/TestParticle.jl) - trace charged particles in electromagnetic fields via Julia.

The Julia package and its dependencies are installed automatically on first
use through [JuliaCall](https://juliacall.readthedocs.io/) (see `juliapkg.json`).

## Installation

The wrapper drives the Julia package `TestParticle.jl` through
[JuliaCall](https://juliacall.readthedocs.io/) and
[juliapkg](https://github.com/JuliaPy/juliapkg). You do **not** need to install
`TestParticle.jl` yourself: the Julia environment (including `TestParticle` and
its dependencies) is resolved from `juliapkg.json` and downloaded/precompiled
automatically on first use.

### Requirements

- Python 3.10 or later
- Julia 1.10 or later, installed and available on your `PATH`

### From PyPI

```bash
python -m pip install testparticle
```

### From a checkout of this repository

```bash
cd bindings/python
python -m pip install -e .
```

By default the wrapper resolves `TestParticle.jl` from the Julia General
registry. To develop and test the wrapper against the **local** checkout of
this repository, register it as a development dependency (this writes an
environment-local override and does not modify the packaged `juliapkg.json`):

```bash
# from bindings/python
python -m juliapkg add TestParticle --dev --path ..
# if the CLI is unavailable: pip install "juliapkg[cli]"
# or equivalently, without the CLI extra:
python -c "from juliapkg import add; add('TestParticle', uuid='953b605b-f162-4481-8f7f-a191c2bb40e3', dev=True, path='..')"
```

The first time you `import testparticle` (or run the tests) juliapkg downloads
and precompiles the Julia environment. This can take several minutes and only
happens once.

## Usage

```python
import testparticle as tp
import numpy as np

# Uniform magnetic field of 10 nT along z.
B = 1e-8

def get_B(x):
    return [0.0, 0.0, B]

# Initial position [m] and velocity [m/s].
x0 = [0.0, 0.0, 0.0]
v0 = [1e5, 0.0, 0.0]

# Trace a proton for 0.1 us with a 0.5 ns step.
res = tp.trace(x0, v0, (0.0, 1e-7), 5e-10, get_B, species="proton")

# `res` exposes numpy arrays.
t = res.t        # time stamps, shape (n,)
x = res.x        # positions, shape (n, 3)
v = res.v        # velocities, shape (n, 3)

# Save the EM fields along the trajectory as well.
res = tp.trace(x0, v0, (0.0, 1e-7), 5e-10, get_B, species="proton", save_fields=True)
E = res.E        # shape (n, 3)
B_sol = res.B    # shape (n, 3)

# Guiding-center tracing.
gc = tp.trace_gc(x0, v0, (0.0, 1e-7), 5e-10, get_B, species="proton")
R = gc.R         # guiding-center position, shape (n, 3)
vpar = gc.vpar   # parallel velocity, shape (n,)

# Magnetic field-line tracing (arc-length parameter).
fl = tp.trace_fieldline([1.0, 0.0, 0.0], get_B, (0.0, 1.0), mode="both")
t_fwd, x_fwd = fl["forward"]
t_bwd, x_bwd = fl["backward"]
```

Field functions may be either time-independent `f(x)` or time-dependent
`f(x, t)`, and must return a length-3 sequence. Alternatively you can pass a
Julia field object directly (e.g. a field built with the `Magnetostatics`
package).

## API

The high-level helpers are:

- `trace(x0, v0, tspan, dt, B, E=None, *, species="proton", trajectories=1,
  savestepinterval=1, save_fields=False, save_work=False, parallel=False,
  alg=None, init_func=None, q=None, m=None)` - Boris pusher.
- `trace_gc(x0, v0, tspan, dt, B, E=None, *, species="proton", ...,
  alg="rk4")` - guiding-center RK4/RK45 solver.
- `trace_fieldline(x0, B, tspan, *, mode="both", alg="Tsit5",
  abstol=1e-8, reltol=1e-8)` - magnetic field lines.

Raw Julia objects (`prepare`, `TraceProblem`, `solve`, `Proton`, `get_gyrofrequency`,
...) are re-exported under their Julia names; functions ending with `!` in Julia
are available with the `_b` suffix (e.g. `trace!` -> `trace_b`).

See the [TestParticle.jl documentation](https://henry2004y.github.io/TestParticle.jl/dev/)
for the full description of the underlying Julia functions.
