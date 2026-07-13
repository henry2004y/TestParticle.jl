# Python Wrapper

TestParticle.jl ships with a Python package, `testparticle`, that lets you trace
charged particles in electromagnetic fields directly from Python. It drives the
Julia package through [JuliaCall](https://juliacall.readthedocs.io/) and
[juliapkg](https://github.com/JuliaPy/juliapkg), so you do **not** need to
install `TestParticle.jl` yourself: the Julia environment (including
`TestParticle` and its dependencies) is resolved and precompiled automatically
on first use.

The source lives in [`bindings/python`](https://github.com/henry2004y/TestParticle.jl/tree/master/bindings/python).

## Installation

### Requirements

- Python 3.10 or later
- Julia 1.11 or later, installed and available on your `PATH`

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
registry. To develop and test the wrapper against the **local** checkout of this
repository, register it as a development dependency (this writes an
environment-local override and does not modify the packaged `juliapkg.json`):

```bash
# from bindings/python
python -m juliapkg add TestParticle --dev --path ..
# if the CLI is unavailable: pip install "juliapkg[cli]"
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
Julia field object directly.

## High-level API

- `trace(x0, v0, tspan, dt, B, E=None, *, species="proton", trajectories=1,
  savestepinterval=1, save_fields=False, save_work=False, parallel=False,
  alg=None, init_func=None, q=None, m=None)` — full-orbit Boris pusher. Returns a
  `TraceResult` (or a list of them when tracing multiple trajectories).
- `trace_gc(x0, v0, tspan, dt, B, E=None, *, species="proton", ...,
  alg="rk4")` — guiding-center RK4/RK45 solver. Returns a `GCTraceResult`.
- `trace_fieldline(x0, B, tspan, *, mode="both", alg="Tsit5",
  abstol=1e-8, reltol=1e-8)` — magnetic field lines. Returns a dict with
  `"forward"` and/or `"backward"` entries, each a `(t, x)` tuple.

The `species` argument accepts the strings `"proton"`, `"electron"`, or `"ion"`,
a `tp.Ion(m, q)`, or any Julia `Species` object.

### Result objects

`TraceResult` exposes the trajectory as numpy arrays via properties: `t`, `x`,
`v`, and — when the corresponding keyword is enabled — `E`, `B` (`save_fields=True`)
and `work` (`save_work=True`). `GCTraceResult` exposes `t`, `R` (guiding-center
position), `vpar` (parallel velocity), and the same optional `E`, `B`, `work`.

## Raw Julia API

Raw Julia objects (`prepare`, `TraceProblem`, `solve`, `get_gyrofrequency`, ...)
are re-exported under their Julia names. Functions ending with `!` in Julia are
available with the `_b` suffix (e.g. `trace!` → `trace_b`). The active Julia
runtime is accessible as `tp.jl` for advanced use.

See the rest of this documentation for a full description of the underlying
Julia functions.
