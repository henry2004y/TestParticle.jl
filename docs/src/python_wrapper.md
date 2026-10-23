# Python Wrapper

TestParticle.jl ships with a Python package, `testparticle`, that lets you trace
charged particles in electromagnetic fields directly from Python. It drives the
Julia packages through [JuliaCall](https://juliacall.readthedocs.io/) and
[juliapkg](https://github.com/JuliaPy/juliapkg), so you do **not** need to
install `TestParticle.jl` yourself: the Julia environment (including
`TestParticle`, its companion package `BorisPushers`, and their dependencies) is
resolved and precompiled automatically on first use.

The wrapper loads both packages, so the whole Boris family (standard,
multistep/Hyper, and adaptive) and the adaptive hybrid solver are available
without writing any Julia code; see [`alg`](#solver-selection) below.

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

By default the wrapper resolves `TestParticle.jl` and `BorisPushers.jl` from the
Julia General registry. To develop and test the wrapper against the **local**
checkout of this repository, register both as development dependencies (this
writes an environment-local override and does not modify the packaged
`juliapkg.json`):

```bash
# from bindings/python
python -m juliapkg add TestParticle --dev --path ..
python -m juliapkg add BorisPushers --dev --path ../../lib/BorisPushers
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

```

Field functions may be either time-independent `f(x)` or time-dependent
`f(x, t)`, and must return a length-3 sequence. Alternatively you can pass a
Julia field object directly.

## High-level API

- `trace(x0, v0, tspan, dt, B, E=None, *, species="proton", trajectories=1,
  saveat=None, save_fields=False, save_work=False, parallel=False,
  alg=None, init_func=None, q=None, m=None)` — full-orbit Boris pusher. Returns a
  `TraceResult` (or a list of them when tracing multiple trajectories). The
  `alg` keyword selects the pusher, see [Solver selection](#solver-selection).
- `trace_gc(x0, v0, tspan, dt, B, E=None, *, species="proton", ...,
  saveat=None, alg="rk4")` — guiding-center RK4/RK45 solver. Returns a `GCTraceResult`.

The `species` argument accepts the strings `"proton"`, `"electron"`, or `"ion"`,
a `tp.Ion(m, q)`, or any Julia `Species` object.

### Solver selection

`trace` takes an `alg` keyword that chooses the full-orbit pusher and defaults to
the standard Boris method. Either give the name of a solver from
`BorisPushers.jl`, or pass a constructed algorithm object:

```python
# By name.
res = tp.trace(x0, v0, tspan, dt, get_B, alg="MultistepBoris6")

# As a Julia algorithm object.
res = tp.trace(x0, v0, tspan, dt, get_B, alg=tp.MultistepBoris6(n=4))
```

The names that resolve directly are `"Boris"`, `"AdaptiveBoris"`,
`"MultistepBoris2"`, `"MultistepBoris4"`, and `"MultistepBoris6"`. Algorithms
whose order is a type parameter without an alias, namely `MultistepBoris{N}` for
other `N` and `AdaptiveMultistepBoris{N}`, are built on the Julia side and then
passed as objects:

```python
alg = tp.jl.seval("AdaptiveMultistepBoris{6}(n=4, safety=0.1)")
res = tp.trace(x0, v0, tspan, dt, get_B, alg=alg)
```

The adaptive hybrid solver needs a `TraceHybridProblem` rather than a
`TraceProblem`, so it is not reachable through `trace`; build it with the raw
Julia API instead, as in the hybrid [examples](@ref).

### Result objects

`TraceResult` exposes the trajectory as numpy arrays via properties: `t`, `x`,
`v`, and — when the corresponding keyword is enabled — `E`, `B` (`save_fields=True`)
and `work` (`save_work=True`). `GCTraceResult` exposes `t`, `R` (guiding-center
position), `vpar` (parallel velocity), and the same optional `E`, `B`, `work`.

## Raw Julia API

Raw Julia objects (`prepare`, `TraceProblem`, `solve`, `get_gyrofrequency`, ...)
are re-exported under their Julia names. A name is looked up in `TestParticle`
first and in `BorisPushers` second, so the pusher algorithms (`tp.Boris`,
`tp.MultistepBoris6`, `tp.AdaptiveHybrid`, ...) are available even though they
live in the companion package. Functions ending with `!` in Julia are available
with the `_b` suffix (e.g. `trace!` → `trace_b`). The active Julia runtime is
accessible as `tp.jl` for advanced use.

See the rest of this documentation for a full description of the underlying
Julia functions.
