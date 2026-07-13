"""Python wrapper for TestParticle.jl.

This wrapper exposes the most useful objects of the Julia package
``TestParticle`` and provides a small set of high-level helpers
(``trace``, ``trace_gc``, ``trace_fieldline``) that accept and return
``numpy`` arrays so that the Boris / guiding-center / field-line solvers
can be driven from Python without touching Julia types directly.

The underlying Julia runtime is managed by ``juliacall``. Functions ending
with ``!`` in Julia are exposed with the ``_b`` suffix (e.g. ``trace!`` is
available as ``trace_b``).
"""

from __future__ import annotations

import inspect
import warnings
from dataclasses import dataclass

import numpy as np
from juliacall import Main as jl

# Load the Julia package. The Julia environment (including TestParticle itself)
# is resolved from ``juliapkg.json`` next to this package.
try:
    jl.seval("using TestParticle")
    jl.seval("using PythonCall")  # exposes `pyconvert` for wrapping Python callables
except Exception as exc:  # pragma: no cover - depends on the user environment
    raise RuntimeError(
        "Failed to `using TestParticle` from Julia. Make sure the Julia "
        "dependencies declared in juliapkg.json are instantiated (this happens "
        "automatically on first `import testparticle`)."
    ) from exc

TP = jl.TestParticle


# Helper that converts a list of TestParticle solutions into Python-friendly
# ``(times, states)`` where ``times`` is a list of 1D arrays and ``states`` a
# list of 2D arrays with shape ``(n_steps, n_vars)``.
jl.seval(
    """
    function _tp_extract(sols)
        sol_list = sols isa TestParticle.EnsembleSolution ? sols.u : sols
        n = length(sol_list)
        ts = Vector{Vector{Float64}}(undef, n)
        us = Vector{Matrix{Float64}}(undef, n)
        for i in 1:n
            s = sol_list[i]
            t = collect(Float64, s.t)
            u = reduce(hcat, [Float64.(uu) for uu in s.u])
            us[i] = permutedims(u)
            ts[i] = t
        end
        return (ts, us)
    end

    function _tp_extract_fieldline(sol)
        _tp_extract((sol,))[1][1], _tp_extract((sol,))[2][1]
    end
    """
)

# ---------------------------------------------------------------------------
# Curated re-exports of the raw Julia API for advanced usage.
#
# Symbols are bound defensively: if a name is absent from the loaded
# ``TestParticle`` module (e.g. due to a version skew or a stale Julia
# precompile cache), we emit a warning and skip it instead of crashing the
# whole import. Functions ending with ``!`` in Julia are exposed with the
# ``_b`` suffix (e.g. ``trace!`` is available as ``trace_b``).
# ---------------------------------------------------------------------------
_REEXPORTS = [
    ("prepare", "prepare"),
    ("prepare_gc", "prepare_gc"),
    ("get_gc", "get_gc"),
    ("get_gc_func", "get_gc_func"),
    ("solve", "solve"),
    ("TraceProblem", "TraceProblem"),
    ("TraceGCProblem", "TraceGCProblem"),
    ("TraceHybridProblem", "TraceHybridProblem"),
    ("TraceFieldlineProblem", "TraceFieldlineProblem"),
    ("Proton", "Proton"),
    ("Electron", "Electron"),
    ("Ion", "Ion"),
    ("ZeroField", "ZeroField"),
    ("Maxwellian", "Maxwellian"),
    ("BiMaxwellian", "BiMaxwellian"),
    ("Kappa", "Kappa"),
    ("BiKappa", "BiKappa"),
    ("AdaptiveBoris", "AdaptiveBoris"),
    ("AdaptiveMultistepBoris", "AdaptiveMultistepBoris"),
    ("AdaptiveHybrid", "AdaptiveHybrid"),
    ("Boris", "Boris"),
    ("MultistepBoris", "MultistepBoris"),
    ("MultistepBoris2", "MultistepBoris2"),
    ("MultistepBoris4", "MultistepBoris4"),
    ("MultistepBoris6", "MultistepBoris6"),
    ("get_gyrofrequency", "get_gyrofrequency"),
    ("get_gyroperiod", "get_gyroperiod"),
    ("get_gyroradius", "get_gyroradius"),
    ("get_velocity", "get_velocity"),
    ("get_energy", "get_energy"),
    ("get_mean_magnitude", "get_mean_magnitude"),
    ("energy2velocity", "energy2velocity"),
    ("get_curvature_radius", "get_curvature_radius"),
    ("get_adiabaticity", "get_adiabaticity"),
    ("adiabaticity_components", "adiabaticity_components"),
    ("sample_unit_sphere", "sample_unit_sphere"),
    ("generate_sphere", "generate_sphere"),
    ("sample_maxwellian", "sample_maxwellian"),
    ("orbit", "orbit"),
    ("monitor", "monitor"),
    ("get_fields", "get_fields"),
    ("get_work", "get_work"),
    ("LazyTimeInterpolator", "LazyTimeInterpolator"),
    ("CartesianGrid", "CartesianGrid"),
    ("RectilinearGrid", "RectilinearGrid"),
    ("StructuredGrid", "StructuredGrid"),
    ("trace_b", "trace!"),
    ("trace_relativistic_b", "trace_relativistic!"),
    ("trace_normalized_b", "trace_normalized!"),
    ("trace_relativistic_normalized_b", "trace_relativistic_normalized!"),
    ("trace_relativistic", "trace_relativistic"),
    ("trace_normalized", "trace_normalized"),
    ("trace_relativistic_normalized", "trace_relativistic_normalized"),
    ("trace_gc_b", "trace_gc!"),
    ("trace_gc_drifts_b", "trace_gc_drifts!"),
    ("trace_gc_flr_b", "trace_gc_flr!"),
    ("trace_gc_exb_b", "trace_gc_exb!"),
    ("trace_fieldline_b", "trace_fieldline!"),
]

_BOUND = []
for _py_name, _jl_name in _REEXPORTS:
    try:
        globals()[_py_name] = getattr(TP, _jl_name)
        _BOUND.append(_py_name)
    except AttributeError:
        warnings.warn(
            f"TestParticle.{_jl_name} is not available in the loaded Julia "
            f"package; '{_py_name}' was not re-exported.",
            RuntimeWarning,
            stacklevel=2,
        )

_SPECIES = {"proton": Proton, "electron": Electron}

_field_counter = 0


def _is_julia(obj):
    return getattr(obj, "__julia__", None) is not None


def _as_species(species):
    if isinstance(species, str):
        return _SPECIES[species.lower()]
    return species


def _wrap_field(func):
    """Wrap a Python callable into a Julia ``Function`` usable by TestParticle.

    The callable may have the signature ``f(x)`` (time-independent) or
    ``f(x, t)`` (time-dependent) and must return a sequence of three numbers.
    """
    global _field_counter
    _field_counter += 1
    base = f"_tp_field_{_field_counter}"
    # Store the raw Python callable; the Julia closure converts its return
    # value (a list or numpy array) into a native SVector{3,Float64} via
    # pyconvert. Returning a Python object directly to Julia would leave a
    # wrapped `Py` value that the solvers cannot operate on.
    setattr(jl, base, func)
    try:
        sig = inspect.signature(func)
        n_pos = sum(
            1
            for p in sig.parameters.values()
            if p.kind in (p.POSITIONAL_ONLY, p.POSITIONAL_OR_KEYWORD)
            and p.default is p.empty
        )
    except (ValueError, TypeError):
        n_pos = 1
    if n_pos >= 2:
        jl.seval(
            f"{base}_jl = (x, t) -> TestParticle.SVector{{3, Float64}}("
            f"pyconvert(Vector{{Float64}}, {base}(Vector(x), float(t))))"
        )
    else:
        jl.seval(
            f"{base}_jl = (x) -> TestParticle.SVector{{3, Float64}}("
            f"pyconvert(Vector{{Float64}}, {base}(Vector(x))))"
        )
    return getattr(jl, f"{base}_jl")


def _wrap_field_or_none(field):
    if field is None:
        return None
    if _is_julia(field):
        return field
    if callable(field):
        return _wrap_field(field)
    raise TypeError("field must be a callable or a Julia object")


def _wrap_init_func(func):
    """Wrap a Python ``init_func(i) -> [x, y, z, vx, vy, vz]`` for ensembles."""
    global _field_counter
    _field_counter += 1
    base = f"_tp_init_{_field_counter}"
    setattr(jl, base, func)
    jl.seval(
        f"{base}_jl = (prob, i, repeat) -> TestParticle.SVector{{6, Float64}}("
        f"pyconvert(Vector{{Float64}}, {base}(int(i))))"
    )
    return getattr(jl, f"{base}_jl")


def _to_svector(arr):
    return TP.SVector(*[float(v) for v in np.asarray(arr, dtype=float).ravel()])


# ---------------------------------------------------------------------------
# Result containers.
# ---------------------------------------------------------------------------
@dataclass
class TraceResult:
    """Trajectory of a full particle traced with the Boris method."""

    t: np.ndarray
    state: np.ndarray
    save_fields: bool = False
    save_work: bool = False

    @property
    def x(self) -> np.ndarray:
        return self.state[:, 0:3]

    @property
    def v(self) -> np.ndarray:
        return self.state[:, 3:6]

    @property
    def E(self) -> np.ndarray:
        if not self.save_fields:
            raise AttributeError("E field was not saved; pass save_fields=True")
        return self.state[:, 6:9]

    @property
    def B(self) -> np.ndarray:
        if not self.save_fields:
            raise AttributeError("B field was not saved; pass save_fields=True")
        return self.state[:, 9:12]

    @property
    def work(self) -> np.ndarray:
        if not self.save_work:
            raise AttributeError("work was not saved; pass save_work=True")
        return self.state[:, 12:16]


@dataclass
class GCTraceResult:
    """Trajectory of a guiding center traced with the RK4/RK45 method."""

    t: np.ndarray
    state: np.ndarray
    save_fields: bool = False
    save_work: bool = False

    @property
    def R(self) -> np.ndarray:
        return self.state[:, 0:3]

    @property
    def vpar(self) -> np.ndarray:
        return self.state[:, 3]

    @property
    def E(self) -> np.ndarray:
        if not self.save_fields:
            raise AttributeError("E field was not saved; pass save_fields=True")
        return self.state[:, 4:7]

    @property
    def B(self) -> np.ndarray:
        if not self.save_fields:
            raise AttributeError("B field was not saved; pass save_fields=True")
        return self.state[:, 7:10]

    @property
    def work(self) -> np.ndarray:
        if not self.save_work:
            raise AttributeError("work was not saved; pass save_work=True")
        return self.state[:, 10:14]


# ---------------------------------------------------------------------------
# High-level API.
# ---------------------------------------------------------------------------
def trace(
    x0,
    v0,
    tspan,
    dt,
    B,
    E=None,
    *,
    species="proton",
    trajectories=1,
    savestepinterval=1,
    save_fields=False,
    save_work=False,
    parallel=False,
    alg=None,
    init_func=None,
    q=None,
    m=None,
):
    """Trace charged particle trajectories with the Boris pusher.

    Parameters
    ----------
    x0, v0 : array_like
        Initial position and velocity, each of length 3 (in metres and
        metres per second).
    tspan : tuple
        ``(t0, t1)`` integration span in seconds.
    dt : float
        Fixed time step in seconds.
    B, E : callable or None
        Magnetic / electric field as a Python callable ``f(x)`` or
        ``f(x, t)`` returning a length-3 sequence, or a Julia field object.
        ``E`` may be ``None`` (no electric field).
    species : str or Julia species, optional
        ``"proton"`` (default), ``"electron"``, a ``tp.Ion(m, q)``, or any
        Julia ``Species`` object.
    trajectories : int, optional
        Number of trajectories. Use ``init_func`` to vary the initial state
        of each trajectory.
    savestepinterval, save_fields, save_work : optional
        Mirrors of the Julia ``solve`` keywords.
    parallel : bool, optional
        Use ``EnsembleThreads`` when ``True``.
    alg : str or Julia object, optional
        Boris solver variant, e.g. ``"Boris"``, ``"MultistepBoris"`` passed
        to the Julia ``solve`` (defaults to the standard Boris method).
    init_func : callable, optional
        ``init_func(i) -> [x, y, z, vx, vy, vz]`` returning the initial state
        of trajectory ``i`` (1-based).
    q, m : float, optional
        Charge and mass overriding those of ``species``.

    Returns
    -------
    TraceResult or list[TraceResult]
        A single result when ``trajectories == 1``, otherwise a list.
    """
    x0 = np.asarray(x0, dtype=float)
    v0 = np.asarray(v0, dtype=float)
    if x0.shape != (3,) or v0.shape != (3,):
        raise ValueError("x0 and v0 must both have length 3")

    species_jl = _as_species(species)
    B_jl = _wrap_field_or_none(B)
    E_jl = _wrap_field_or_none(E)

    prepare_kw = {}
    if q is not None:
        prepare_kw["q"] = float(q)
    if m is not None:
        prepare_kw["m"] = float(m)

    if E_jl is None:
        param = TP.prepare(B_jl, species=species_jl, **prepare_kw)
    else:
        param = TP.prepare(E_jl, B_jl, species=species_jl, **prepare_kw)

    u0 = _to_svector(np.concatenate([x0, v0]))
    if init_func is not None:
        prob_func = _wrap_init_func(init_func)
        prob = TP.TraceProblem(u0, tuple(tspan), param, prob_func=prob_func)
    else:
        prob = TP.TraceProblem(u0, tuple(tspan), param)

    solve_kw = dict(
        dt=float(dt),
        trajectories=int(trajectories),
        savestepinterval=int(savestepinterval),
        save_fields=bool(save_fields),
        save_work=bool(save_work),
    )
    if alg is None:
        alg_jl = TP.Boris()
    elif _is_julia(alg):
        alg_jl = alg
    else:
        alg_jl = getattr(TP, alg)()
    if parallel:
        sols = TP.solve(prob, alg_jl, TP.EnsembleThreads(), **solve_kw)
    else:
        sols = TP.solve(prob, alg_jl, **solve_kw)

    ts, us = jl._tp_extract(sols)
    ts_py = [np.asarray(t) for t in ts]
    us_py = [np.asarray(u) for u in us]

    if trajectories == 1 and init_func is None:
        return TraceResult(ts_py[0], us_py[0], save_fields, save_work)
    return [TraceResult(t, u, save_fields, save_work) for t, u in zip(ts_py, us_py)]


def trace_gc(
    x0,
    v0,
    tspan,
    dt,
    B,
    E=None,
    *,
    species="proton",
    trajectories=1,
    savestepinterval=1,
    save_fields=False,
    save_work=False,
    alg="rk4",
    q=None,
    m=None,
):
    """Trace guiding-center trajectories with the RK4/RK45 method.

    Parameters mirror :func:`trace`, with the addition of ``alg``
    (``"rk4"`` or ``"rk45"``). The returned :class:`GCTraceResult` exposes the
    guiding-center position ``R`` and parallel velocity ``vpar``.
    """
    x0 = np.asarray(x0, dtype=float)
    v0 = np.asarray(v0, dtype=float)
    if x0.shape != (3,) or v0.shape != (3,):
        raise ValueError("x0 and v0 must both have length 3")

    species_jl = _as_species(species)
    B_jl = _wrap_field_or_none(B)
    E_jl = _wrap_field_or_none(E)
    E_arg = E_jl if E_jl is not None else TP.ZeroField()

    prepare_kw = {}
    if q is not None:
        prepare_kw["q"] = float(q)
    if m is not None:
        prepare_kw["m"] = float(m)

    stateinit_gc, param = TP.prepare_gc(
        _to_svector(np.concatenate([x0, v0])), E_arg, B_jl,
        species=species_jl, **prepare_kw
    )
    u0 = _to_svector(stateinit_gc)
    prob = TP.TraceGCProblem(u0, tuple(tspan), param)

    solve_kw = dict(
        dt=float(dt),
        trajectories=int(trajectories),
        savestepinterval=int(savestepinterval),
        save_fields=bool(save_fields),
        save_work=bool(save_work),
        alg=jl.Symbol(alg),
    )
    sols = TP.solve(prob, **solve_kw)

    ts, us = jl._tp_extract(sols)
    ts_py = [np.asarray(t) for t in ts]
    us_py = [np.asarray(u) for u in us]

    if trajectories == 1:
        return GCTraceResult(ts_py[0], us_py[0], save_fields, save_work)
    return [GCTraceResult(t, u, save_fields, save_work) for t, u in zip(ts_py, us_py)]


def trace_fieldline(x0, B, tspan, *, mode="both", alg="Tsit5", abstol=1e-8, reltol=1e-8):
    """Trace magnetic field lines.

    Parameters
    ----------
    x0 : array_like
        Initial position of length 3.
    B : callable or Julia object
        Magnetic field as ``f(x)`` returning a length-3 sequence, or a Julia
        field object.
    tspan : tuple
        Arc-length span ``(s0, s1)``.
    mode : str, optional
        ``"forward"``, ``"backward"`` or ``"both"`` (default).
    alg : str or Julia object, optional
        OrdinaryDiffEq integrator name (e.g. ``"Tsit5"``) or a Julia alg.

    Returns
    -------
    dict
        Dictionary with ``"forward"`` and/or ``"backward"`` entries, each a
        ``(t, x)`` tuple of numpy arrays (``x`` has shape ``(n_steps, 3)``).
    """
    x0 = np.asarray(x0, dtype=float)
    if x0.shape != (3,):
        raise ValueError("x0 must have length 3")

    B_jl = _wrap_field_or_none(B)
    u0 = _to_svector(x0)
    probs = TP.TraceFieldlineProblem(u0, B_jl, tuple(tspan), mode=jl.Symbol(mode))

    jl.seval("using OrdinaryDiffEq")
    if _is_julia(alg):
        alg_jl = alg
    else:
        alg_jl = getattr(jl, alg)()

    def _solve_one(prob):
        sol = jl.solve(prob, alg_jl, abstol=float(abstol), reltol=float(reltol))
        t, u = jl._tp_extract_fieldline(sol)
        return np.asarray(t), np.asarray(u)

    out = {}
    if mode == "both":
        probs_py = list(probs)
        out["forward"] = _solve_one(probs_py[0])
        out["backward"] = _solve_one(probs_py[1])
    else:
        out[mode] = _solve_one(probs)
    return out


__version__ = "0.1.0"

__all__ = [
    "jl",
    *_BOUND,
    "trace", "trace_gc", "trace_fieldline",
    "TraceResult", "GCTraceResult",
    "__version__",
]
