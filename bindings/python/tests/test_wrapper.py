import numpy as np
import pytest


def test_import():
    import testparticle as tp

    assert hasattr(tp, "trace")
    assert hasattr(tp, "trace_gc")
    assert hasattr(tp, "Proton")
    assert hasattr(tp, "Electron")
    # The Boris family lives in the companion package BorisPushers.jl.
    assert hasattr(tp, "Boris")
    assert hasattr(tp, "MultistepBoris6")
    assert hasattr(tp, "AdaptiveHybrid")


def test_basic_types():
    import testparticle as tp

    # `Species` stores mass and charge in the `m` and `q` fields.
    m_p = tp.Proton.m
    q_p = tp.Proton.q
    assert m_p > 0
    assert q_p > 0


def test_gyrofrequency():
    import testparticle as tp

    # get_gyrofrequency(B; q, m): omega = q * B / m.
    # With B = q = m = 1.0 the result is exactly 1.0.
    omega = tp.get_gyrofrequency(1.0, q=1.0, m=1.0)
    assert abs(float(omega) - 1.0) < 1e-12


def test_trace_gyroradius_circle():
    import testparticle as tp

    # Trace a proton in a uniform magnetic field along z and check that it
    # follows a circular orbit of the expected gyroradius.
    B = 1.0e-8  # [T]
    v0 = 1.0e5  # [m/s], along x
    x0 = [0.0, 0.0, 0.0]
    v = [v0, 0.0, 0.0]

    def get_B(x):
        return [0.0, 0.0, B]

    # Gyroradius r = v_perp / omega, with omega = |q| * B / m.
    q, m = tp.Proton.q, tp.Proton.m
    omega = abs(q) * B / m
    r_expected = v0 / omega

    # Trace over ~one full gyroperiod and a bit more so the orbit closes.
    t_end = 2.0 * np.pi / omega * 1.05
    dt = t_end / 2000
    res = tp.trace(x0, v, (0.0, t_end), dt, get_B, species="proton")

    # Over a full orbit the mean position is the circle center; the distance
    # from that center equals the gyroradius everywhere.
    center = res.x[:, 0:2].mean(axis=0)
    radius = np.sqrt((res.x[:, 0] - center[0]) ** 2 + (res.x[:, 1] - center[1]) ** 2)
    assert np.all(np.abs(radius - r_expected) / r_expected < 0.05)


def test_trace_save_fields():
    import testparticle as tp

    def get_B(x):
        return [0.0, 0.0, 1.0e-8]

    res = tp.trace(
        [0.0, 0.0, 0.0],
        [1.0e5, 0.0, 0.0],
        (0.0, 1.0e-7),
        5.0e-10,
        get_B,
        species="proton",
        save_fields=True,
    )
    # With save_fields the state columns are x(3), v(3), E(3), B(3).
    assert res.state.shape[1] == 12
    assert res.B.shape == (res.t.size, 3)
    assert np.allclose(res.B, [0.0, 0.0, 1.0e-8], atol=1e-12)


def test_trace_gc():
    import testparticle as tp

    def get_B(x):
        return [0.0, 0.0, 1.0e-8]

    res = tp.trace_gc(
        [0.0, 0.0, 0.0],
        [1.0e5, 0.0, 0.0],
        (0.0, 1.0e-7),
        5.0e-10,
        get_B,
        species="proton",
    )
    # Guiding-center state is (R_x, R_y, R_z, v_par).
    assert res.state.shape[1] == 4
    assert res.R.shape == (res.t.size, 3)


def test_trace_ensemble_init_func():
    import testparticle as tp

    # Ensemble tracing with a per-trajectory initial-state function. This
    # exercises the Julia-side init wrapper, which used to call the undefined
    # `int(i)` and crash the whole ensemble simulation.
    def get_B(x):
        return [0.0, 0.0, 1.0e-8]

    def init_func(i):
        return [float(i), 0.0, 0.0, 1.0e5, 0.0, 0.0]

    results = tp.trace(
        [0.0, 0.0, 0.0],
        [1.0e5, 0.0, 0.0],
        (0.0, 1.0e-7),
        5.0e-10,
        get_B,
        species="proton",
        trajectories=2,
        init_func=init_func,
    )
    assert isinstance(results, list)
    assert len(results) == 2
    # init_func(i) places trajectory i at x = i initially.
    assert np.isclose(results[0].x[0, 0], 1.0, atol=1e-9)
    assert np.isclose(results[1].x[0, 0], 2.0, atol=1e-9)


def test_trace_alg_as_julia_object():
    import testparticle as tp

    # `alg` accepts a BorisPushers algorithm object and not only its name, and
    # both must select the same solver. A Julia value used to fall through to
    # the string branch and raise TypeError.
    def get_B(x):
        return [0.0, 0.0, 1.0e-8]

    args = ([0.0, 0.0, 0.0], [1.0e5, 0.0, 0.0], (0.0, 1.0e-7), 5.0e-10, get_B)
    by_name = tp.trace(*args, species="proton", alg="MultistepBoris6")
    by_object = tp.trace(*args, species="proton", alg=tp.MultistepBoris6(n=4))
    assert by_name.t.size == by_object.t.size
    assert np.allclose(by_name.x, by_object.x)


def test_trace_julia_field_object():
    import testparticle as tp

    # A field given as a Julia callable is passed through instead of being
    # treated as a Python callable.
    jl_B = tp.jl.seval("x -> TestParticle.SVector{3, Float64}(0.0, 0.0, 1.0e-8)")
    res = tp.trace(
        [0.0, 0.0, 0.0],
        [1.0e5, 0.0, 0.0],
        (0.0, 1.0e-7),
        5.0e-10,
        jl_B,
        species="proton",
    )
    assert res.x.shape[0] > 1
    assert np.allclose(res.v[0], [1.0e5, 0.0, 0.0], rtol=1.0e-6, atol=1.0)


def test_field_time_dependence_detected():
    import testparticle as tp

    # A field declared as ``def B(x, t=0.0)`` can accept a time argument and
    # must be treated as time-dependent: the solver must pass the time
    # parameter through instead of always calling it with the default t = 0.
    times = []

    def get_B(x, t=0.0):
        times.append(t)
        return [0.0, 0.0, 1.0e-8]

    tp.trace(
        [0.0, 0.0, 0.0],
        [1.0e5, 0.0, 0.0],
        (0.0, 1.0e-7),
        5.0e-10,
        get_B,
        species="proton",
    )
    # If the field were misclassified as time-independent the solver would only
    # ever invoke it with t = 0.0; non-zero recorded times prove the time
    # parameter is actually honoured.
    assert len(times) > 0
    assert any(t != 0.0 for t in times)
