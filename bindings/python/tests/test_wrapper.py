import numpy as np
import pytest


def test_import():
    import testparticle as tp

    assert hasattr(tp, "trace")
    assert hasattr(tp, "trace_gc")
    assert hasattr(tp, "trace_fieldline")
    assert hasattr(tp, "Proton")
    assert hasattr(tp, "Electron")


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

    # Over a full orbit the mean position is the circle centre; the distance
    # from that centre equals the gyroradius everywhere.
    centre = res.x[:, 0:2].mean(axis=0)
    radius = np.sqrt((res.x[:, 0] - centre[0]) ** 2 + (res.x[:, 1] - centre[1]) ** 2)
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
