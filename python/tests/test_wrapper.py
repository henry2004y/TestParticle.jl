import pytest

def test_import():
    import testparticle as tp
    assert hasattr(tp, "trace")
    assert hasattr(tp, "trace_b")
    assert hasattr(tp, "Proton")
    assert hasattr(tp, "Electron")

def test_basic_types():
    import testparticle as tp
    from juliacall import Main as jl
    
    # Check constants
    m_p = tp.Proton.mass
    q_p = tp.Proton.charge
    assert m_p > 0
    assert q_p > 0

def test_gyrofrequency():
    import testparticle as tp
    
    # B = 1.0 (scalar or magnitude), q = 1.0, m = 1.0
    # The signature in Julia is usually get_gyrofrequency(B, species) or (B, q, m)
    # Let's check (B, q, m)
    
    omega = tp.get_gyrofrequency(1.0, 1.0, 1.0)
    assert abs(omega - 1.0) < 1e-6

