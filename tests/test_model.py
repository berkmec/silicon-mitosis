import numpy as np

from qcell import model as m

P = m.CellParams()


def test_charging_energy_scales_inverse_size():
    e0, e1 = m.charging_energy(P, 0), m.charging_energy(P, 1)
    assert np.isclose(e1 / e0, 1 / P.r)


def test_level_spacing_scales_inverse_square():
    s0, s1 = m.level_spacing(P, 0), m.level_spacing(P, 1)
    assert np.isclose(s1 / s0, 1 / P.r**2)


def test_regime_monotonic_in_temperature():
    g300 = m.regime_generation(P, 300)[0]
    g4 = m.regime_generation(P, 4)[0]
    assert g4 is not None and (g300 is None or g4 <= g300)


def test_current_zero_at_zero_bias():
    I = m.single_dot_current(1e-18, 1e-18, 1e-18, np.linspace(0, 0.3, 20), [0.0], T=4)
    assert np.allclose(I, 0, atol=1e-18)


def test_coulomb_blockade_and_peaks():
    Cs = Cd = Cg = 1e-18
    Ec = m.E**2 / (2 * (Cs + Cd + Cg))
    Vg = np.linspace(0, 0.6, 400)
    low = np.abs(m.single_dot_current(Cs, Cd, Cg, Vg, [0.1 * Ec / m.E], T=1.0))[:, 0]
    assert low.max() > 50 * np.median(low)  # sharp peaks on a blockaded floor
    n_peaks = ((low[1:-1] > low[:-2]) & (low[1:-1] > low[2:]) & (low[1:-1] > 0.5 * low.max())).sum()
    assert n_peaks >= 2


def test_spawn_voltage_rises_with_depth():
    v = [m.spawn_voltage(P, g) for g in range(5)]
    assert all(b > a for a, b in zip(v, v[1:]))


def test_cascade_stops_by_breakdown_for_large_pulse():
    g, why = m.cascade_depth(P, v_pulse=50.0)
    assert why == "breakdown"
