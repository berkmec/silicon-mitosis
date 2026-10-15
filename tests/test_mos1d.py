import numpy as np

from qcell import mos1d as M


def test_triangular_well_matches_airy():
    # linear potential V = q F z with an infinite wall at z=0: E_0 = 2.338 (hbar^2 q^2 F^2 / 2 m)^(1/3)
    mos = M.Mos(zmax_nm=40.0, dz_nm=0.05)
    F = 3e7  # V/m = 0.03 V/nm
    psi = -F * mos.z * 1e-9  # Ec = Eg/2 + phif - psi = const + F z
    bands = mos._subbands(psi + 0.0)
    for w, _, v in bands:
        e_airy = 2.33811 * (M.HB_EVNM2 * (0.03) ** 2 / v["mz"]) ** (1 / 3)
        e_num = w[0] - (M.EG / 2 + mos.phif)
        assert abs(e_num - e_airy) / e_airy < 5e-3


def test_depletion_approximation_at_small_bias():
    mos = M.Mos(zmax_nm=140.0, dz_nm=0.2)
    s = mos.solve(0.45)
    W = np.sqrt(2 * M.EPS_SI * M.EPS0 * 0.45 / (M.Q * mos.na))
    psi_dep = 0.45 * np.clip(1 - mos.z * 1e-9 / W, 0, None) ** 2
    mid = (mos.z > 5) & (mos.z < 0.8 * W * 1e9)
    assert np.max(np.abs(s["psi"] - psi_dep)[mid]) < 0.05  # exponential hole tail ~ Vt
    assert s["ns_cm2"] < 1e10


def test_inversion_charge_increases_with_gate_voltage():
    mos = M.Mos(zmax_nm=100.0, dz_nm=0.2)
    a, b = mos.solve(0.9), mos.solve(1.0)
    assert b["Vg"] > a["Vg"] and b["ns_cm2"] > a["ns_cm2"]
