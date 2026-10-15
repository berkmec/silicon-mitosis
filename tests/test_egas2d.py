import numpy as np

from qcell import egas2d as G
from qcell.model import E, EPS0, HBAR, M0


def test_exchange_unscreened_matches_analytic():
    # per particle of one species: -(4/(3 pi)) e^2 kF / (4 pi eps0 eps)
    eps, n, g = 7.8, 1e16, 2
    ex = G.exchange_per_area(n, g, G.V_uniform(eps))
    kF = np.sqrt(4 * np.pi * n / g)
    ref = -(4 / (3 * np.pi)) * E**2 * kF / (4 * np.pi * EPS0 * eps) * n
    assert abs(ex / ref - 1) < 1e-8


def test_exchange_in_hartree_units_matches_paper_form():
    # Attaccalite et al.: eps_x(rs, 0) = -4 sqrt(2) / (3 pi rs) Hartree
    m, eps, n = 0.19, 7.8, 3e15
    Ha, _ = G.units(m, eps)
    rs = G.rs_of(n, m, eps)
    ex = G.exchange_per_area(n, 2, G.V_uniform(eps)) / n / Ha
    assert abs(ex / (-4 * np.sqrt(2) / (3 * np.pi * rs)) - 1) < 1e-8


def test_hf_critical_rs_single_valley():
    # HF, g=2: K = 0 at rs = pi/sqrt(2)
    m, eps = 0.19, 7.8
    _, aB = G.units(m, eps)
    rs_c = np.pi / np.sqrt(2)
    n_c = 1 / (np.pi * aB**2 * rs_c**2)
    kw = dict(g=2, m_eff=m, V_x=G.V_uniform(eps), eps_c=eps, with_corr=False)
    assert G.inverse_compressibility(n_c * 1.05, **kw) > 0 > G.inverse_compressibility(n_c * 0.95, **kw)


def test_gate_kernel_limits():
    V = G.V_gate(11.7, 3.9, 4.148e-9)
    q_big = 5e10
    assert abs(V(q_big) / G.V_uniform((11.7 + 3.9) / 2)(q_big) - 1) < 1e-6  # t_ox -> infinity limit
    q_small = 1e5
    cap = E**2 * 4.148e-9 / (EPS0 * 3.9)  # e^2 / C_ox per area
    first_order = 1 / (1 + 11.7 * q_small * 4.148e-9 / 3.9)  # V -> cap / (1 + eps_si q t / eps_ox)
    assert abs(V(q_small) / (cap * first_order) - 1) < 1e-6


def test_attaccalite_total_energy_against_paper_table():
    # Table I (FN-DMC, zeta=0): rs=2: -0.255802, 5: -0.149134, 10: -0.0852706, 20: -0.046241 (Hartree)
    for rs, ref in ((2, -0.255802), (5, -0.149134), (10, -0.0852706), (20, -0.046241)):
        t = 1 / (2 * rs**2)
        ex = -4 * np.sqrt(2) / (3 * np.pi * rs)
        tot = t + ex + G.alpha0(rs)
        assert abs(tot / ref - 1) < 0.015, (rs, tot, ref)


def test_kinetic_unit_consistency():
    m, n, g = 0.19, 1e16, 2
    ek = G.kinetic_per_area(n, g, m) / n
    kF = np.sqrt(4 * np.pi * n / g)
    assert abs(ek / (HBAR**2 * kF**2 / (2 * m * M0) / 2) - 1) < 1e-12
