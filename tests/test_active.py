import numpy as np

from qcell.active import ActiveCH, rates_charging, rates_zwicker
from qcell.grayscott import GrayScott, rhs_kinetics, uniform_states


def test_leaky_charging_is_single_line_and_zwicker_rates():
    m = ActiveCH(n=16, **rates_charging(J=2e-3, lam=1e-2))
    c = np.linspace(-0.2, 1.2, 57)
    assert np.allclose(m.s(c), 2e-3 - 1e-2 * c, atol=1e-15)  # dnu = 0: one straight line through both phases
    r = rates_zwicker(eps=0.19, dnu=-1.1e-3)
    assert np.isclose(r["nu_p"], 1.9e-3) and np.isclose(r["nu_m"], 7e-3)  # Zwicker Fig. 1C values


def test_gray_scott_uniform_states_and_spot_count():
    for F, k in ((0.02, 0.05), (0.04, 0.06), (0.1, 0.03)):
        for u, v in uniform_states(F, k):
            du, dv = rhs_kinetics(u, v, F, k)
            assert abs(du) < 1e-14 and abs(dv) < 1e-14
    assert len(uniform_states(0.02, 0.059)) == 1  # Pearson: spots live where (1, 0) is the only uniform state
    gs = GrayScott(n=64)
    V = np.zeros((64, 64))
    V[5:9, 5:9] = 0.3
    V[30:34, 0:3] = 0.3
    V[30:34, 62:64] = 0.3  # wraps across the periodic edge: still one spot
    V[50, :] = 0.3  # a stripe that wraps the box
    out = gs.spots(V)
    assert out["n"] == 2 and out["perc"]


def test_droplet_count_periodic_and_percolation():
    m = ActiveCH(n=64, L=32.0)
    c = m.disk(4.0, center=(0.5, 16.0))  # straddles the periodic edge
    c = np.maximum(c, m.disk(3.0, center=(16.0, 16.0)))
    out = m.droplets(c)
    assert out["n"] == 2 and not out["perc"]
    assert abs(sorted(out["radii"])[1] - 4.0) < 0.3
    c[:, 40:44] = 1.0
    assert m.droplets(c)["perc"]


def test_reaction_flux_continuous_with_paper_branches():
    m = ActiveCH(n=16, nu_p=2e-3, nu_m=7e-3, k_p=1e-2, k_m=1e-2)
    eps = 1e-7
    for x in (m.cc_p, m.cc_m):
        assert abs(m.s(np.array([x - eps]))[0] - m.s(np.array([x + eps]))[0]) < 1e-8
        d_lo = (m.s(np.array([x - eps]))[0] - m.s(np.array([x - 2 * eps]))[0]) / eps
        d_hi = (m.s(np.array([x + 2 * eps]))[0] - m.s(np.array([x + eps]))[0]) / eps
        assert abs(d_lo - d_hi) < 1e-5
    assert np.isclose(m.s(np.array([0.1]))[0], 2e-3 - 1e-2 * 0.1)
    assert np.isclose(m.s(np.array([0.9]))[0], -7e-3 - 1e-2 * (0.9 - 1.0))


def test_flat_interface_profile_and_surface_tension():
    # passive (s = 0): stationary profile c = 1/2 + 1/2 tanh(x / w), w = 1; gamma = sqrt(kappa b)/6 = 1/12 (b=1, kappa=1/4)
    m = ActiveCH(n=256, L=64.0)
    x = (np.arange(256) + 0.5) * m.dx
    prof = 0.5 * (1 + np.tanh((x - 16) / 1.7)) - 0.5 * (1 + np.tanh((x - 48) / 1.7))  # wrong width on purpose
    c = np.tile(prof[:, None], (1, 256))[:, :256]
    m2 = ActiveCH(n=256, L=64.0)  # passive: no chemistry, no Coulomb
    assert np.all(m2.s(np.linspace(0, 1, 11)) == 0)
    for _ in range(4000):
        c = m2.step(c, 0.1)
    p = c[:, 0]
    i0 = np.argmin(np.abs(p[:128] - 0.5))
    xs = x[i0 - 6:i0 + 7]
    exact = 0.5 * (1 + np.tanh(xs - np.interp(0.5, p[i0 - 6:i0 + 7], xs)))
    assert np.max(np.abs(p[i0 - 6:i0 + 7] - exact)) < 0.02
    dcdx = np.gradient(p, m.dx)
    energy = np.sum(0.5 * p**2 * (p - 1) ** 2 + 0.125 * dcdx**2) * m.dx
    assert abs(energy / 2 - 1 / 12) / (1 / 12) < 0.02


def test_far_field_relaxes_to_c_infinity():
    m = ActiveCH(n=32, nu_p=2e-3, nu_m=7e-3, k_p=1e-2, k_m=1e-2)
    c = np.zeros((32, 32))
    for _ in range(20000):
        c = m.step(c, 0.1)
    assert abs(c.mean() - 2e-3 / 1e-2) < 1e-3  # c_inf = nu+ / k+


def test_mass_conserved_without_reactions():
    m = ActiveCH(n=64, A=0.02, d=3.0)
    rng = np.random.default_rng(0)
    c = 0.4 + 0.05 * rng.standard_normal((64, 64))
    c1 = c.copy()
    for _ in range(200):
        c1 = m.step(c1, 0.1)
    assert abs(c1.mean() - c.mean()) < 1e-12
