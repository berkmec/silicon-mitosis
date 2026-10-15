"""Exp14 checks: d-dimensional sharp-interface theory, axisymmetric (z, r) solver and lean Cartesian d-dim solver."""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "experiments" / "exp14_dimension"))


@pytest.mark.parametrize("d", [2.0, 3.0, 4.5, 6.0])
def test_linear_theory_mode_identities(d):
    from qcell import active_sharp as a

    eps, k = 0.2, 0.01
    nu_m = -0.001 + k - eps * k
    R = a.stable_radius(d, eps, nu_m, k)
    h = 1e-5
    dv = (a.velocity(R + h, d, eps, nu_m, k) - a.velocity(R - h, d, eps, nu_m, k)) / (2 * h)
    assert abs(a.sigma(R, d, 1, 0.0, eps, nu_m, k)) < 1e-12              # translation is neutral
    assert a.sigma(R, d, 0, 0.0, eps, nu_m, k) == pytest.approx(dv, rel=1e-6)


@pytest.mark.parametrize("d", [3.0, 4.5, 6.0])
def test_passive_limit_gives_rayleigh_plateau_threshold(d):
    """No chemistry (k -> 0), far field at the Gibbs-Thomson value: the tube's varicose rate changes sign at kz R = sqrt(d-2)."""
    from qcell import active_sharp as a

    R, k = 5.0, 1e-10
    eps = a.LG * (d - 2.0) / R
    kz = np.linspace(0.02, 2.0, 4000) / R
    s = a.sigma(R, d - 1.0, 0, kz, eps, 0.0, k)
    kc = kz[np.nonzero(np.diff(np.sign(s)))[0][0]]
    assert kc * R == pytest.approx(np.sqrt(d - 2.0), abs=2e-3)


def test_nonlinear_theory_mode_identities():
    from qcell import active_sharp_nl as nl

    ch = nl.Chem(0.17, 0.0083, 0.01)
    R = nl.stable_radius(2.0, ch)
    b = nl.base(R, 2.0, ch)
    h = 1e-3
    dv = (nl.velocity(R + h, 2.0, ch) - nl.velocity(R - h, 2.0, ch)) / (2 * h)
    assert abs(nl.sigma(R, 2.0, 1, 0.0, ch, b=b)) < 1e-8
    assert nl.sigma(R, 2.0, 0, 0.0, ch, b=b) == pytest.approx(dv, rel=1e-3)
    assert R == pytest.approx(4.69, rel=0.03)                                # exp11 2D phase-field droplet


def test_axisymmetric_solver_conserves_mass():
    from qcell.active_axi import AxiCH

    m = AxiCH(4.5, 64, 24, h=0.5, dt=0.5)
    c = m.tube(4.0, amp=0.1)
    M0 = m.mass(c)
    c = m.step(c, 50)
    assert abs(m.mass(c) - M0) / M0 < 1e-11


def test_axisymmetric_d2_converges_to_mirror_symmetric_cartesian_2d():
    """d = 2 axisymmetric = mirror-symmetric 2D problem; the r discretisation is 2nd order, so halving h must cut the
    difference to the spectral Cartesian solution by ~4."""
    from qcell.active import ActiveCH
    from qcell.active_axi import AxiCH

    err = []
    for n, h in ((64, 0.5), (128, 0.25)):
        ax = AxiCH(2.0, n, n // 2, h=h, dt=0.1, nu_p=0.002, nu_m=0.007, k_p=0.01, k_m=0.01)
        c0 = ax.droplet(6.0, ecc=0.2)
        full = np.concatenate([c0[:, ::-1], c0], axis=1)                       # mirror image in y: periodic 2D box
        ca = ActiveCH(n=n, L=n * h, dim=2, S=1.0, nu_p=0.002, nu_m=0.007, k_p=0.01, k_m=0.01).run(full, 8.0, 0.1)
        err.append(np.abs(ca[:, n // 2:] - ax.step(c0, 80)).max())
    assert err[0] < 1e-2 and err[0] / err[1] > 3.0


def test_lean_cartesian_solver_matches_activech():
    import cart_run as C
    from qcell.active import ActiveCH, rates_zwicker

    r = rates_zwicker(0.2, -0.001)
    m = C.CartCH((16, 16, 16), 0.75, 0.2, S=1.0, **r)
    c0 = m.droplet(4.0, 0.05)
    ca = ActiveCH(n=16, L=12.0, dim=3, S=1.0, dtype=np.float32, **r).run(c0.copy(), 4.0, 0.2)
    assert np.abs(ca - m.step(c0.copy(), 20)).max() < 1e-4


def test_periodic_components_merge_across_corner():
    import cart_run as C

    m = C.CartCH((10, 10, 10, 10), 1.0, 0.2)
    x = np.roll(m.droplet(2.5, 0.0), (5, 5, 5, 5), axis=(0, 1, 2, 3))
    assert len(C.components(x, 1.0)) == 1
