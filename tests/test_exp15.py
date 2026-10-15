"""Exp15 PIC hizli testleri (< 30 s): yukleme, Poisson, oz-kuvvet, momentum, Boris, dispersiyon, yapi takibi."""
import os
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "experiments" / "exp15_phase_space"))
import pic  # noqa: E402
import structures as S  # noqa: E402


@pytest.mark.parametrize("d", [1, 2, 3])
def test_deposit_conserves_charge_and_uniform_lattice(d):
    N = (8,) * d
    L = (4.0,) * d
    g = [((np.arange(n * 2) + 0.5) / (n * 2)) * l for n, l in zip(N, L)]
    mesh = np.meshgrid(*g, indexing="ij")
    x = np.stack([m.ravel() for m in mesh])
    sim = pic.PIC(L, N, x, np.zeros((d, x.shape[1])), 0.1)
    ne, phi, E = sim.fields()
    assert np.isclose(ne.mean(), 1.0)
    assert np.allclose(ne, 1.0, atol=1e-12)
    assert all(np.abs(e).max() < 1e-12 for e in E)


def test_poisson_sinusoid_spectral():
    Ng, L, a = 64, 2 * np.pi * 2, 0.01
    k = 2 * np.pi / L
    sim = pic.PIC([L], [Ng], np.zeros((1, 4)), np.zeros((1, 4)), 0.1)
    xg = np.arange(Ng) * L / Ng
    phi, E = sim.solve(1.0 + a * np.cos(k * xg))
    assert np.allclose(E[0], -(a / k) * np.sin(k * xg), atol=1e-12)
    assert np.allclose(phi, -(a / k ** 2) * np.cos(k * xg), atol=1e-12)


@pytest.mark.parametrize("d", [1, 2, 3])
def test_zero_self_force(d):
    rng = np.random.default_rng(1)
    L, N = (5.0,) * d, (16,) * d
    x = rng.random((d, 1)) * 5.0
    sim = pic.PIC(L, N, x, np.zeros((d, 1)), 0.1)
    ne, phi, E = sim.fields()
    Ep = sim.gather(E)
    assert np.abs(Ep).max() < 1e-10


def test_momentum_conservation_2d():
    rng = np.random.default_rng(2)
    L, N = (6.0, 6.0), (16, 16)
    x, v = pic.beams(L, 4000, 1.0, 0.2, 2, rng)
    sim = pic.PIC(L, N, x, v, 0.1)
    sim.init_leapfrog()
    p0 = sim.v.sum(axis=1)
    for _ in range(30):
        sim.step()
    assert np.allclose(sim.v.sum(axis=1), p0, atol=1e-9 * sim.np)


def test_boris_gyration_conserves_vperp():
    x = np.zeros((2, 1)) + 1.0
    v = np.array([[0.3], [0.5], [0.0]])
    sim = pic.PIC((4.0, 4.0), (8, 8), x, v, 0.1, B=2.0)
    sim.push(np.zeros((2, 1)))  # yalniz donme
    vp = np.hypot(sim.v[1], sim.v[2])
    assert np.isclose(vp[0], 0.5, rtol=1e-12)
    assert np.isclose(sim.v[0, 0], 0.3)
    ang = np.arctan2(sim.v[2, 0], sim.v[1, 0])
    assert np.isclose(abs(ang), 2 * np.arctan(2.0 * 0.05), rtol=1e-10)


def test_dispersion_roots():
    assert np.isclose(pic.cold_two_stream_growth(np.sqrt(3 / 8)), 1 / (2 * np.sqrt(2)), rtol=1e-10)
    assert pic.cold_two_stream_growth(1.1) <= 1e-12
    w = pic.landau_root(0.5)
    assert abs(w.real - 1.4157) < 2e-4 and abs(w.imag + 0.1534) < 2e-4
    # sicak demetler soguk sinirda soguga yakinsar
    assert np.isclose(pic.warm_two_stream_growth(0.6, 1.0, 0.01), pic.cold_two_stream_growth(0.6), rtol=1e-3)


def test_short_langmuir_frequency():
    Ng, ppc, L, dt = 32, 16, 4 * np.pi, 0.1
    N = Ng * ppc
    x0 = (np.arange(N) + 0.5) * L / N
    k = 2 * np.pi / L
    sim = pic.PIC([L], [Ng], (x0 + 1e-3 / k * np.sin(k * x0))[None], np.zeros((1, N)), dt)
    sim.init_leapfrog()
    s = []
    for _ in range(400):
        o = sim.step(diag=True)
        s.append(np.fft.rfft(o["E"][0])[1].imag)
    s = np.array(s)
    zc = np.where(np.diff(np.sign(s)) != 0)[0]
    w = np.pi / np.mean(np.diff(zc * dt))
    _, wp = pic.langmuir_freq_pred(k, L / Ng, dt)
    assert abs(w / wp - 1) < 5e-3


def test_periodic_labeling_and_events():
    m = np.zeros((10, 6), bool)
    m[0, 2:4] = True
    m[9, 2:4] = True          # sinirin ote yaninda ayni blob
    lab, n = S.label_periodic(m)
    assert n == 1
    a = np.zeros((20,), bool); a[5:15] = True
    b = np.zeros((20,), bool); b[5:9] = True; b[11:15] = True
    la, na = S.label_periodic(a)
    lb, nb = S.label_periodic(b)
    e = S.overlap_events(la, na, lb, nb)
    assert e["splits"] == 1 and e["merges"] == 0
    e2 = S.overlap_events(lb, nb, la, na)
    assert e2["merges"] == 1 and e2["splits"] == 0


def test_persistent_split_vs_flicker():
    one = np.zeros(40); one[10:30] = 1.0
    two = np.zeros(40); two[10:18] = 1.0; two[22:30] = 1.0
    real = np.array([one] * 3 + [two] * 8)            # kalici bolunme
    flick = np.array([one] * 3 + [two] + [one] * 7)   # tek karelik titresim
    _, ev = S.track(real, 0.5, 2, sigma=0.3, P=5)
    assert ev["splits"].sum() == 1 and ev["splits_persistent"].sum() == 1
    _, ev2 = S.track(flick, 0.5, 2, sigma=0.3, P=5)
    assert ev2["splits"].sum() == 1 and ev2["splits_persistent"].sum() == 0


def test_coherence_planar_vs_2d():
    x = np.arange(32) * 0.5
    planar = np.cos(2 * np.pi * x / 16)[:, None] * np.ones((1, 32))
    c, _, _ = S.coherence(planar, (0.5, 0.5))
    assert c > 0.999
    blob = np.cos(2 * np.pi * x / 16)[:, None] * np.cos(2 * np.pi * x / 16)[None, :]
    c2, _, _ = S.coherence(blob, (0.5, 0.5))
    assert c2 < 1e-6
