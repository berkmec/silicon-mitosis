import numpy as np

from qcell import cluster
from qcell.phasefield import OKField


def test_cluster_gradient_matches_finite_difference():
    rng = np.random.default_rng(0)
    n = 6
    x = rng.normal(size=3 * n) * 1.2
    for g, d in ((0.0, np.inf), (0.3, np.inf), (0.3, 2.0)):
        e, gr = cluster.energy_grad(x, n, g, d)
        num = np.zeros_like(x)
        for i in range(len(x)):
            h = 1e-6
            xp, xm = x.copy(), x.copy()
            xp[i] += h
            xm[i] -= h
            num[i] = (cluster.energy_grad(xp, n, g, d)[0] - cluster.energy_grad(xm, n, g, d)[0]) / (2 * h)
        assert np.allclose(gr, num, rtol=1e-5, atol=1e-6)


def test_lj_global_minima():
    for n, ref in ((7, -16.505384), (13, -44.326801)):
        e, _ = cluster.minimise(n, niter=40, seed=1)
        assert abs(e - ref) < 1e-4


def test_phasefield_conserves_mass_and_linear_growth():
    f = OKField(nx=64, ny=4, Lx=100.0, Ly=1.0, A=0.1)
    x = np.arange(64) * 100.0 / 64
    phi = np.tile((1e-7 * np.cos(2 * np.pi * 5 * x / 100.0))[:, None], (1, 4))
    out, _ = f.run(phi.copy(), t_end=4.0, dt=0.002)
    assert abs(out.sum() - phi.sum()) < 1e-12
    q = 2 * np.pi * 5 / 100.0
    sig = np.log(abs(np.fft.fft(out[:, 0])[5]) / abs(np.fft.fft(phi[:, 0])[5])) / 4.0
    assert abs(sig - f.growth_rate(q)) / abs(f.growth_rate(q)) < 5e-3
