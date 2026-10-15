"""Charged Lennard-Jones clusters with gate-screened Coulomb repulsion (classical particles, reduced units eps=sigma=1).

E = 4 sum (r^-12 - r^-6) + g sum ( 1/r - 1/sqrt(r^2 + 4 d^2) )      (d = inf -> plain 1/r)
g = q^2 / (4 pi eps sigma eps_LJ): Coulomb strength; d = gate distance (grounded gate image, in units of sigma).
"""
import numpy as np
from scipy.optimize import basinhopping, minimize


def energy_grad(x, n, g, d):
    X = x.reshape(n, 3)
    diff = X[:, None, :] - X[None, :, :]
    r2 = np.sum(diff**2, axis=-1)
    iu = np.triu_indices(n, 1)
    r2u = np.maximum(r2[iu], 1e-12)
    inv6 = r2u**-3
    e = 4 * np.sum(inv6**2 - inv6)
    dedr2 = 4 * (-6 * inv6**2 / r2u + 3 * inv6 / r2u)  # d/d(r^2) of 4(r^-12 - r^-6)
    if g != 0.0:
        r = np.sqrt(r2u)
        if np.isinf(d):
            e += g * np.sum(1 / r)
            dedr2 += g * (-0.5) * r**-3
        else:
            s = np.sqrt(r2u + 4 * d * d)
            e += g * np.sum(1 / r - 1 / s)
            dedr2 += g * (-0.5 * r**-3 + 0.5 * s**-3)
    grad = np.zeros_like(X)
    coef = np.zeros((n, n))
    coef[iu] = 2 * dedr2
    coef = coef + coef.T
    grad = np.sum(coef[:, :, None] * diff, axis=1)
    return e, grad.ravel()


def _fun(x, n, g, d):
    return energy_grad(x, n, g, d)


def minimise(n, g=0.0, d=np.inf, niter=60, seed=0, x0=None, T=1.0, stepsize=0.4):
    """Basin-hopping global minimisation (best of random restarts + hopping)."""
    rng = np.random.default_rng(seed)
    if n == 1:
        return 0.0, np.zeros((1, 3))
    best = (np.inf, None)
    for rep in range(2):
        if x0 is None or rep > 0:
            R = 0.7 * n ** (1 / 3) + 0.5
            pts = rng.normal(size=(n, 3))
            pts = pts / np.linalg.norm(pts, axis=1)[:, None] * (rng.random(n) ** (1 / 3))[:, None] * R
            start = pts.ravel()
        else:
            start = x0.ravel()
        res = basinhopping(lambda x: _fun(x, n, g, d)[0], start, niter=niter, T=T, stepsize=stepsize,
                           minimizer_kwargs=dict(method="L-BFGS-B", jac=lambda x: _fun(x, n, g, d)[1]), seed=int(rng.integers(1 << 30)))
        if res.fun < best[0]:
            best = (res.fun, res.x.reshape(n, 3))
    return float(best[0]), best[1]
