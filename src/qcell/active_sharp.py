"""Sharp-interface (effective droplet) theory of the chemically active droplets of qcell.active in d dimensions.

Quasi-static limit of ActiveCH in Zwicker et al. (2017) units (w = 1; D = f''(0) = f''(1) = 1; Delta c = 1):
    outside (c ~ 0):  0 = lap c + k (eps - c)          (s = nu+ - k c,  eps = nu+ / k)
    inside  (c ~ 1):  0 = lap c - nu- - k (c - 1)
    interface:        c_out = lg H,  c_in = 1 + lg H    (Gibbs-Thomson; lg = gamma / (f'' Delta c^2) = 1/12 for f = c^2(1-c)^2/2,
                                                        kappa = 1/4; H = sum of the principal curvatures)
                      v_n = d_n c_out - d_n c_in         (mass balance, n = outward normal)
Geometry: an m-dimensional ball of radius R times R^(d-m), i.e. m = d is a droplet, m = d-1 a tube (hypercylinder),
m = 1 in d = 2 a strip. Perturbation r = R + delta Y_l(angles on S^(m-1)) exp(i kz z), kz along the d-m flat
directions. Radial functions are r^-nu I_(nu+l)(Q r) inside and r^-nu K_(nu+l)(Q r) outside, nu = m/2 - 1,
Q^2 = k + kz^2; the curvature perturbation is dH = delta [ (l (l+m-2) - (m-1)) / R^2 + kz^2 ]. m, d may be real.
Checks (tests/test_exp14.py): sigma(l=1, kz=0) = 0 (translation), sigma(l=0, kz=0) = dv/dR, and in the passive limit
the tube's varicose threshold is kz R = sqrt(m-1) = sqrt(d-2) (Rayleigh-Plateau in d dimensions).
"""
import numpy as np
from scipy import optimize, special

LG = 1.0 / 12.0


def _rI(mu, x):
    """I_(mu+1)(x) / I_mu(x)."""
    return special.ive(mu + 1, x) / special.ive(mu, x)


def _rK(mu, x):
    """K_(mu+1)(x) / K_mu(x)."""
    return special.kve(mu + 1, x) / special.kve(mu, x)


def base(R, m, eps, nu_m, k, lg=LG):
    """Radially symmetric quasi-static state of an m-ball of radius R: values at the interface
    (U_out, c_out', c_out'', U_in, c_in', c_in''), U_out = c_out - eps, U_in = c_in - (1 - nu_m / k)."""
    q, nu, H0 = np.sqrt(k), m / 2.0 - 1.0, (m - 1.0) / R
    uo, ui = lg * H0 - eps, lg * H0 + nu_m / k
    dco, dci = -q * _rK(nu, q * R) * uo, q * _rI(nu, q * R) * ui
    return uo, dco, k * uo - (m - 1.0) / R * dco, ui, dci, k * ui - (m - 1.0) / R * dci


def velocity(R, m, eps, nu_m, k, lg=LG):
    """Interface velocity dR/dt of the unperturbed m-ball (tube for m = d-1)."""
    b = base(R, m, eps, nu_m, k, lg)
    return b[1] - b[4]


def sigma(R, m, l, kz, eps, nu_m, k, lg=LG):
    """Linear growth rate of the shape mode (l, kz) of an m-ball / tube of radius R (exact when velocity(R) = 0)."""
    Q, mu = np.sqrt(k + kz**2), m / 2.0 - 1.0 + l
    _, dco, d2co, _, dci, d2ci = base(R, m, eps, nu_m, k, lg)
    dH = (l * (l + m - 2.0) - (m - 1.0)) / R**2 + kz**2
    a, b = lg * dH - dco, lg * dH - dci
    return d2co - d2ci + a * (l / R - Q * _rK(mu, Q * R)) - b * (l / R + Q * _rI(mu, Q * R))


def radii(m, eps, nu_m, k, lg=LG, rmin=0.3, rmax=400.0, n=4000):
    """All roots of velocity(R) = 0 on [rmin, rmax] as (R, stable) with stable = dv/dR < 0."""
    R = np.geomspace(rmin, rmax, n)
    v = velocity(R, m, eps, nu_m, k, lg)
    out = []
    for i in np.nonzero(np.sign(v[:-1]) != np.sign(v[1:]))[0]:
        r = optimize.brentq(velocity, R[i], R[i + 1], args=(m, eps, nu_m, k, lg), xtol=1e-12)
        out.append((r, bool(v[i] > 0 > v[i + 1])))
    return out


def stable_radius(m, eps, nu_m, k, lg=LG):
    s = [r for r, st in radii(m, eps, nu_m, k, lg) if st]
    return s[0] if s else np.nan


def shape_threshold(d, l, eps, nu_m, k, lg=LG, rmax=400.0):
    """Smallest R at which the droplet mode l becomes unstable (sigma_l(R) = 0 crossing from below), or nan."""
    R = np.geomspace(0.5, rmax, 3000)
    s = sigma(R, d, l, 0.0, eps, nu_m, k, lg)
    i = np.nonzero((s[:-1] < 0) & (s[1:] > 0))[0]
    if len(i) == 0:
        return np.nan
    return optimize.brentq(lambda r: sigma(r, d, l, 0.0, eps, nu_m, k, lg), R[i[0]], R[i[0] + 1], xtol=1e-12)


def tube_spectrum(d, eps, nu_m, k, lg=LG, l=0, kz=None, R=None):
    """Growth rates sigma(kz) of a tube (m = d-1) at its stable radius (or at R). Returns (R, kz, sigma)."""
    R = stable_radius(d - 1.0, eps, nu_m, k, lg) if R is None else R
    kz = np.linspace(1e-4, 2.0, 2000) if kz is None else np.asarray(kz)
    return R, kz, sigma(R, d - 1.0, l, kz, eps, nu_m, k, lg)
