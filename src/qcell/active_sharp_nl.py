"""Nonlinear sharp-interface theory of the active droplets of qcell.active in d dimensions (chemical-potential form).

The linear effective-droplet model (active_sharp) uses D = f''(0) = 1 in the dilute phase. For eps ~ 0.2 the far field
sits next to the spinodal (c_sp = 0.211, f''(0.2) = 0.04), so it overestimates the supply and R* by ~30 % (checked
against the 2D phase-field droplets of exp11). Here each bulk phase is solved in the chemical potential mu = f'(c):
    0 = lap mu + s(c(mu))         c(mu): inverse of f'(c) = 2c^3 - 3c^2 + c on the dilute (c < c_sp1) / dense (c > c_sp2)
                                  branch; s linear on both branches (s = nu+ - k c, s = -nu- - k (c-1))
    interface: mu = lg H (lg = gamma / Delta c = 1/12),  Delta c v_n = d_n mu_out - d_n mu_in
Geometry and perturbations as in active_sharp (m-ball x R^(d-m)). Perturbations obey
    dmu'' + (m-1)/r dmu' - [l(l+m-2)/r^2 + kz^2 + k / f''(c0(r))] dmu = 0,
integrated numerically (regular inside, decaying outside). If the interior would have to go below the dense spinodal
(mu < mu_sp_dense = -0.0962) no dense core exists: the droplet hollows (liquid shell) - returned as nan / flagged.
"""
import numpy as np
from scipy import integrate, optimize, special

from qcell import active_sharp as lin

LG = 1.0 / 12.0
C_SP1, C_SP2 = (3 - np.sqrt(3)) / 6, (3 + np.sqrt(3)) / 6


def fp(c):
    return c * (c - 1.0) * (2.0 * c - 1.0)


def fpp(c):
    return 6.0 * c * c - 6.0 * c + 1.0


MU_SP = float(fp(C_SP1))            # +0.0962: upper end of the dilute branch; -MU_SP: lower end of the dense branch
_cd = np.linspace(-0.6, C_SP1, 200001)
_ce = np.linspace(C_SP2, 1.6, 200001)
_md, _me = fp(_cd), fp(_ce)          # monotone tables mu(c) on each branch (computed once)


def c_dilute(mu):
    return np.interp(mu, _md, _cd)


def c_dense(mu):
    return np.interp(mu, _me, _ce)


class Chem:
    def __init__(self, eps, nu_m, k):
        self.eps, self.nu_p, self.nu_m, self.k = eps, eps * k, nu_m, k
        self.mu_inf = float(fp(eps))

    def s_out(self, mu):
        return self.nu_p - self.k * c_dilute(mu)

    def s_in(self, mu):
        return -self.nu_m - self.k * (c_dense(mu) - 1.0)


RTOL = 1e-8


def _inside(R, m, ch, mu_I):
    """Regular interior solution with mu(R) = mu_I: returns (mu'(R), mu(0), sol) or None if no dense core exists."""
    def shoot(mu0, dense_out=False):
        r0 = 1e-6 * R
        a = -float(ch.s_in(mu0)) / m
        f = lambda r, y: [y[1], -(m - 1.0) / r * y[1] - float(ch.s_in(y[0]))]
        return integrate.solve_ivp(f, [r0, R], [mu0 + 0.5 * a * r0**2, a * r0], rtol=RTOL, atol=1e-12,
                                   dense_output=dense_out)

    g = lambda mu0: shoot(mu0).y[0, -1] - mu_I
    lo, hi = -MU_SP + 1e-9, mu_I + 0.5
    if g(lo) > 0:
        return None                                   # even the lowest admissible core value overshoots: hollow
    mu0 = optimize.brentq(g, lo, hi, xtol=1e-12)
    sol = shoot(mu0, True)
    return float(sol.y[1, -1]), float(mu0), sol


def _outside(R, m, ch, mu_I):
    """Outer solution with mu(R) = mu_I and mu -> mu_inf, by inward shooting from r_f = R + span: the far field is
    mu_inf + delta * (decaying modified-Bessel mode), delta < 0 found by root finding in log|delta|.
    Returns (mu'(R), dense solution on [R, r_f])."""
    kap = np.sqrt(ch.k / max(float(fpp(ch.eps)), 1e-6))
    span = 12.0 / kap + 4.0 / np.sqrt(ch.k)
    rf, nu = R + span, m / 2.0 - 1.0
    gK = -kap * special.kve(nu + 1, kap * rf) / special.kve(nu, kap * rf)       # d ln(r^-nu K_nu(kap r)) / dr
    f = lambda r, y: [y[1], -(m - 1.0) / r * y[1] - float(ch.s_out(min(y[0], MU_SP - 1e-12)))]

    def shoot(ld, dense_out=False):
        dlt = -np.exp(ld)
        return integrate.solve_ivp(f, [rf, R], [ch.mu_inf + dlt, gK * dlt], rtol=RTOL, atol=1e-14,
                                   dense_output=dense_out)

    g = lambda ld: shoot(ld).y[0, -1] - mu_I
    lo, hi = -60.0, np.log(max(ch.mu_inf - mu_I, 1e-12))
    if g(lo) * g(hi) > 0:
        return None
    ld = optimize.brentq(g, lo, hi, xtol=1e-10)
    sol = shoot(ld, True)
    return float(sol.y[1, -1]), sol


def base(R, m, ch, lg=LG):
    mu_I = lg * (m - 1.0) / R
    ins = _inside(R, m, ch, mu_I)
    if ins is None:
        return None
    dmi, mu0, si = ins
    out = _outside(R, m, ch, mu_I)
    if out is None:
        return None
    dmo, so = out
    dc = float(c_dense(mu_I) - c_dilute(mu_I))
    return dict(mu_I=mu_I, dmi=dmi, dmo=dmo, mu0=mu0, c0=float(c_dense(mu0)), dc=dc, sin=si, sout=so,
                d2mi=-(m - 1.0) / R * dmi - float(ch.s_in(mu_I)), d2mo=-(m - 1.0) / R * dmo - float(ch.s_out(mu_I)))


def velocity(R, m, ch, lg=LG):
    b = base(R, m, ch, lg)
    return np.nan if b is None else (b["dmo"] - b["dmi"]) / b["dc"]


def stable_radius(m, ch, lg=LG, guess=None):
    """Stable stationary radius (v > 0 below, < 0 above), bracketed around the linear-theory value."""
    g = lin.stable_radius(m, ch.eps, ch.nu_m, ch.k, lg) if guess is None else guess
    if not np.isfinite(g):
        return np.nan
    R = np.array([0.3, 0.45, 0.6, 0.75, 0.9, 1.1, 1.4]) * g
    v = [velocity(r, m, ch, lg) for r in R]
    for i in range(len(R) - 1):
        if np.isfinite(v[i]) and np.isfinite(v[i + 1]) and v[i] > 0 > v[i + 1]:
            return optimize.brentq(lambda r: velocity(r, m, ch, lg), R[i], R[i + 1], xtol=1e-4)
    return np.nan


def sigma(R, m, l, kz, ch, lg=LG, b=None):
    """Growth rate of the shape mode (l, kz) of the m-ball / tube of radius R (stationary state assumed)."""
    b = base(R, m, ch, lg) if b is None else b
    if b is None:
        return np.nan
    lam = l * (l + m - 2.0)
    # inside: regular solution ~ r^l
    si = b["sin"]

    def fin(r, y):
        c0 = c_dense(si.sol(r)[0])
        return [y[1], -(m - 1.0) / r * y[1] + (lam / r**2 + kz**2 + ch.k / fpp(c0)) * y[0]]
    r0 = 1e-3 * R
    yi = integrate.solve_ivp(fin, [r0, R], [r0**l, l * r0 ** (l - 1) if l > 0 else 0.0], rtol=1e-9, atol=1e-14)
    Fi = yi.y[1, -1] / yi.y[0, -1]
    # outside: decaying solution integrated inward from the far end of the base solution
    so = b["sout"]
    rf = so.t[0]
    Qf = np.sqrt(kz**2 + ch.k / fpp(ch.eps))

    def fout(r, y):
        c0 = c_dilute(min(so.sol(r)[0], MU_SP - 1e-12))
        return [y[1], -(m - 1.0) / r * y[1] + (lam / r**2 + kz**2 + ch.k / fpp(c0)) * y[0]]
    yo = integrate.solve_ivp(fout, [rf, R], [1e-8, -Qf * 1e-8], rtol=1e-9, atol=1e-30)
    Fo = yo.y[1, -1] / yo.y[0, -1]
    dH = (lam - (m - 1.0)) / R**2 + kz**2
    a, bb = lg * dH - b["dmo"], lg * dH - b["dmi"]
    return float((b["d2mo"] - b["d2mi"] + a * Fo - bb * Fi) / b["dc"])
