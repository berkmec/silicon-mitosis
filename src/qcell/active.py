"""Unified phase-field model for the "which drives division harder?" comparison (exp11).

Base: Zwicker, Seyboldt, Weber, Hyman, Juelicher, Nat. Phys. 13, 408 (2017), supplement Eqs. (S.1)-(S.5), in their units
(concentration / Delta c with c+(0) = 0, c-(0) = 1; length / w, w = 2 sqrt(kappa/b); time / t0 = w^2/D, D = m b):

    d_t c = lap[ c (c-1)(2c-1) ] - (1/4) lap^2 c  +  s(c)            (driven chemistry / charging, s != 0)
            + lap[ G * (c - <c>) ]                                     (gate-screened Coulomb, A != 0; 2D only)
    i.e. d_t c = lap[ mu + G * (c - <c>) ] + s(c): a repulsive long-range term (review 2026-10-11: the sign was written
    as '-' here before; the code always implemented '+', see `lin` in __init__).

s(c) = nu+ - k+ c                     for c < c_c+ (= 0.25)
     = -nu- - k- (c - 1)              for c > c_c- (= 0.75)
     = cubic Hermite in between (value and slope continuous), as in the paper.
G(q) = A (1 - exp(-2 q d)) / q  (2D Coulomb screened by a grounded gate at distance d; d = inf -> A / q), G(0) = 0.

Note: with k+ = k- = k and nu- = k - nu+ (dnu = nu- + nu+ - k = 0) s(c) = nu+ - k c is one straight line, and the
dynamics is the H^-1 gradient flow of an Ohta-Kawasaki functional with kernel k / q^2 (Glotzer et al. 1995 mapping).
dnu != 0 is the part that cannot be written as such a gradient flow.
Periodic box, semi-implicit Fourier scheme (stiff linear part implicit, stabiliser S), real FFTs.
"""
import numpy as np
from scipy import fft as sfft
from scipy import ndimage


def rates_zwicker(eps, dnu, k=1e-2):
    """Linear chemistry with supersaturation eps = nu+/k and excess turnover dnu = nu- + nu+ - k (Delta c = 1)."""
    nu_p = eps * k
    return dict(nu_p=nu_p, nu_m=dnu + k - nu_p, k_p=k, k_m=k)


def rates_charging(J, lam=0.0):
    """Uniform charge injection J with linear leakage lam * c, i.e. s = J - lam c (lies on the dnu = 0 line)."""
    return dict(nu_p=J, nu_m=lam - J, k_p=lam, k_m=lam)


class ActiveCH:
    def __init__(self, n=128, L=64.0, dim=2, nu_p=0.0, nu_m=0.0, k_p=0.0, k_m=0.0, cc_p=0.25, cc_m=0.75,
                 A=0.0, d=np.inf, S=2.0, workers=1, dtype=np.float64):
        """Defaults are PASSIVE (no chemistry, no Coulomb). dtype=np.float32 runs the time loop in single precision."""
        self.dtype = np.dtype(dtype).type
        if A != 0.0 and dim != 2:
            raise ValueError("Coulomb kernel A/q is the 2D-layer kernel; use dim=2")
        self.n, self.L, self.dim, self.dx = n, L, dim, L / n
        self.nu_p, self.nu_m, self.k_p, self.k_m = nu_p, nu_m, k_p, k_m
        self.cc_p, self.cc_m = cc_p, cc_m
        self.A, self.d, self.S, self.workers = A, d, S, workers
        freqs = [2 * np.pi * np.fft.fftfreq(n, d=self.dx)] * (dim - 1) + [2 * np.pi * np.fft.rfftfreq(n, d=self.dx)]
        grids = np.meshgrid(*freqs, indexing="ij")
        self.q2 = sum(g**2 for g in grids)
        q = np.sqrt(self.q2)
        g = np.zeros_like(q)
        nz = q > 0
        if A != 0.0:
            g[nz] = A / q[nz] if np.isinf(d) else A * (-np.expm1(-2 * q[nz] * d)) / q[nz]
        self.g = g
        self.lin = self.q2 * (S + g) + 0.25 * self.q2**2
        self.active = any(v != 0.0 for v in (nu_p, nu_m, k_p, k_m))
        self._cubic = self._hermite()

    # ---------------------------------------------------------------- chemistry
    def _hermite(self):
        x0, x1 = self.cc_p, self.cc_m
        y0, y1 = self.nu_p - self.k_p * x0, -self.nu_m - self.k_m * (x1 - 1.0)
        m0, m1 = -self.k_p, -self.k_m
        M = np.array([[1, x0, x0**2, x0**3], [1, x1, x1**2, x1**3], [0, 1, 2 * x0, 3 * x0**2], [0, 1, 2 * x1, 3 * x1**2]])
        return np.linalg.solve(M, np.array([y0, y1, m0, m1]))

    def s(self, c):
        ty = c.dtype.type  # keep the caller's precision (no float64 upcast of float32 fields)
        a = self._cubic.astype(c.dtype)
        out = a[0] + c * (a[1] + c * (a[2] + c * a[3]))
        out = np.where(c < ty(self.cc_p), ty(self.nu_p) - ty(self.k_p) * c, out)
        return np.where(c > ty(self.cc_m), -ty(self.nu_m) - ty(self.k_m) * (c - ty(1.0)), out)

    # ---------------------------------------------------------------- dynamics
    def run(self, c, T, dt, every=None, cb=None):
        """Integrate to time T. cb(t, c) is called every `every` steps; returning False stops the run."""
        w, ty = self.workers, self.dtype
        c = np.asarray(c, dtype=ty)
        q2, den, S, h = self.q2.astype(ty), (1.0 + dt * self.lin).astype(ty), ty(self.S), ty(dt)
        one, two = ty(1.0), ty(2.0)
        ch = sfft.rfftn(c, workers=w)
        for i in range(1, int(round(T / dt)) + 1):
            nl = c * (c - one) * (two * c - one)
            rhs = ch - h * q2 * (sfft.rfftn(nl, workers=w) - S * ch)
            if self.active:
                rhs += h * sfft.rfftn(self.s(c), workers=w)
            ch = rhs / den
            c = sfft.irfftn(ch, s=c.shape, workers=w)
            if cb is not None and i % every == 0 and cb(i * dt, c) is False:
                break
        return c

    def step(self, c, dt):
        return self.run(c, dt, dt)

    # ---------------------------------------------------------------- analysis
    def droplets(self, c, thr=0.5, r_min=1.5):
        """Droplets (c > thr) on the periodic grid. Returns dict(n, radii, perc, aspect); perc = a component wraps
        the box; aspect = sqrt(largest/smallest eigenvalue of the second-moment tensor) of the largest droplet."""
        n, dim = self.n, self.dim
        big = np.tile(c > thr, (3,) * dim)
        lab, k = ndimage.label(big)
        if k == 0:
            return dict(n=0, radii=[], perc=False, aspect=1.0)
        idx = range(1, k + 1)
        com = ndimage.center_of_mass(big, lab, idx)
        size = ndimage.sum(big, lab, idx)
        slices = ndimage.find_objects(lab)
        radii, perc, best = [], False, (0.0, None)
        for i in range(k):
            if not all(n <= x < 2 * n for x in com[i]):
                continue
            if any(sl.stop - sl.start >= 2 * n for sl in slices[i]):
                perc = True
                continue
            vol = size[i] * self.dx**dim
            r = float(np.sqrt(vol / np.pi)) if dim == 2 else float((3 * vol / (4 * np.pi)) ** (1 / 3))
            if r >= r_min:
                radii.append(r)
                if r > best[0]:
                    best = (r, i)
        aspect = 1.0
        if best[1] is not None:
            pts = np.argwhere(lab[slices[best[1]]] == best[1] + 1).astype(float)
            ev = np.linalg.eigvalsh(np.cov(pts.T))
            aspect = float(np.sqrt(ev[-1] / max(ev[0], 1e-12)))
        return dict(n=len(radii), radii=radii, perc=perc, aspect=aspect)

    def disk(self, R, center=None, c_out=None, ecc=0.0):
        """Droplet with a tanh interface of width w = 1; ecc > 0 gives semi-axes R(1+ecc), R(1-ecc) (2D) or a
        prolate spheroid (3D)."""
        n, dx, dim = self.n, self.dx, self.dim
        x = (np.arange(n) + 0.5) * dx
        X = np.meshgrid(*([x] * dim), indexing="ij")
        ctr = [self.L / 2] * dim if center is None else center
        a = [R * (1 + ecc)] + [R * (1 - ecc)] * (dim - 1)
        L = self.L  # minimum-image (periodic) offsets
        rho = np.sqrt(sum((((Xi - ci + L / 2) % L - L / 2) / ai) ** 2 for Xi, ci, ai in zip(X, ctr, a)))
        dist = (rho - 1.0) * R
        c_far = (self.nu_p / self.k_p if self.k_p > 0 else 0.0) if c_out is None else c_out
        return c_far + (1.0 - c_far) * 0.5 * (1 - np.tanh(dist / 1.0))
