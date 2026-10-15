"""Axisymmetric reduction of the active Cahn-Hilliard model (qcell.active) to the (z, r) half-plane in d dimensions.

For shapes with O(d-1) symmetry about the z axis the d-dimensional equation
    d_t c = lap[ c(c-1)(2c-1) - lap c / 4 ] + s(c),    lap = d_zz + d_rr + (d-2)/r d_r
is exact on the half-plane, and d may be any real number >= 2 (dimensional continuation). d = 2 is the mirror-symmetric
(varicose-only) 2D problem.
Discretisation: z periodic (Fourier), r cell-centred finite volumes on [0, r_max] with face areas r^(d-2) and exact cell
volumes (r_+^(d-1) - r_-^(d-1)) / (d-1), no flux through the axis and the outer wall, so mass sum(V c) is conserved
exactly by the transport part. Semi-implicit stabilised step as in ActiveCH:
    (1 - dt S L_k + dt/4 L_k^2) c^(n+1) = c^n + dt L_k (f'(c^n) - S c^n) + dt s(c^n),   L_k = L_r - kz^2,
solved per z-Fourier mode with precomputed dense inverses.
"""
import numpy as np
from scipy import ndimage, special

from qcell.active import ActiveCH


def ball_volume(d):
    return np.pi ** (d / 2) / special.gamma(d / 2 + 1)


def sphere_area(n):
    """Area of the unit n-sphere S^n (n = d-2 for the transverse sphere)."""
    return 2 * np.pi ** ((n + 1) / 2) / special.gamma((n + 1) / 2)


class AxiCH:
    def __init__(self, d, nz, nr, h=0.5, dt=0.2, S=1.0, nu_p=0.0, nu_m=0.0, k_p=0.0, k_m=0.0):
        self.d, self.nz, self.nr, self.h, self.dt, self.S = float(d), nz, nr, h, dt, S
        self.Lz, self.rmax = nz * h, nr * h
        rf = np.arange(nr + 1) * h
        self.r = (np.arange(nr) + 0.5) * h
        self.z = (np.arange(nz) + 0.5) * h
        A = rf ** (self.d - 2.0)
        A[0] = A[-1] = 0.0                                     # axis (regularity / mirror) and outer wall: no flux
        self.V = (rf[1:] ** (self.d - 1.0) - rf[:-1] ** (self.d - 1.0)) / (self.d - 1.0)
        lo, up = A[:-1] / (h * self.V), A[1:] / (h * self.V)
        Lr = np.diag(-(lo + up)) + np.diag(lo[1:], -1) + np.diag(up[:-1], 1)
        self.Lr = Lr
        self.kz = 2 * np.pi * np.fft.rfftfreq(nz, d=h)
        I = np.eye(nr)
        Lk = Lr[None] - (self.kz**2)[:, None, None] * I[None]
        self.Lk = Lk
        self.Minv = np.linalg.inv(I[None] - dt * S * Lk + 0.25 * dt * Lk @ Lk)
        self.chem = ActiveCH(n=4, L=4.0, dim=1, nu_p=nu_p, nu_m=nu_m, k_p=k_p, k_m=k_m)   # reuse s(c) only
        self.active = self.chem.active
        self.t = 0.0

    # ------------------------------------------------------------------ dynamics
    def step(self, c, n=1):
        dt, S = self.dt, self.S
        for _ in range(n):
            N = c * (c - 1.0) * (2.0 * c - 1.0) - S * c
            Nh = np.fft.rfft(N, axis=0)
            rhs = np.fft.rfft(c, axis=0) + dt * (Nh @ self.Lr.T - (self.kz**2)[:, None] * Nh)
            if self.active:
                rhs += dt * np.fft.rfft(self.chem.s(c), axis=0)
            c = np.fft.irfft(np.einsum("kij,kj->ki", self.Minv, rhs), n=self.nz, axis=0)
            self.t += dt
        return c

    # ------------------------------------------------------------------ initial states
    def droplet(self, R, ecc=0.0, zc=None, c_far=None):
        zc = self.Lz / 2 if zc is None else zc
        Z, Rr = np.meshgrid(self.z, self.r, indexing="ij")
        rho = np.sqrt(((Z - zc) / (R * (1 + ecc))) ** 2 + (Rr / (R * (1 - ecc))) ** 2)
        cf = self._cfar() if c_far is None else c_far
        return cf + (1.0 - cf) * 0.5 * (1 - np.tanh((rho - 1.0) * R))

    def tube(self, R, amp=0.0, m=1, c_far=None):
        """Tube of radius R(z) = R (1 + amp cos(2 pi m z / Lz)) spanning the periodic z box."""
        Z, Rr = np.meshgrid(self.z, self.r, indexing="ij")
        Rz = R * (1 + amp * np.cos(2 * np.pi * m * Z / self.Lz))
        cf = self._cfar() if c_far is None else c_far
        return cf + (1.0 - cf) * 0.5 * (1 - np.tanh(Rr - Rz))

    def _cfar(self):
        return self.chem.nu_p / self.chem.k_p if self.chem.k_p > 0 else 0.0

    # ------------------------------------------------------------------ analysis
    def mass(self, c):
        return float((c * self.V[None]).sum() * self.h * sphere_area(self.d - 2.0))

    def bodies(self, c, thr=0.5):
        """Connected dense regions on the periodic half-plane. Each: dict(axis=touches r=0 (a body on the axis; else a
        ring / hollow shell), vol (d-dim volume), R_eq, zc, length, rmax, neck = min radius over the interior z range)."""
        b = c > thr
        nz = self.nz
        big = np.concatenate([b, b, b], axis=0)
        lab, k = ndimage.label(big)
        out = []
        if k == 0:
            return out
        sl = ndimage.find_objects(lab)
        for i in range(k):
            zs = sl[i][0]
            if zs.stop - zs.start >= 3 * nz:          # wraps the box: a tube spanning z
                if not any(o.get("span") for o in out):
                    msk = lab[nz:2 * nz] == i + 1
                    out.append(self._describe(msk, span=True))
                continue
            ctr = 0.5 * (zs.start + zs.stop - 1)
            if not (nz <= ctr < 2 * nz):
                continue
            msk = lab[zs] == i + 1
            out.append(self._describe(msk, span=False, z0=zs.start - nz))
        return out

    def bubbles(self, c, thr=0.5):
        """Number of enclosed dilute regions (c < thr not connected to the outer wall r = r_max): a hollow droplet /
        liquid shell has one."""
        lab, k = ndimage.label(np.concatenate([c, c, c], axis=0) < thr)
        if k == 0:
            return 0
        outer = set(np.unique(lab[:, -1])) - {0}
        mid = lab[self.nz:2 * self.nz]
        return len(set(np.unique(mid)) - {0} - outer)

    def _describe(self, msk, span, z0=0):
        h, d = self.h, self.d
        vol = float((msk * self.V[None]).sum() * h * sphere_area(d - 2.0))
        axis = bool(msk[:, 0].any())
        rows = np.nonzero(msk.any(axis=1))[0]
        prof = np.array([(np.nonzero(msk[j])[0].max() + 1) * h if msk[j].any() else 0.0 for j in range(msk.shape[0])])
        inner = prof[rows[0]:rows[-1] + 1]
        neck = float(inner[len(inner) // 8: len(inner) - len(inner) // 8].min()) if len(inner) >= 8 else float(inner.min())
        return dict(axis=axis, span=span, vol=vol, R_eq=float((vol / ball_volume(d)) ** (1 / d)),
                    zc=float(((rows.mean() + z0) % self.nz + 0.5) * h), length=float(len(rows) * h),
                    rmax=float(prof.max()), neck=neck)
