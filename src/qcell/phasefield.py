"""Conserved phase-field model with a gate-screened Coulomb (Ohta-Kawasaki type) frustration term.

d_t phi = M lap( mu ),   mu = -phi + phi^3 - kappa lap(phi) + F^-1[ g(q) phi_hat ]
g(q) = A (1 - exp(-2 q d)) / q     (2D Coulomb kernel A/q screened by a grounded gate at distance d; d = inf -> A/q)
Free energy:  F = int [ -phi^2/2 + phi^4/4 + kappa/2 |grad phi|^2 ] + 1/2 int phi G phi.
Linear growth rate about a uniform state phi0:  sigma(q) = - M q^2 [ -1 + 3 phi0^2 + kappa q^2 + g(q) ].
Semi-implicit (linearly stabilised) Fourier pseudo-spectral time stepping; total phi is conserved exactly (q=0 mode untouched).
"""
import numpy as np


class OKField:
    def __init__(self, nx=128, ny=128, Lx=100.0, Ly=100.0, kappa=1.0, A=0.1, d=np.inf, M=1.0, S=2.0):
        self.nx, self.ny, self.Lx, self.Ly = nx, ny, Lx, Ly
        self.kappa, self.A, self.d, self.M, self.S = kappa, A, d, M, S
        qx = 2 * np.pi * np.fft.fftfreq(nx, d=Lx / nx)
        qy = 2 * np.pi * np.fft.fftfreq(ny, d=Ly / ny)
        self.qx, self.qy = qx, qy
        QX, QY = np.meshgrid(qx, qy, indexing="ij")
        self.q2 = QX**2 + QY**2
        self.q = np.sqrt(self.q2)
        self.g = self.kernel(self.q)

    def kernel(self, q):
        q = np.asarray(q, float)
        out = np.zeros_like(q)
        nz = q > 0
        if np.isinf(self.d):
            out[nz] = self.A / q[nz]
        else:
            out[nz] = self.A * (-np.expm1(-2 * q[nz] * self.d)) / q[nz]
            out[~nz] = 2 * self.A * self.d
        return out

    def growth_rate(self, q, phi0=0.0):
        q = np.asarray(q, float)
        return -self.M * q**2 * (-1 + 3 * phi0**2 + self.kappa * q**2 + self.kernel(q))

    def step(self, phi, dt):
        nl = -phi + phi**3
        ph = np.fft.fft2(phi)
        nlh = np.fft.fft2(nl)
        lin = self.kappa * self.q2 + self.g
        ph_new = (ph - self.M * dt * self.q2 * (nlh - self.S * ph)) / (1 + self.M * dt * self.q2 * (lin + self.S))
        return np.real(np.fft.ifft2(ph_new))

    def run(self, phi, t_end, dt=0.05, snap_every=None):
        n = int(round(t_end / dt))
        snaps = []
        for i in range(n):
            phi = self.step(phi, dt)
            if snap_every and (i + 1) % snap_every == 0:
                snaps.append((i + 1) * dt)
                snaps.append(phi.copy())
        return phi, snaps

    def structure_peak(self, phi):
        """Mean wavenumber weighted by the structure factor (excluding q=0)."""
        s = np.abs(np.fft.fft2(phi - phi.mean())) ** 2
        s = s * (self.q > 0)
        return float((s * self.q).sum() / s.sum())

    def critical_A(self, phi0=0.0):
        """Smallest A... returns the largest A for which the uniform state phi0 is unstable (numerical, over a fine q grid)."""
        qs = np.linspace(1e-3, 3.0, 20000)
        base = -1 + 3 * phi0**2
        lo, hi = 0.0, 5.0
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            test = OKField.__new__(OKField)
            test.A, test.d = mid, self.d
            test.kappa = self.kappa
            if np.min(base + self.kappa * qs**2 + OKField.kernel(test, qs)) < 0:
                lo = mid
            else:
                hi = mid
        return 0.5 * (lo + hi)
