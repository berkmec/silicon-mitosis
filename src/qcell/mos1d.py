"""1D self-consistent Poisson-Schroedinger solver for an n+poly / SiO2 / p-Si MOS stack (electron inversion layer).

Units: z in nm inside the Schroedinger solver, SI in Poisson. Energies in eV relative to the Fermi level.
Valleys: unprimed (m_z=0.916, g=2, m_d=0.19), primed (m_z=0.19, g=4, m_d=0.417); spin 2 included in the 2D DOS.
Approximations: infinite barrier at the oxide interface, no poly depletion, no interface charge, T=300 K, Boltzmann holes.
"""
import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve

Q = 1.602176634e-19
KB = 1.380649e-23
EPS0 = 8.8541878128e-12
HBAR = 1.054571817e-34
M0 = 9.1093837015e-31
HB_EVNM2 = 0.0380998212  # hbar^2/(2 m0) in eV nm^2
NI = 1.0e16  # m^-3 (1e10 cm^-3)
EPS_SI, EPS_OX = 11.7, 3.9
EG = 1.12
VALLEYS = [dict(mz=0.916, g=2, md=0.19), dict(mz=0.19, g=4, md=0.417)]


class Mos:
    def __init__(self, na_cm3=1.7e17, tox_nm=4.148, zmax_nm=140.0, dz_nm=0.1, T=300.0):
        self.na = na_cm3 * 1e6
        self.tox = tox_nm * 1e-9
        self.T = T
        self.vt = KB * T / Q
        self.dz = dz_nm
        self.z = np.arange(0.0, zmax_nm + 1e-9, dz_nm)
        self.phif = self.vt * np.log(self.na / NI)
        self.vfb = -(EG / 2 + self.phif)  # n+ poly
        self.cox = EPS_OX * EPS0 / self.tox

    # --- Schroedinger ---
    def _subbands(self, psi, nlev=6):
        """Return list of (E_j eV rel. E_F, |chi_j|^2 in 1/nm, g*md, ) for each valley."""
        Ec = EG / 2 + self.phif - psi  # eV, E_F = 0
        out = []
        for v in VALLEYS:
            t = HB_EVNM2 / v["mz"] / self.dz**2
            diag = Ec[1:-1] + 2 * t
            off = -t * np.ones(len(diag) - 1)
            w, U = eigh_tridiagonal(diag, off, select="i", select_range=(0, nlev - 1))
            chi2 = np.zeros((nlev, len(self.z)))
            chi2[:, 1:-1] = U.T**2 / self.dz  # normalised: sum chi2 * dz = 1  (1/nm)
            out.append((w, chi2, v))
        return out

    def electron_density(self, psi):
        """n(z) in m^-3 and areal ns in cm^-2, plus subband list."""
        kT_eV = KB * self.T / Q
        n = np.zeros_like(self.z)
        bands = self._subbands(psi)
        ns = 0.0
        for w, chi2, v in bands:
            for E, c2 in zip(w, chi2):
                occ_m2 = v["g"] * (v["md"] * M0 * KB * self.T / (np.pi * HBAR**2)) * np.log1p(np.exp(-E / kT_eV))  # 1/m^2
                n += occ_m2 * c2 * 1e9  # chi2 in 1/nm -> 1/m
                ns += occ_m2
        return n, ns * 1e-4, bands

    # --- Poisson ---
    def _poisson(self, psi_s, psi0, n_prev):
        z = self.z * 1e-9
        h = self.dz * 1e-9
        N = len(z)
        psi = psi0.copy()
        for _ in range(60):
            p = self.na * np.exp(-psi / self.vt)
            n = n_prev * np.exp((psi - psi0) / self.vt)
            rho = Q * (p - self.na - n)
            F = np.zeros(N)
            F[1:-1] = (psi[2:] - 2 * psi[1:-1] + psi[:-2]) / h**2 + rho[1:-1] / (EPS_SI * EPS0)
            dr = Q * (-p / self.vt - n / self.vt)
            main = -2 / h**2 + dr[1:-1] / (EPS_SI * EPS0)
            A = diags([np.ones(N - 3) / h**2, main, np.ones(N - 3) / h**2], [-1, 0, 1], format="csc")
            d = spsolve(A, -F[1:-1])
            lam = 1.0
            mx = np.max(np.abs(d))
            if mx > 0.05:
                lam = 0.05 / mx
            psi[1:-1] += lam * d
            psi[0], psi[-1] = psi_s, 0.0
            if mx < 1e-10:
                break
        return psi

    def solve(self, psi_s, iters=200, tol=1e-7):
        z = self.z * 1e-9
        W = np.sqrt(2 * EPS_SI * EPS0 * max(psi_s, 1e-3) / (Q * self.na))
        psi = psi_s * np.clip(1 - z / W, 0, None) ** 2
        n = np.zeros_like(psi)
        for it in range(iters):
            n_new, ns, bands = self.electron_density(psi)
            n = n_new if it == 0 else 0.5 * n + 0.5 * n_new
            psi_new = self._poisson(psi_s, psi, n + 1e-30)
            err = np.max(np.abs(psi_new - psi))
            psi = psi + 0.5 * (psi_new - psi)
            if err < tol:
                break
        n_fin, ns, bands = self.electron_density(psi)
        p = self.na * np.exp(-psi / self.vt)
        h = self.dz * 1e-9
        Es = Q / (EPS_SI * EPS0) * np.trapezoid(self.na + n_fin - p, dx=h)  # surface field (V/m)
        vg = self.vfb + psi_s + EPS_SI * Es * self.tox / EPS_OX
        return dict(psi_s=psi_s, psi=psi, n=n_fin, ns_cm2=ns, Es=Es, Vg=vg, bands=bands, iters=it + 1, err=err)


def vt_from_curve(vg, ns_cm2):
    """Linear extrapolation of ns(Vg) in strong inversion to ns=0, and constant-charge Vg at 1e11 / 1e12 cm^-2."""
    vg, ns = np.asarray(vg), np.asarray(ns_cm2)
    sel = ns > 3e12
    a, b = np.polyfit(vg[sel], ns[sel], 1)
    vt_lin = -b / a
    out = dict(VT_linear_extrap=float(vt_lin), cinv_uF_cm2=float(a * Q * 1e6))
    for target in (1e11, 1e12):
        out[f"Vg_at_ns_{target:.0e}"] = float(np.interp(np.log(target), np.log(ns), vg))
    return out
