"""Uniform 2D electron gas energetics for the physical meaning of phi: kinetic + exchange (Hartree-Fock, any kernel)
+ correlation (Attaccalite et al., PRL 88, 256601 (2002), zeta=0, single valley), and the local inverse
compressibility K(n) = d^2(E/A)/dn^2.  SI units throughout (n in m^-2, energies per area in J/m^2).

Coulomb kernels (2D Fourier transform of the in-plane interaction):
  unscreened, uniform eps:          V(q) = e^2 / (2 eps0 eps q)
  interface Si | SiO2 slab | gate:  V(q) = e^2 / (eps0 q (eps_Si + eps_ox coth(q t_ox)))   (sheet at the interface)
Optional Fang-Howard form factor F(q) = (1 + 9/8 x + 3/8 x^2)/(1+x)^3, x=q/b, b = 3/z_avg (from memory, sensitivity only).
"""
import numpy as np
from scipy.integrate import quad

from .model import E, EPS0, HBAR, M0

# Attaccalite et al. 2002, Table II, i=0 (zeta=0), Hartree units; D0 = -A0*H0
ATT = dict(A=-0.1925, B=0.0863136, C=0.0572384, E=1.0022, F=-0.02069, G=0.33997, H=1.747e-2)


def alpha0(rs):
    a = ATT
    D = -a["A"] * a["H"]
    den = a["E"] * rs + a["F"] * rs**1.5 + a["G"] * rs**2 + a["H"] * rs**3
    return a["A"] + (a["B"] * rs + a["C"] * rs**2 + D * rs**3) * np.log(1 + 1 / den)


def units(m_eff, eps):
    """Effective Hartree (J) and Bohr radius (m)."""
    aB = 4 * np.pi * EPS0 * eps * HBAR**2 / (m_eff * M0 * E**2)
    Ha = E**2 / (4 * np.pi * EPS0 * eps * aB)
    return Ha, aB


def rs_of(n, m_eff, eps):
    _, aB = units(m_eff, eps)
    return 1.0 / (aB * np.sqrt(np.pi * n))


def V_uniform(eps):
    return lambda q: E**2 / (2 * EPS0 * eps * q)


def V_gate(eps_si, eps_ox, t_ox):
    def V(q):
        x = q * t_ox
        coth = 1.0 / np.tanh(x) if x < 350 else 1.0
        return E**2 / (EPS0 * q * (eps_si + eps_ox * coth))
    return V


def fang_howard(z_avg):
    b = 3.0 / z_avg
    return lambda q: (1 + 9 / 8 * q / b + 3 / 8 * (q / b) ** 2) / (1 + q / b) ** 3


def kinetic_per_area(n, g, m_eff):
    return np.pi * HBAR**2 * n**2 / (g * m_eff * M0)


def exchange_per_area(n, g, V, form=None):
    """Hartree-Fock exchange of g equal species (density n/g each) with interaction kernel V(q) [J m^2]."""
    kF = np.sqrt(4 * np.pi * n / g)

    def integrand(q):
        x = q / (2 * kF)
        overlap = 2 * kF**2 * (np.arccos(x) - x * np.sqrt(max(0.0, 1 - x * x)))
        f = 1.0 if form is None else form(q)
        return q * V(q) * f * overlap

    val, _ = quad(integrand, 0.0, 2 * kF, epsabs=0.0, epsrel=1e-11, limit=400)
    return -g * 0.5 * (2 * np.pi) / (2 * np.pi) ** 4 * val


def correlation_per_area(n, m_eff, eps):
    """Attaccalite zeta=0 correlation per particle times n (single-valley parametrization; for g_v=2 only a proxy)."""
    Ha, _ = units(m_eff, eps)
    return n * alpha0(rs_of(n, m_eff, eps)) * Ha


def energy_per_area(n, g, m_eff, V_x, eps_c, with_corr=True, form=None):
    e = kinetic_per_area(n, g, m_eff) + exchange_per_area(n, g, V_x, form)
    if with_corr:
        e += correlation_per_area(n, m_eff, eps_c)
    return e


def inverse_compressibility(n, h=1e-2, **kw):
    """K(n) = d^2(E/A)/dn^2 by a central difference in log n (J m^2)."""
    ep = energy_per_area(n * (1 + h), **kw)
    e0 = energy_per_area(n, **kw)
    em = energy_per_area(n * (1 - h), **kw)
    # non-uniform 3-point second derivative
    n1, n2 = n * (1 - h), n * (1 + h)
    return 2 * ((ep - e0) / (n2 - n) - (e0 - em) / (n - n1)) / (n2 - n1)
