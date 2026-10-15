"""Single-phonon (longitudinal acoustic, deformation potential) golden-rule rates between two FEM states.

Gamma_em(Delta) = Xi^2 q0^3 <|F(q0)|^2>_angles (n+1) / (2 pi hbar rho v^2),  q0 = Delta/(hbar v),
F(q) = <psi2| exp(i q.r) |psi1>.  Gamma1 = Gamma_em + Gamma_abs = Gamma0 (2n+1).
Positions in nm inside the FEM; converted to metres here.
"""
import numpy as np
from skfem import Basis, ElementTetP2

from .model import E, HBAR, KB

# Silicon (bulk) parameters
XI = 9.0 * E          # effective LA deformation potential (J); literature range ~ 5-10 eV (assumption)
RHO = 2329.0          # kg/m^3
V_LA = 8430.0         # m/s


def fib_sphere(n):
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    th = np.pi * (1 + 5**0.5) * i
    return np.stack([np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)], axis=1)


def refined_tet_quadrature(order=6, levels=1):
    """Quadrature on the reference tet obtained by midpoint-subdividing it `levels` times (8^levels children)."""
    from skfem.quadrature import get_quadrature_tet
    X, W = get_quadrature_tet(order)  # X: (3, n), W sums to 1/6
    tets = [np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], float)]
    for _ in range(levels):
        new = []
        for v in tets:
            v0, v1, v2, v3 = v
            m01, m02, m03, m12, m13, m23 = (v0 + v1) / 2, (v0 + v2) / 2, (v0 + v3) / 2, (v1 + v2) / 2, (v1 + v3) / 2, (v2 + v3) / 2
            new += [np.array(t) for t in (
                (v0, m01, m02, m03), (v1, m01, m12, m13), (v2, m02, m12, m23), (v3, m03, m13, m23),
                (m02, m13, m01, m03), (m02, m13, m03, m23), (m02, m13, m23, m12), (m02, m13, m12, m01))]
        tets = new
    Xs, Ws = [], []
    for v in tets:
        J = (v[1:] - v[0]).T  # (3,3)
        Xs.append(J @ X + v[0][:, None])
        Ws.append(W * abs(np.linalg.det(J)))
    return np.concatenate(Xs, axis=1), np.concatenate(Ws)


class Overlap:
    """Quadrature data for F(q) of two FEM states (positions in nm)."""

    def __init__(self, basis, psi1, psi2, order=6, levels=1):
        # re-integrate on the same P2 mesh with a refined rule: exp(i q.x) oscillates inside elements
        basis = Basis(basis.mesh, ElementTetP2(), quadrature=refined_tet_quadrature(order, levels))
        x = basis.global_coordinates().value  # (3, nel, nqp), nm
        self.x = x.reshape(3, -1) * 1e-9
        w = basis.dx.reshape(-1)
        # psi are M-normalised in nm^-3/2: product is nm^-3, weights dx are nm^3 -> dimensionless
        self.f = (basis.interpolate(psi1).value * basis.interpolate(psi2).value).reshape(-1) * w

    def F(self, qvec):
        return np.sum(self.f * np.exp(1j * (qvec @ self.x)))

    def mean_F2(self, q0, ndir=48):
        dirs = fib_sphere(ndir)
        return float(np.mean([abs(self.F(q0 * d)) ** 2 for d in dirs]))

    def dipole(self, axis=0):
        """x12 = <psi2|x_axis|psi1> in nm (small-q check)."""
        return float(np.sum(self.f * self.x[axis])) * 1e9


def rate_em_zero_T(delta_J, mean_F2):
    q0 = delta_J / (HBAR * V_LA)
    return XI**2 * q0**3 * mean_F2 / (2 * np.pi * HBAR * RHO * V_LA**2)


def coherence(delta_J, mean_F2, T):
    """T1, T2 (upper bound: single-phonon only, T2 = 2 T1) and thermal polarisation."""
    n = 1.0 / np.expm1(delta_J / (KB * T)) if T > 0 else 0.0
    g0 = rate_em_zero_T(delta_J, mean_F2)
    g1 = g0 * (2 * n + 1)
    T1 = 1.0 / g1 if g1 > 0 else np.inf
    return dict(Gamma0=g0, Gamma1=g1, T1=T1, T2_max=2 * T1, n_bose=n, P=float(np.tanh(delta_J / (2 * KB * T))))
