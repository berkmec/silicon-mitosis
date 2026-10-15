"""Compact nFET model for exp12 (kink / floating-body feedback), per unit channel width.

Forms follow BSIM4 (the model used by the official SKY130 card):
  V_T(V_bs) = V_T0 + K1 (sqrt(Phi_s - V_bs) - sqrt(Phi_s)) - K2 V_bs                      (body effect)
  I_d       = W mu C_ox / L * (V_ov V_de - V_de^2 / 2) * (1 + lam (V_ds - V_de)),             V_ov = V_gs - V_T
  V_dsat    = V_ov E_sat L / (V_ov + E_sat L),  E_sat = 2 v_sat / mu  (velocity saturation)
  V_de      = smooth min(V_ds, V_dsat)  (BSIM4 'delta' smoothing)
  I_ii      = (ALPHA0 + ALPHA1 L) / L * (V_ds - V_de) * exp(-BETA0 / (V_ds - V_de)) * I_d   (BSIM4 Isub form)
Subthreshold conduction is not modelled (V_ov is softplus-smoothed with slope n kT/e).
Parameters for 4 K are NOT measured: they are assumptions documented by the caller.
"""
from dataclasses import dataclass, replace

import numpy as np

EPS0, Q, KB = 8.8541878128e-12, 1.602176634e-19, 1.380649e-23


@dataclass(frozen=True)
class FetParams:
    L: float = 0.15e-6 - 2 * 1.1932e-8      # effective channel length (SKY130 lint), m
    tox: float = 4.148e-9                   # SKY130 toxe, m
    vt0: float = 0.52                       # V (fitted to ngspice at 300 K by the caller)
    k1: float = 0.54086565                  # SKY130 k1, V^0.5
    k2: float = -0.026724591                # SKY130 k2
    phis: float = 0.86                      # surface potential 2 phi_F (exp5: phi_F = 0.430 V)
    mu: float = 0.0318614                   # SKY130 u0, m^2/Vs (effective, fitted by caller)
    vsat: float = 8.0e4                     # SKY130 vsat, m/s
    lam: float = 0.05                       # output conductance (1/V), fitted by caller
    delta: float = 0.01                     # SKY130 delta (V_de smoothing)
    alpha0: float = 3.0e-8                  # SKY130 alpha0, m/V
    alpha1: float = 0.85                    # SKY130 alpha1, 1/V
    beta0: float = 13.86                    # SKY130 beta0, V
    T: float = 300.0                        # K (only used for the softplus width)
    nss: float = 1.5                        # subthreshold slope factor for the smoothing

    @property
    def cox(self):
        return EPS0 * 3.9 / self.tox

    def at(self, **kw):
        return replace(self, **kw)


def vt(p: FetParams, vbs):
    vbs = np.minimum(vbs, p.phis - 1e-3)
    return p.vt0 + p.k1 * (np.sqrt(p.phis - vbs) - np.sqrt(p.phis)) - p.k2 * vbs


def _softplus(x, s):
    return s * np.logaddexp(0.0, x / s)


def vdsat(p: FetParams, vgs, vbs):
    vov = _softplus(vgs - vt(p, vbs), p.nss * KB * p.T / Q)
    esl = 2 * p.vsat / p.mu * p.L
    return vov * esl / (vov + esl), vov


def vde(p: FetParams, vds, vdsat_):
    """BSIM4 effective Vds: smooth min(V_ds, V_dsat)."""
    a = vdsat_ - vds - p.delta
    return vdsat_ - 0.5 * (a + np.sqrt(a * a + 4 * p.delta * vdsat_))


def id_per_w(p: FetParams, vgs, vds, vbs):
    vds_sat, vov = vdsat(p, vgs, vbs)
    v = vde(p, vds, vds_sat)
    esl = 2 * p.vsat / p.mu * p.L
    core = p.mu * p.cox / p.L * (vov * v - 0.5 * v * v) / (1 + v / esl)
    return core * (1 + p.lam * (vds - v))


def iii_per_w(p: FetParams, vgs, vds, vbs):
    vds_sat, _ = vdsat(p, vgs, vbs)
    x = np.maximum(vds - vde(p, vds, vds_sat), 1e-9)
    return (p.alpha0 + p.alpha1 * p.L) / p.L * x * np.exp(-p.beta0 / x) * id_per_w(p, vgs, vds, vbs)
