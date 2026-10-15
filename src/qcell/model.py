"""Recursive-cell model: size scaling, charging energy, confinement, cascade rule,
and orthodox (sequential tunnelling) single-dot transport.

Approximations (see research/notes/model_assumptions.md): spherical self-capacitance,
cubic-box confinement, constant interaction model, carrier-agnostic gate polarity.
"""
from dataclasses import dataclass

import numpy as np

E = 1.602176634e-19
KB = 1.380649e-23
EPS0 = 8.8541878128e-12
HBAR = 1.054571817e-34
M0 = 9.1093837015e-31
SI_LATTICE = 0.543e-9


@dataclass(frozen=True)
class CellParams:
    L0: float = 100e-9        # root cell edge (m)
    r: float = 0.5            # size ratio child/parent
    eps_r: float = 11.7       # silicon surroundings
    eps_ox: float = 3.9       # gate dielectric
    t_ox_frac: float = 0.1    # dielectric thickness / cell edge
    boost: float = 1.0        # "supercapacitor" multiplier on gate capacitance
    m_eff: float = 0.26       # electron effective mass / m0
    e_break: float = 1.0e9    # dielectric breakdown field (V/m), ~10 MV/cm
    margin: float = 10.0      # E >= margin * kT counts as "regime reached"


def size(p: CellParams, g):
    return p.L0 * p.r ** np.asarray(g, dtype=float)


def c_gate(p: CellParams, g):
    L = size(p, g)
    return p.boost * EPS0 * p.eps_ox * L / p.t_ox_frac


def c_total(p: CellParams, g):
    """C_sigma = C_source + C_drain + C_gate, with C_s = C_d = C_self/2."""
    L = size(p, g)
    c_self = 2.0 * np.pi * EPS0 * p.eps_r * L
    return c_self + c_gate(p, g)


def charging_energy(p: CellParams, g):
    return E**2 / (2.0 * c_total(p, g))


def level_spacing(p: CellParams, g):
    """First excitation of a cubic box: 3 * hbar^2 pi^2 / (2 m L^2)."""
    L = size(p, g)
    return 3.0 * HBAR**2 * np.pi**2 / (2.0 * p.m_eff * M0 * L**2)


def regime_generation(p: CellParams, T: float, g_max: int = 60):
    """First generation where E_C and level spacing both exceed margin*kT.

    Returns (g_coulomb, g_quantum, g_both); None if never reached before the
    cell edge falls below one silicon lattice constant.
    """
    out = {"coulomb": None, "confinement": None, "both": None}
    for g in range(g_max):
        if size(p, g) < SI_LATTICE:
            break
        cb = charging_energy(p, g) >= p.margin * KB * T
        qc = level_spacing(p, g) >= p.margin * KB * T
        if cb and out["coulomb"] is None:
            out["coulomb"] = g
        if qc and out["confinement"] is None:
            out["confinement"] = g
        if cb and qc and out["both"] is None:
            out["both"] = g
    return out["coulomb"], out["confinement"], out["both"]


def spawn_voltage(p: CellParams, g, carriers: float = 1.0):
    """Minimum pulse amplitude for generation g to load `carriers` into its child.

    Rule: stored energy 1/2 C_g V^2 >= carriers * E_C(child)  =>  V = e*sqrt(carriers/(C_g*C_{g+1})).
    """
    return E * np.sqrt(carriers / (c_total(p, g) * c_total(p, g + 1)))


def breakdown_voltage(p: CellParams, g):
    return p.e_break * p.t_ox_frac * size(p, g)


def cascade_depth(p: CellParams, v_pulse: float, g_max: int = 60):
    """Deepest generation reached by a fixed-amplitude pulse, and the stop reason."""
    g = 0
    while g < g_max:
        if size(p, g + 1) < SI_LATTICE:
            return g, "lattice"
        if v_pulse > breakdown_voltage(p, g):
            return g, "breakdown"
        if v_pulse < spawn_voltage(p, g):
            return g, "energy"
        g += 1
    return g, "g_max"


def _rate(dE, T, R):
    """Orthodox tunnelling rate for free-energy change dE (J) at temperature T (K)."""
    kT = KB * max(T, 1e-3)
    x = np.clip(dE / kT, -700.0, 700.0)
    with np.errstate(over="ignore", invalid="ignore"):
        f = np.where(np.abs(x) < 1e-8, 1.0, x / np.expm1(x))
    return kT / (E**2 * R) * f


def single_dot_current(C_s, C_d, C_g, Vg, Vsd, T, R_s=1e6, R_d=1e6, n_states=24):
    """Stationary sequential-tunnelling current I(Vg, Vsd) for one dot (A).

    Vg: 1-D array, Vsd: 1-D array -> returns shape (len(Vg), len(Vsd)).
    Source grounded, drain at Vsd. Carrier-agnostic gate polarity.
    """
    Vg = np.asarray(Vg, float)[:, None, None]
    Vd = np.asarray(Vsd, float)[None, :, None]
    Cs = C_s + C_d + C_g
    Ec = E**2 / (2 * Cs)
    ng = (C_g * Vg + C_d * Vd) / E
    n0 = int(np.floor(np.min(ng))) - n_states // 2
    N = (n0 + np.arange(n_states))[None, None, :]
    U = lambda n: Ec * (n - ng) ** 2
    up_s = _rate(U(N + 1) - U(N), T, R_s)
    up_d = _rate(U(N + 1) - U(N) - E * Vd, T, R_d)
    dn_s = _rate(U(N - 1) - U(N), T, R_s)
    dn_d = _rate(U(N - 1) - U(N) + E * Vd, T, R_d)
    up = up_s + up_d
    dn = dn_s + dn_d
    tiny = 1e-300
    logr = np.log(up[..., :-1] + tiny) - np.log(dn[..., 1:] + tiny)
    logP = np.concatenate([np.zeros(up.shape[:2] + (1,)), np.cumsum(logr, axis=-1)], axis=-1)
    logP -= logP.max(axis=-1, keepdims=True)
    P = np.exp(logP)
    P /= P.sum(axis=-1, keepdims=True)
    return E * (P * (dn_d - up_d)).sum(axis=-1)
