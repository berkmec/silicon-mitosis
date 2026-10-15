"""Sequential-tunnelling (orthodox) master equation for a chain of islands.

Islands i=0..N-1 are connected in series: S - 0 - 1 - ... - N-1 - D.
Capacitances come from a Maxwell matrix (FEM). Carrier-agnostic polarity (+e carriers).
"""
import itertools

import numpy as np

from .model import E, KB


def _rate(dE, T, R):
    kT = KB * max(T, 1e-3)
    x = np.clip(dE / kT, -700.0, 700.0)
    f = 1.0 if abs(x) < 1e-8 else x / np.expm1(x)
    return kT / (E**2 * R) * f


class Chain:
    """C_II: island block (Maxwell, positive diagonal). c_lead: dict lead -> positive coupling vector (-C_Ik)."""

    def __init__(self, C_II, c_S, c_D, c_G, R=1e6, nmax=3):
        self.Cinv = np.linalg.inv(np.asarray(C_II, float))
        self.cS, self.cD, self.cG = (np.asarray(v, float) for v in (c_S, c_D, c_G))
        self.N = len(self.cS)
        self.R = R
        self.states = np.array(list(itertools.product(range(nmax + 1), repeat=self.N)))
        self.index = {tuple(s): k for k, s in enumerate(self.states)}

    def U(self, n, q0):
        d = n * E - q0
        return 0.5 * d @ self.Cinv @ d

    def stationary(self, Vg, Vsd, T):
        """Returns mean occupations, source current (A) for gate Vg, drain bias Vsd (source grounded)."""
        q0 = self.cG * Vg + self.cD * Vsd
        ns = len(self.states)
        W = np.zeros((ns, ns))
        i_s_in = np.zeros(ns)  # net electron rate from the source into the chain, per state
        Uc = np.array([self.U(s, q0) for s in self.states])
        N = self.N
        for k, s in enumerate(self.states):
            for i in range(N):
                for sign in (+1, -1):
                    s2 = s.copy()
                    s2[i] += sign
                    k2 = self.index.get(tuple(s2))
                    if k2 is None:
                        continue
                    # lead hops only at the chain ends
                    if i == 0:
                        dE = Uc[k2] - Uc[k] - (E * 0.0 if sign > 0 else -E * 0.0)
                        r = _rate(dE, T, self.R)  # source at 0 V
                        W[k2, k] += r
                        W[k, k] -= r
                        i_s_in[k] += sign * r
                    if i == N - 1:
                        dE = Uc[k2] - Uc[k] - sign * E * Vsd
                        r = _rate(dE, T, self.R)
                        W[k2, k] += r
                        W[k, k] -= r
            for i in range(N - 1):
                for sign in (+1, -1):  # hop i -> i+1 (sign=+1) or i+1 -> i
                    s2 = s.copy()
                    s2[i] -= sign
                    s2[i + 1] += sign
                    if s2.min() < 0:
                        continue
                    k2 = self.index.get(tuple(s2))
                    if k2 is None:
                        continue
                    r = _rate(Uc[k2] - Uc[k], T, self.R)
                    W[k2, k] += r
                    W[k, k] -= r
        A = W.copy()
        A[-1, :] = 1.0
        b = np.zeros(ns)
        b[-1] = 1.0
        p = np.linalg.solve(A, b)
        p = np.clip(p, 0, None)
        p /= p.sum()
        return (p @ self.states), E * float(p @ i_s_in)
