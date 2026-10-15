"""Exp14: Cartesian d-dimensional active Cahn-Hilliard on a periodic box (up to the 6D hypercube), lean float32 version of
qcell.active.ActiveCH (same equation, same semi-implicit stabilised Fourier step, same s(c)) that never builds meshgrids,
so a 16^6 box needs < 1 GB. Droplets are counted with periodic union-find labelling.
usage: python cart_run.py <name> <d> <n_long> <n_other> <dx> <eps> <nu_m/k> <k> <R0> <T> <dt> [ecc] [seed]"""
import os
import sys

os.environ.setdefault("OMP_NUM_THREADS", "1")
if sys.platform == "win32":
    import ctypes
    from ctypes import wintypes

    _k = ctypes.windll.kernel32
    _k.GetCurrentProcess.restype = wintypes.HANDLE
    _k.SetPriorityClass.argtypes = [wintypes.HANDLE, wintypes.DWORD]
    _k.SetPriorityClass(_k.GetCurrentProcess(), 0x00004000)
import json  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
from scipy import fft as sfft  # noqa: E402
from scipy import ndimage, special  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from qcell.active import ActiveCH  # noqa: E402

OUT = ROOT / "results" / "exp14_dimension" / "cart"
F32 = np.float32


def _b(v, i, d):
    return v.reshape([-1 if j == i else 1 for j in range(d)])


class CartCH:
    def __init__(self, shape, dx, dt, S=1.0, nu_p=0.0, nu_m=0.0, k_p=0.0, k_m=0.0):
        self.shape, self.d, self.dx, self.dt, self.S = tuple(shape), len(shape), dx, dt, S
        ks = [2 * np.pi * np.fft.fftfreq(n, d=dx) for n in shape[:-1]] + [2 * np.pi * np.fft.rfftfreq(shape[-1], d=dx)]
        q2 = np.zeros([len(k) for k in ks], F32)
        for i, k in enumerate(ks):
            q2 += _b((k**2).astype(F32), i, self.d)
        self.q2 = q2
        self.den = (1.0 + dt * (S * q2 + 0.25 * q2 * q2)).astype(F32)
        self.chem = ActiveCH(n=4, L=4.0, dim=1, nu_p=nu_p, nu_m=nu_m, k_p=k_p, k_m=k_m)
        self.t = 0.0

    def step(self, c, n=1):
        h, S = F32(self.dt), F32(self.S)
        for _ in range(n):
            nl = c * (c - F32(1)) * (F32(2) * c - F32(1))
            ch = sfft.rfftn(c)
            rhs = ch - h * self.q2 * (sfft.rfftn(nl) - S * ch)
            del nl
            rhs += h * sfft.rfftn(self.chem.s(c))
            rhs /= self.den
            c = sfft.irfftn(rhs, s=self.shape).astype(F32, copy=False)
            self.t += self.dt
        return c

    def droplet(self, R, ecc=0.05):
        """Prolate droplet (long axis = axis 0) at the box centre on top of the far field c = eps."""
        r2 = np.zeros(self.shape, F32)
        for i, n in enumerate(self.shape):
            x = (np.arange(n) + 0.5) * self.dx - n * self.dx / 2
            a = R * (1 + ecc) if i == 0 else R * (1 - ecc)
            r2 += _b(((x / a) ** 2).astype(F32), i, self.d)
        cf = self.chem.nu_p / self.chem.k_p if self.chem.k_p > 0 else 0.0
        rho = np.sqrt(r2, out=r2)
        return (cf + (1.0 - cf) * 0.5 * (1.0 - np.tanh((rho - 1.0) * R))).astype(F32)


def components(c, dx, thr=0.5):
    """Periodic connected components of c > thr: list of dict(vol, R_eq, extent (per axis, in grid points), span)."""
    b = c > thr
    d = b.ndim
    lab, n = ndimage.label(b)
    if n == 0:
        return []
    parent = np.arange(n + 1)

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for ax in range(d):
        a0, a1 = np.take(lab, 0, axis=ax), np.take(lab, -1, axis=ax)
        m = (a0 > 0) & (a1 > 0)
        if m.any():
            for p, q in np.unique(np.stack([a0[m], a1[m]], 1), axis=0):
                rp, rq = find(p), find(q)
                if rp != rq:
                    parent[rq] = rp
    root = np.array([find(i) for i in range(n + 1)])
    lab = root[lab]
    ids = np.unique(lab[lab > 0])
    vols = np.bincount(lab.ravel(), minlength=n + 1)
    out = []
    unit = np.pi ** (d / 2) / special.gamma(d / 2 + 1)
    for i in ids:
        msk = lab == i
        ext = [int(np.any(msk, axis=tuple(j for j in range(d) if j != ax)).sum()) for ax in range(d)]
        vol = float(vols[i]) * dx**d
        out.append(dict(vol=vol, R_eq=float((vol / unit) ** (1 / d)), extent=ext,
                        span=[e == s for e, s in zip(ext, b.shape)]))
    return sorted(out, key=lambda o: -o["vol"])


def main():
    a = sys.argv[1:]
    name, d, n_long, n_other, dx = a[0], int(a[1]), int(a[2]), int(a[3]), float(a[4])
    eps, nmk, k, R0, T, dt = map(float, a[5:11])
    ecc = float(a[11]) if len(a) > 11 else 0.05
    seed = int(a[12]) if len(a) > 12 else 1
    shape = (n_long,) + (n_other,) * (d - 1)
    m = CartCH(shape, dx, dt, nu_p=eps * k, nu_m=nmk * k, k_p=k, k_m=k)
    c = m.droplet(R0, ecc)
    c += (1e-3 * np.random.default_rng(seed).standard_normal(shape)).astype(F32)
    hist, snaps, t0 = [], {}, time.time()
    every = max(1, int(round(25.0 / dt)))
    OUT.mkdir(parents=True, exist_ok=True)
    mid = tuple(s // 2 for s in shape[2:])
    while m.t < T - 1e-9:
        c = m.step(c, every)
        comp = components(c, dx)
        big = [o for o in comp if o["R_eq"] >= 1.5]
        hist.append(dict(t=round(m.t, 2), n=len(big), comps=comp[:12], cmin=float(c.min()), cmax=float(c.max()),
                         cmean=float(c.mean()), wall=round(time.time() - t0, 1)))
        tt = int(round(m.t))
        if tt % 100 == 0:
            snaps[f"slice_t{tt}"] = c[(slice(None), slice(None)) + mid].astype(np.float16)
            snaps[f"proj_t{tt}"] = c.reshape(shape[0], shape[1], -1).max(axis=2).astype(np.float16)
            print(f"t={m.t:7.1f} N={len(big)} R_eq={[round(o['R_eq'], 2) for o in big[:6]]} "
                  f"ext0={[o['extent'][0] for o in big[:6]]} cmin={c.min():.2f} ({time.time() - t0:.0f} s)", flush=True)
            (OUT / f"{name}.json").write_text(json.dumps(dict(name=name, d=d, shape=shape, dx=dx, eps=eps, nmk=nmk, k=k,
                                                              R0=R0, dt=dt, ecc=ecc, seed=seed, hist=hist)),
                                              encoding="utf-8")
            np.savez_compressed(OUT / f"{name}.npz", **snaps)
        if not comp:
            break
    print("TAMAM", flush=True)


if __name__ == "__main__":
    main()
