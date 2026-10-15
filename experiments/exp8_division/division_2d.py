"""2D division experiment: one large domain (phi=+1 disk in a phi=-1 sea) evolves under Cahn-Hilliard + gate-screened Coulomb frustration.
Scan A (Coulomb strength) and gate distance d; count droplets, droplet radius and the structure-factor wavenumber."""
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from qcell.phasefield import OKField  # noqa: E402

OUT = ROOT / "results" / "exp8"
OUT.mkdir(parents=True, exist_ok=True)
N, L = 128, 100.0
DX = L / N


def droplets(phi):
    """Number of connected positive domains (periodic) and their equivalent radii (grid units -> length)."""
    big = np.tile(phi > 0, (3, 3))
    lab, n = ndimage.label(big)
    cx = ndimage.center_of_mass(big, lab, range(1, n + 1))
    areas = ndimage.sum(big, lab, range(1, n + 1))
    keep = [i for i, c in enumerate(cx) if N <= c[0] < 2 * N and N <= c[1] < 2 * N]
    r = [float(np.sqrt(areas[i] * DX * DX / np.pi)) for i in keep]
    return len(keep), r


def make_initial(R0):
    x = (np.arange(N) + 0.5) * DX
    X, Y = np.meshgrid(x, x, indexing="ij")
    r = np.hypot(X - L / 2, Y - L / 2)
    return np.tanh((R0 - r) / 1.4)


def run_case(A, d, R0=18.0, t_end=4000.0, dt=0.1, seed=0):
    f = OKField(nx=N, ny=N, Lx=L, Ly=L, kappa=1.0, A=A, d=d)
    rng = np.random.default_rng(seed)
    phi = make_initial(R0) + 1e-3 * rng.standard_normal((N, N))
    phi = np.clip(phi, -1.2, 1.2)
    mass0 = phi.sum()
    hist = []
    snap = {}
    steps = int(round(t_end / dt))
    marks = set(int(round(x / dt)) for x in (0, 100, 300, 1000, 2000, t_end))
    for i in range(steps + 1):
        if i % int(round(100 / dt)) == 0:
            n, r = droplets(phi)
            hist.append((i * dt, n, float(np.mean(r)) if r else 0.0))
        if i in marks:
            snap[i * dt] = phi.copy()
        if i < steps:
            phi = f.step(phi, dt)
    n, r = droplets(phi)
    qpk = f.structure_peak(phi)
    return dict(A=A, d=("inf" if np.isinf(d) else d), n_final=n, r_mean=float(np.mean(r)) if r else 0.0,
                r_min=float(min(r)) if r else 0.0, r_max=float(max(r)) if r else 0.0, q_peak=qpk,
                phi_bar=float(phi.mean()), mass_drift=float(abs(phi.sum() - mass0) / N**2), hist=hist), snap


if __name__ == "__main__":
    results = []
    snaps = {}
    for A in (0.05, 0.1, 0.2):
        for d in (np.inf, 12.0, 6.0, 3.0):
            res, sn = run_case(A, d)
            results.append(res)
            snaps[(A, d)] = sn
            print(f"A={A:<5} d={str(res['d']):>4}: n_final={res['n_final']:>3}  r_mean={res['r_mean']:5.2f} (min {res['r_min']:.2f}, max {res['r_max']:.2f})  "
                  f"q_peak={res['q_peak']:.3f}  phi_bar={res['phi_bar']:.3f}  mass drift={res['mass_drift']:.1e}", flush=True)
    (OUT / "division_scan.json").write_text(json.dumps(results, indent=2))
    np.savez_compressed(OUT / "division_snaps.npz", **{f"A{A}_d{d}_t{t}": v for (A, d), s in snaps.items() for t, v in s.items()})
    print("saved")
