"""Pilot: (B) Zwicker active droplet in 2D with the paper's rates, single droplet; (GS) Gray-Scott at Pearson's F=0.02, k=0.059."""
import json
import sys
from pathlib import Path

import numpy as np
from scipy import ndimage

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from qcell.active import ActiveCH  # noqa: E402

OUT = ROOT / "results" / "exp11"
OUT.mkdir(parents=True, exist_ok=True)
which = sys.argv[1]

if which == "B":
    m = ActiveCH(n=128, L=64.0, nu_p=2e-3, nu_m=7e-3, k_p=1e-2, k_m=1e-2)
    c = m.disk(6.0) + 1e-3 * np.random.default_rng(1).standard_normal((128, 128))
    hist = []
    dt, T = 0.2, 30000.0
    for i in range(int(T / dt) + 1):
        if i % int(500 / dt) == 0:
            n, r = m.droplets(c)
            hist.append((i * dt, n, float(np.mean(r)) if r else 0.0))
            print(f"t={i*dt:7.0f}  damla={n:3d}  ort. yarıçap={hist[-1][2]:5.2f}", flush=True)
        if i in (int(5000 / dt), int(15000 / dt), int(30000 / dt)):
            np.save(OUT / f"pilotB_t{int(i*dt)}.npy", c.astype(np.float32))
        c = m.step(c, dt)
    (OUT / "pilotB.json").write_text(json.dumps(hist))

if which == "GS":
    # Pearson 1993: Du=2e-5, Dv=1e-5, 2.5 x 2.5, 256^2, dt=1, periodic, 20x20 square perturbed to (1/2, 1/4) with +-1% noise
    N, Lb = 256, 2.5
    h = Lb / N
    Du, Dv, F, k = 2e-5, 1e-5, 0.02, 0.059
    rng = np.random.default_rng(0)
    U = np.ones((N, N)); V = np.zeros((N, N))
    s = slice(N // 2 - 10, N // 2 + 10)
    U[s, s], V[s, s] = 0.5, 0.25
    U *= 1 + 0.01 * (2 * rng.random((N, N)) - 1); V *= 1 + 0.01 * (2 * rng.random((N, N)) - 1)

    def lap(Z):
        return (np.roll(Z, 1, 0) + np.roll(Z, -1, 0) + np.roll(Z, 1, 1) + np.roll(Z, -1, 1) - 4 * Z) / h**2

    def spots(V):
        big = np.tile(V > 0.1, (3, 3))
        lab, kk = ndimage.label(big)
        if kk == 0:
            return 0
        com = ndimage.center_of_mass(big, lab, range(1, kk + 1))
        return sum(1 for cc in com if N <= cc[0] < 2 * N and N <= cc[1] < 2 * N)

    hist = []
    for t in range(40001):
        if t % 1000 == 0:
            hist.append((t, spots(V)))
            print(f"GS t={t:6d} leke={hist[-1][1]}", flush=True)
        if t in (4000, 12000, 40000):
            np.save(OUT / f"pilotGS_t{t}.npy", V.astype(np.float32))
        uvv = U * V * V
        U, V = U + (Du * lap(U) - uvv + F * (1 - U)), V + (Dv * lap(V) + uvv - (F + k) * V)
    (OUT / "pilotGS.json").write_text(json.dumps(hist))
