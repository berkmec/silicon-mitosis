"""Negative-result record: pre-formed single domains (disk, rod, ellipse) under Coulomb frustration do not split in this 2D model.
Records droplet count and mean equivalent radius up to t=3000 (d = inf)."""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(Path(__file__).parent))
import division_2d as D  # noqa: E402
from qcell.phasefield import OKField  # noqa: E402

x = (np.arange(D.N) + 0.5) * D.DX
X, Y = np.meshgrid(x, x, indexing="ij")


def rod(Lr=60.0, w=12.0):
    dx = np.clip(np.abs(X - D.L / 2) - (Lr / 2 - w), 0, None)
    return np.tanh((w - np.hypot(dx, Y - D.L / 2)) / 1.4)


def ell(a=30.0, b=9.0):
    r = np.sqrt(((X - D.L / 2) / a) ** 2 + ((Y - D.L / 2) / b) ** 2)
    return np.tanh((1 - r) * b / 1.4)


cases = {"disk R0=18": (D.make_initial(18.0), (0.0, 0.004, 0.016, 0.05)),
         "rod 60x12": (rod(), (0.0, 0.004, 0.01, 0.02)),
         "ellipse 30x9": (ell(), (0.0, 0.004, 0.01, 0.02))}
rng = np.random.default_rng(0)
out = []
for name, (ic, As) in cases.items():
    for A in As:
        f = OKField(nx=D.N, ny=D.N, Lx=D.L, Ly=D.L, kappa=1.0, A=A, d=np.inf)
        phi = ic + 1e-3 * rng.standard_normal((D.N, D.N))
        t, ser = 0.0, []
        for t_next in (100, 300, 1000, 3000):
            phi, _ = f.run(phi, t_end=float(t_next - t), dt=0.25)
            t = t_next
            n, r = D.droplets(phi)
            ser.append(dict(t=t, n=n, r_mean=float(np.mean(r)) if r else 0.0, pos_frac=float((phi > 0).mean())))
        out.append(dict(initial=name, A=A, series=ser, phi_bar=float(phi.mean())))
        s = ser[-1]
        print(f"{name:<13} A={A:<6}: t=3000 n={s['n']} r={s['r_mean']:.1f} pos={s['pos_frac']:.3f}", flush=True)
(ROOT / "results" / "exp8" / "division_diag.json").write_text(json.dumps(out, indent=2))
print("done")
