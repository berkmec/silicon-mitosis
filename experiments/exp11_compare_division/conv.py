"""Convergence of the 2D active-droplet runs: time step, grid spacing and box size (paper rates, eps=0.19, dnu=-1.1e-3)."""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from qcell.active import ActiveCH, rates_zwicker  # noqa: E402

OUT = ROOT / "results" / "exp11"
OUT.mkdir(parents=True, exist_ok=True)
cases = [dict(n=128, L=64.0, dt=0.5), dict(n=128, L=64.0, dt=0.2), dict(n=192, L=64.0, dt=0.5),
         dict(n=256, L=128.0, dt=0.5)]
res = []
for cs in cases:
    m = ActiveCH(n=cs["n"], L=cs["L"], **rates_zwicker(0.19, -1.1e-3))
    c = m.disk(8.0, ecc=0.05)
    rec = []

    def cb(t, c):
        d = m.droplets(c)
        rec.append((t, d["n"], max(d["radii"]) if d["radii"] else 0.0))

    m.run(c, 4000.0, cs["dt"], every=int(round(500 / cs["dt"])), cb=cb)
    res.append(dict(**cs, R_t=rec))
    print(cs, " ".join(f"{t:.0f}:{r:.3f}" for t, _, r in rec), flush=True)
ref = np.array([r for _, _, r in res[0]["R_t"]])
for r in res[1:]:
    print(r["n"], r["L"], r["dt"], "max |ΔR| vs (128,64,0.5):", float(np.max(np.abs(np.array([x for _, _, x in r["R_t"]]) - ref))))
(OUT / "convergence.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
