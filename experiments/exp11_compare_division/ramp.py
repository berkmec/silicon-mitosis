"""Exp11 exploratory run: does ramping the drive make cells divide into progressively SMALLER cells?
At each stage all rates (nu+, nu-, k) are doubled, so l = sqrt(D/k) shrinks by sqrt(2) while eps = nu+/k and
dnu/k stay fixed. usage: python ramp.py eps dnu_over_k [stages] [T_stage]"""
import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
import json  # noqa: E402
import sys  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from qcell.active import ActiveCH, rates_zwicker  # noqa: E402

OUT = ROOT / "results" / "exp11"
eps, dnu_k = float(sys.argv[1]), float(sys.argv[2])
stages = int(sys.argv[3]) if len(sys.argv) > 3 else 4
T_stage = float(sys.argv[4]) if len(sys.argv) > 4 else 8000.0
k0, n, L, dt = 1e-2, 256, 128.0, 0.2
c, hist, snaps = None, [], {}
for s in range(stages):
    k = k0 * 2**s
    m = ActiveCH(n=n, L=L, S=1.0, **rates_zwicker(eps, dnu_k * k, k))
    if c is None:
        c = m.disk(8.0, ecc=0.05) + 1e-3 * np.random.default_rng(1).standard_normal((n, n))
    t_off = s * T_stage

    def cb(t, c, m=m, k=k, t_off=t_off):
        d = m.droplets(c)
        hist.append(dict(t=t_off + t, stage=s, k=k, l=float(np.sqrt(1 / k)), n=d["n"], perc=d["perc"],
                         rmean=float(np.mean(d["radii"])) if d["radii"] else 0.0,
                         rmax=max(d["radii"]) if d["radii"] else 0.0))

    c = m.run(c, T_stage, dt, every=int(round(250 / dt)), cb=cb)
    snaps[f"stage{s}"] = c.astype(np.float16)
    h = hist[-1]
    print(f"evre {s}: k={k:.3g}, l={np.sqrt(1/k):.2f}w -> N={h['n']}, ort. R={h['rmean']:.2f}w, şerit={h['perc']}",
          flush=True)
tag = f"ramp_e{eps:.2f}_dk{dnu_k:+.2f}"
np.savez_compressed(OUT / f"{tag}.npz", **snaps)
(OUT / f"{tag}.json").write_text(json.dumps(dict(eps=eps, dnu_over_k=dnu_k, T_stage=T_stage, hist=hist)),
                                 encoding="utf-8")
