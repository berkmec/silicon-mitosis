"""Review (2026-10-11): reproducible versions of two checks that were first run inline (no artifact before the audit).
(1) float32 vs float64 time loop of ActiveCH (256^2, paper rates, t = 500): max |dc| and droplet radius.
(2) Gray-Scott robustness: two replicating (F, k) points x seeds 2, 3, 4 with 5 % initial noise (128^2, 20000 steps).
Writes results/exp11/checks_review.json."""
import os
import sys

os.environ.setdefault("OMP_NUM_THREADS", "1")
import json  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from qcell.active import ActiveCH, rates_zwicker  # noqa: E402
from qcell.grayscott import GrayScott  # noqa: E402

out = {}
res = {}
for name, dt_ in (("float64", np.float64), ("float32", np.float32)):
    m = ActiveCH(n=256, L=128.0, S=1.0, dtype=dt_, **rates_zwicker(0.19, -0.0011))
    c = m.run(m.disk(8.0, ecc=0.05), 500.0, 0.2)
    res[name] = (np.asarray(c, dtype=np.float64), m.droplets(np.asarray(c, dtype=np.float64))["radii"])
out["float32_vs_float64"] = dict(t=500.0, max_abs_dc=float(np.abs(res["float64"][0] - res["float32"][0]).max()),
                                 radii64=res["float64"][1], radii32=res["float32"][1])
print(out["float32_vs_float64"], flush=True)
gs_out = {}
for F, k in ((0.026, 0.059), (0.034, 0.063)):
    for seed in (2, 3, 4):
        gs = GrayScott(n=128, F=F, k=k)
        U, V = gs.seed(np.random.default_rng(seed), noise=0.05)
        rec = []
        gs.run(U, V, 20000, every=1000, cb=lambda t, U, V: rec.append(gs.spots(V)))
        gs_out[f"F{F}_k{k}_s{seed}"] = dict(n=[r["n"] for r in rec], perc_end=rec[-1]["perc"],
                                            n_max=max(r["n"] for r in rec))
        print(F, k, seed, gs_out[f"F{F}_k{k}_s{seed}"]["n_max"], rec[-1]["perc"], flush=True)
out["gray_scott_robustness_5pct_noise"] = gs_out
(ROOT / "results" / "exp11" / "checks_review.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
print("TAMAM", flush=True)
