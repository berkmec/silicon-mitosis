"""Exp14: why the division window closes at high d (nonlinear sharp-interface theory, Zwicker chemistry).
For each d: R_hollow (largest radius with a dense core; above it the droplet must become a liquid shell) and, at several
radii, the growth velocity v, the l = 2 rate sigma_2 and the relative elongation rate sigma_2 - v / R of a growing droplet.
Writes results/exp14_dimension/window_theory.json."""
import json
import os
import sys
import warnings
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np  # noqa: E402

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from qcell import active_sharp_nl as nl  # noqa: E402

EPS, K = 0.20, 0.01


def r_hollow(d, ch, lo=6.0, hi=60.0):
    if nl.base(lo, d, ch) is None:
        return None
    if nl.base(hi, d, ch) is not None:
        return f">{hi}"
    while hi - lo > 0.05:
        mid = 0.5 * (lo + hi)
        if nl.base(mid, d, ch) is None:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


out = {}
for dnu in (-0.001, 0.004):
    ch = nl.Chem(EPS, dnu + K - EPS * K, K)
    for d in (2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0):
        rh = r_hollow(d, ch)
        rows = []
        for R in (6.0, 8.0, 10.0, 12.0):
            b = nl.base(R, d, ch)
            if b is None:
                rows.append(dict(R=R, hollow=True))
                continue
            v = (b["dmo"] - b["dmi"]) / b["dc"]
            s2 = nl.sigma(R, d, 2, 0.0, ch, b=b)
            rows.append(dict(R=R, v=v, sigma2=s2, rel=s2 - v / R, c0=b["c0"]))
        out[f"{dnu}_{d}"] = dict(dnu=dnu, d=d, R_hollow=rh, rows=rows)
        print(f"dnu={dnu} d={d}: R_hollow={rh} | " + " ".join(
            f"R={r['R']:.0f}:" + ("hollow" if r.get("hollow") else f"v={r['v']:.1e},s2={r['sigma2']:.1e},rel={r['rel']:+.1e}")
            for r in rows), flush=True)
        (ROOT / "results" / "exp14_dimension" / "window_theory.json").write_text(json.dumps(out, indent=1, default=str),
                                                                               encoding="utf-8")
