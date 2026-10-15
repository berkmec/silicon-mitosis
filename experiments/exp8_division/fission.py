"""Coulomb fission of charged Lennard-Jones clusters with gate screening.
For each (d, g): E_min(n) for n=1..NMAX by basin hopping (warm-started along g), then compare E_single(N) with the best split
E_min(N1)+E_min(N-N1) (fragments infinitely far apart, valid for screened and for unscreened repulsion)."""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from qcell import cluster  # noqa: E402

OUT = ROOT / "results" / "exp8"
OUT.mkdir(parents=True, exist_ok=True)
NMAX = 26
NS = (10, 13, 16, 19, 22, 26)
GS = np.geomspace(0.02, 1.2, 12)
DS = (np.inf, 4.0, 2.0)
TAG = ""
if len(sys.argv) > 1:  # single gate distance, separate output file (parallel runs)
    DS = (float(sys.argv[1]),)
    TAG = "_d" + sys.argv[1]

# ---- validation at g = 0 -----------------------------------------------------------
ref = {7: -16.505384, 13: -44.326801, 19: -72.659782}
val = {}
for n, e_ref in ref.items():
    e, _ = cluster.minimise(n, niter=60, seed=2)
    val[n] = (e, e_ref)
    print(f"LJ{n}: E = {e:.6f}  literature {e_ref:.6f}  diff {abs(e-e_ref):.2e}", flush=True)

res = {"validation": {str(k): dict(E=v[0], E_ref=v[1]) for k, v in val.items()}, "scan": []}
base = {n: cluster.minimise(n, niter=25, seed=n)[1] for n in range(2, NMAX + 1)}
for d in DS:
    prev = dict(base)
    for g in GS:
        E = {1: 0.0}
        cur = {}
        for n in range(2, NMAX + 1):
            e, x = cluster.minimise(n, g=g, d=d, niter=25, seed=int(1000 * g) + n, x0=prev[n])
            E[n], cur[n] = e, x
        prev = cur
        for N in NS:
            splits = [(E[a] + E[N - a], a) for a in range(1, N // 2 + 1)]
            e_split, a_best = min(splits)
            res["scan"].append(dict(d=("inf" if np.isinf(d) else d), g=float(g), N=N, E_single=E[N], E_split=e_split, N1=a_best,
                                    gain=float(E[N] - e_split), fissions=bool(e_split < E[N] - 1e-6)))
        print(f"d={d} g={g:.3f}: " + " ".join(f"N={N}:{'S' if E[N] > min(E[a]+E[N-a] for a in range(1,N//2+1)) + 1e-6 else '-'}" for N in NS), flush=True)
        (OUT / f"fission_scan{TAG}.json").write_text(json.dumps(res, indent=2))

# ---- thresholds g_c(N, d) and scaling --------------------------------------------------
summary = {}
for d in DS:
    key = "inf" if np.isinf(d) else d
    for N in NS:
        rows = sorted([r for r in res["scan"] if r["d"] == key and r["N"] == N], key=lambda r: r["g"])
        gc = None
        for r in rows:
            if r["fissions"]:
                gc = r["g"]
                break
        summary[f"d={key},N={N}"] = gc
res["g_c"] = summary
for d in DS:
    key = "inf" if np.isinf(d) else d
    xs = [N for N in NS if summary[f"d={key},N={N}"] is not None]
    ys = [summary[f"d={key},N={N}"] for N in xs]
    if len(xs) >= 3:
        slope = np.polyfit(np.log(xs), np.log(ys), 1)[0]
        res[f"slope_gc_vs_N_d={key}"] = float(slope)
        print(f"d={key}: g_c(N) = {dict(zip(xs, [round(v,3) for v in ys]))}   log-log slope = {slope:.2f}  (Rayleigh/liquid-drop: -1)")
(OUT / f"fission_scan{TAG}.json").write_text(json.dumps(res, indent=2))
print("done")
