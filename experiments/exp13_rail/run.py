"""Exp13: subthreshold floating-body latch in a wide SKY130 nFET at 4 K WITH A LOCAL INHIBITOR (resistive drain).
Answers the external review (inceleme.txt): exp12 already had the multi-cell subthreshold model, but its only inhibitor was
the global load line, which by S-NDC theory gives ONE domain. Gray-Scott-like multiplication needs a LOCAL inhibitor:
here the drain voltage drop along a resistive drain (impact ionisation M(V_ds) is very steep, so a conducting filament
suppresses its neighbours). Solver: experiments/exp12_filament/spatial.py (rail_voltage validated: centre drop
r_D I W / 8 exact; ballast decay length sqrt(r_v/r_D) exact on a large domain).
Runs:  A robustness of the exp12 single-filament result over several bistable parameter sets (global coupling only)
       B end-fed resistive drain (contacts at both ends), r_D scan, W = 100 / 300 um
       C distributed ballast (every cell tied to the common node through r_v; inhibition length lambda = sqrt(r_v/r_D)),
         W = 300 um, lambda = 5..40 um, with and without the global load
       D convergence of the R_sq = 1e10 case (dy = 0.1, 0.05; global only)
Lineage per V_DD step (overlap): births, splits, merges, deaths."""
import os
import sys

os.environ.setdefault("OMP_NUM_THREADS", "1")
if sys.platform == "win32":
    import ctypes
    from ctypes import wintypes

    _k32 = ctypes.windll.kernel32
    _k32.GetCurrentProcess.restype = wintypes.HANDLE
    _k32.SetPriorityClass.argtypes = [wintypes.HANDLE, wintypes.DWORD]
    _k32.SetPriorityClass(_k32.GetCurrentProcess(), 0x00004000)
import json  # noqa: E402
from multiprocessing import Pool  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "experiments" / "exp12_filament"))
OUT = ROOT / "results" / "exp13"
OUT.mkdir(parents=True, exist_ok=True)

BASE = dict(ss=0.02, vs=0.005, fii=10, vov0=-0.1, g_sub=1e-13, von=1.05)


def jobs():
    js = []
    # A: other bistable sets with a closing hysteresis window (from exp12 reduced_bistability.json)
    rb = json.loads((ROOT / "results" / "exp12" / "reduced_bistability.json").read_text(encoding="utf-8"))["rows"]
    cand = [r for r in rb if r["hysteresis"] and r["V_down"] is not None and 0.1 < r["width"] < 1.0]
    picked, seen = [], set()
    for r in sorted(cand, key=lambda r: (r["ss"], r["vs"], r["fii"], r["vov0"], r["g_sub"])):
        key = (r["ss"], r["vov0"], r["fii"])
        if key not in seen:
            seen.add(key)
            picked.append(r)
    for i, r in enumerate(picked[:5]):
        p = dict(ss=r["ss"], vs=r["vs"], fii=r["fii"], vov0=r["vov0"], g_sub=r["g_sub"], von=1.05)
        js.append(dict(name=f"A{i}_global", kw=dict(R_sq=1e12, W=100.0, dy=0.5, **p), set=p,
                       window=[r["V_down"], r["V_up"]]))
    # B: end-fed resistive drain
    for W in (100.0, 300.0):
        for rD in (30.0, 100.0, 300.0, 1000.0):
            js.append(dict(name=f"B_W{W:g}_rD{rD:g}", kw=dict(R_sq=1e12, W=W, dy=0.5, r_D=rD, contacts="both", **BASE)))
    # C: distributed ballast, position-independent inhibition length lambda = sqrt(r_v / r_D)
    for lam in (5.0, 10.0, 20.0, 40.0):
        for RL in (0.0, 5e3):
            js.append(dict(name=f"C_lam{lam:g}_RL{RL:g}",
                           kw=dict(R_sq=1e12, W=300.0, dy=0.5, r_D=10.0, contacts="distributed", r_v=10.0 * lam**2,
                                   R_L=RL, **BASE)))
    # D: convergence of the R_sq = 1e10 case
    for dy in (0.1, 0.05):
        js.append(dict(name=f"D_Rsq1e10_dy{dy:g}", kw=dict(R_sq=1e10, W=100.0, dy=dy, **BASE)))
    return js


def run(job):
    import spatial as sp

    out = OUT / f"{job['name']}.json"
    if out.exists():
        return job["name"], "cached"
    kw = dict(job["kw"])
    kw.setdefault("vdd_up", (0.6, 3.4))
    kw.setdefault("dvdd", 0.04)
    r = sp.simulate(store_at=[("up", v) for v in (1.4, 1.8, 2.2, 2.6, 3.0)], **kw)
    up = [h for h in r["hist"] if h["dir"] == "up"]
    nmax_h = max(up, key=lambda h: h["n_fil"])
    cen = sorted(nmax_h["centres_um"])
    summ = dict(n_max_up=nmax_h["n_fil"], vdd_at_nmax=nmax_h["vdd"], widths_at_nmax=nmax_h["widths_um"],
                spacing_at_nmax=float(np.mean(np.diff(cen))) if len(cen) > 1 else None, events=r["events"],
                n_vs_vdd_up={f"{h['vdd']:.2f}": h["n_fil"] for h in up})
    out.write_text(json.dumps(dict(job=job, summary=summ, result=r)), encoding="utf-8")
    return job["name"], json.dumps(summ["events"]) + f" N_max={summ['n_max_up']}"


if __name__ == "__main__":
    js = jobs()
    print(len(js), "iş", flush=True)
    with Pool(4) as pool:
        for i, (name, msg) in enumerate(pool.imap_unordered(run, js), 1):
            print(f"[{i}/{len(js)}] {name}: {msg}", flush=True)
    print("TAMAM", flush=True)
