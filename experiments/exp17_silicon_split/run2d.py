"""Exp17: 2D current spots in an array of nFET finger segments with a 2D drain sheet + distributed ballast
(spatial2d.py). Prediction (isolated disk, centre quench): Delta_u(R) = r_sqD lam^2 J x K1(x) (I0(x) - 1), x = R / lam,
J = j_on / dx; split (centre quench -> ring) when Delta_u = V_M - V_down.
usage: python run2d.py <name> [nproc]  -> results/exp17/runs2d/<name>.json|npz"""
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
from multiprocessing import Pool  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
from scipy import optimize, special  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import law  # noqa: E402
import spatial2d as s2  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "exp17" / "runs2d"


def r_crit(lam, r_sqD, dx, c):
    J = c["j_on"] / dx
    y = c["dV_c"] / (r_sqD * lam**2 * J)
    g = lambda x: x * special.k1(x) * (special.i0(x) - 1) - y
    xs = np.linspace(1e-3, 40, 40000)
    v = g(xs)
    i = np.nonzero(np.sign(v[:-1]) != np.sign(v[1:]))[0]
    return float(optimize.brentq(g, xs[i[0]], xs[i[0] + 1]) * lam) if len(i) else float("inf")


def run(job):
    t0 = time.time()
    c = law.constants()
    sim = s2.simulate2d_sm if job.get("sm") else s2.simulate2d
    kw = dict(max_it=job.get("max_it", 600)) if job.get("sm") else {}
    r = sim(Lx=job["L"], Ly=job["L"], dx=job["dx"], R_sq=job.get("R_sq", 5e10), R_L=job["R_L"], r_sqD=job["r_sqD"],
            lam=job["lam"], vdd=(job.get("vdd_min", 0.6), job["vdd_max"]), dvdd=job["dvdd"], seed=job["seed"],
            snap_every=job.get("snap_every", 0.1), **kw)
    OUT.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(OUT / f"{job['name']}.npz", **r.pop("snaps"))
    res = dict(job=job, prediction=dict(R_c=r_crit(job["lam"], job["r_sqD"], job["dx"], c)), **r, wall=time.time() - t0)
    (OUT / f"{job['name']}.json").write_text(json.dumps(res), encoding="utf-8")
    nmax = max(h["n"] for h in r["hist"])
    rings = sum(1 for h in r["hist"] for s in h["spots"] if s["holes"])
    return f"{job['name']}: R_c(öngörü)={res['prediction']['R_c']:.1f} n_max={nmax} halka-adımı={rings} ({res['wall']:.0f} s)"


JOBS = {
    "S64": dict(lam=64.0), "S45": dict(lam=45.0), "S128": dict(lam=128.0), "S25": dict(lam=25.0),
    # converged solver (Sherman-Morrison load line), large domain L = 600 um, dx = 6 um, slower body coupling R_sq
    "T64": dict(lam=64.0, sm=True, L=600.0, dx=6.0, R_sq=2.5e10, R_L=1000.0, vdd_min=1.3, vdd_max=16.0, dvdd=0.1,
                snap_every=0.5),
    "T30": dict(lam=30.0, sm=True, L=600.0, dx=6.0, R_sq=2.5e10, R_L=1000.0, vdd_min=1.3, vdd_max=16.0, dvdd=0.1,
                snap_every=0.5),
}
BASE = dict(L=400.0, dx=4.0, R_L=300.0, r_sqD=50.0, vdd_max=7.0, dvdd=0.02, seed=0)

if __name__ == "__main__":
    names = sys.argv[1].split(",")
    nproc = int(sys.argv[2]) if len(sys.argv) > 2 else 2
    js = [dict(BASE, **dict(JOBS[n], name=n)) for n in names]
    with Pool(nproc) as p:
        for line in p.imap_unordered(run, js):
            print(line, flush=True)
    print("TAMAM", flush=True)
