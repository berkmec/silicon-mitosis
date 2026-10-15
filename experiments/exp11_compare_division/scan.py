"""Exp11 scans. usage: python scan.py {B|A|GS|ALL} [n_workers]
B  : dissipative (Zwicker active droplet), 2D, grid in (eps, dnu), A = 0.
A  : supercapacitor idea, 2D: gate-screened Coulomb (A, d) + electrical charging J (no leak) or J - lam c (leaky),
     plus equilibrium Coulomb baseline.
GS : Gray-Scott (Pearson 1993) (F, k) grid, single seed, spot count vs time.
Each job writes results/exp11/runs/<name>.json (+ .npz snapshots, float16)."""
import os
import sys

os.environ.setdefault("OMP_NUM_THREADS", "1")
if sys.platform == "win32":  # below-normal priority (inherited by pool workers): keep the user's desktop responsive
    import ctypes
    from ctypes import wintypes

    _k32 = ctypes.windll.kernel32
    _k32.GetCurrentProcess.restype = wintypes.HANDLE  # 64-bit pseudo-handle; default c_int truncates it
    _k32.SetPriorityClass.argtypes = [wintypes.HANDLE, wintypes.DWORD]
    _k32.SetPriorityClass(_k32.GetCurrentProcess(), 0x00004000)
import json  # noqa: E402
import time  # noqa: E402
from multiprocessing import Pool  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
from scipy import ndimage  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from qcell.active import ActiveCH, rates_charging, rates_zwicker  # noqa: E402
from qcell.grayscott import GrayScott  # noqa: E402

RUNS = ROOT / "results" / "exp11" / "runs"
RUNS.mkdir(parents=True, exist_ok=True)
SPINODAL = (3 - np.sqrt(3)) / 6  # f''(c) = 0 for f = c^2 (c-1)^2 / 2  -> 0.2113


def far_field_max(m, c, dist_min=6.0):
    """Largest c farther than dist_min from any droplet (c > 0.5); flags possible spinodal birth (c > 0.211)."""
    n = m.n
    big = np.tile(c > 0.5, (3, 3))
    if not big.any():
        return float(c.max())
    dist = ndimage.distance_transform_edt(~big)[n:2 * n, n:2 * n] * m.dx
    sel = dist > dist_min
    return float(c[sel].max()) if sel.any() else float("nan")


def job_ch(spec):
    name = spec["name"]
    out = RUNS / f"{name}.json"
    if out.exists():
        return name, "cached"
    t0 = time.time()
    m = ActiveCH(n=spec["n"], L=spec["L"], S=1.0, A=spec.get("A", 0.0), d=spec.get("d", np.inf),
                 dtype=np.float32, **spec["rates"])
    c = m.disk(spec.get("R0", 8.0), ecc=0.05, c_out=spec.get("c_out"))
    c = c + spec.get("noise", 1e-3) * np.random.default_rng(spec.get("seed", 1)).standard_normal(c.shape)
    T, dt, rec_dt, snap_dt = spec.get("T", 20000.0), 0.2, 250.0, 2500.0
    hist, snaps, state = [], {}, dict(dead=0, perc=0, why="T sonu")
    win = int(round(4000 / rec_dt))  # stationarity window: 4000 time units

    def cb(t, c):
        d = m.droplets(c)
        hist.append(dict(t=t, n=d["n"], rmax=max(d["radii"]) if d["radii"] else 0.0,
                         rmean=float(np.mean(d["radii"])) if d["radii"] else 0.0, perc=d["perc"],
                         aspect=d["aspect"], cmean=float(c.mean()), cfar=far_field_max(m, c)))
        if abs(t / snap_dt - round(t / snap_dt)) < 1e-9:
            snaps[f"t{int(round(t))}"] = c.astype(np.float16)
        state["dead"] = state["dead"] + 1 if (d["n"] == 0 and not d["perc"]) else 0
        state["perc"] = state["perc"] + 1 if d["perc"] else 0
        stop = None
        if state["dead"] >= 4:
            stop = "çözündü"
        elif state["perc"] >= 8:
            stop = "şerit"
        elif d["n"] >= 12:
            stop = "N>=12"
        elif t >= 4000 and len(hist) > win:
            old, new = hist[-1 - win], hist[-1]
            if old["n"] == new["n"] >= 1 and not new["perc"] \
                    and abs(new["rmax"] - old["rmax"]) <= 5e-3 * new["rmax"] \
                    and abs(new["rmean"] - old["rmean"]) <= 5e-3 * new["rmean"] \
                    and abs(new["aspect"] - old["aspect"]) <= 0.01:
                stop = "durağan"
        if stop:
            state["why"] = stop
            snaps[f"t{int(round(t))}"] = c.astype(np.float16)
            return False

    m.run(c, T, dt, every=int(round(rec_dt / dt)), cb=cb)
    np.savez_compressed(RUNS / f"{name}.npz", **snaps)
    res = dict(spec={k: (str(v) if isinstance(v, float) and np.isinf(v) else v) for k, v in spec.items()},
               hist=hist, wall_s=time.time() - t0, spinodal=SPINODAL, stop=state["why"], dtype="float32")
    out.write_text(json.dumps(res), encoding="utf-8")
    return name, f"{time.time() - t0:.0f}s"


def job_gs(spec):
    name = spec["name"]
    out = RUNS / f"{name}.json"
    if out.exists():
        return name, "cached"
    t0 = time.time()
    gs = GrayScott(n=128, F=spec["F"], k=spec["k"])
    U, V = gs.seed(np.random.default_rng(0))
    hist, snaps = [], {}

    def cb(t, U, V):
        hist.append(dict(t=t, **gs.spots(V)))
        if t % 5000 == 0:
            snaps[f"t{t}"] = V.astype(np.float16)

    gs.run(U, V, 20000, every=500, cb=cb)
    np.savez_compressed(RUNS / f"{name}.npz", **snaps)
    out.write_text(json.dumps(dict(spec=spec, hist=hist, wall_s=time.time() - t0)), encoding="utf-8")
    return name, f"{time.time() - t0:.0f}s"


def jobs_B():
    js = []
    for eps in (0.08, 0.11, 0.14, 0.17, 0.20):
        for dnu in (-0.006, -0.004, -0.002, 0.0, 0.002, 0.004):
            js.append(dict(kind="ch", name=f"B_e{eps:.2f}_d{dnu:+.3f}", n=256, L=128.0, eps=eps, dnu=dnu,
                           rates=rates_zwicker(eps, dnu)))
    for n, L in ((256, 128.0), (384, 192.0)):  # Zwicker Fig. 1C / 3A rates
        js.append(dict(kind="ch", name=f"B_paper_L{int(L)}", n=n, L=L, eps=0.19, dnu=-0.0011,
                       rates=rates_zwicker(0.19, -0.0011)))
    return js


def jobs_A():
    js = []
    for A, d in ((0.0, np.inf), (0.02, np.inf), (0.02, 3.0), (0.02, 1.0), (0.08, np.inf), (0.08, 3.0), (0.08, 1.0)):
        ds = "inf" if np.isinf(d) else f"{d:g}"
        js.append(dict(kind="ch", name=f"A_charge_A{A:g}_d{ds}", n=256, L=128.0, A=A, d=d, J=2e-5, lam=0.0,
                       rates=rates_charging(2e-5, 0.0), c_out=0.0))
        if A > 0:
            js.append(dict(kind="ch", name=f"A_leaky_A{A:g}_d{ds}", n=256, L=128.0, A=A, d=d, J=1.7e-3, lam=1e-2,
                           rates=rates_charging(1.7e-3, 1e-2)))
    for A in (0.02, 0.08):  # equilibrium baseline (no drive): one big charged droplet, conserved charge
        js.append(dict(kind="ch", name=f"A_eq_A{A:g}_dinf", n=256, L=128.0, A=A, d=np.inf, R0=16.0,
                       rates=rates_charging(0.0, 0.0), c_out=0.0))
    return js


def jobs_A2():
    """Weak-Coulomb regime (most favourable to the supercapacitor idea): the charged cell must first survive.
    Splitting into far-apart halves lowers Coulomb + line energy once R > 1.32 / sqrt(A) (d = inf estimate)."""
    js = []
    for A in (0.001, 0.003, 0.01):
        for d in (np.inf, 3.0, 1.0):
            ds = "inf" if np.isinf(d) else f"{d:g}"
            js.append(dict(kind="ch", name=f"A2_charge_A{A:g}_d{ds}", n=256, L=128.0, A=A, d=d, J=2e-5, lam=0.0,
                           rates=rates_charging(2e-5, 0.0), c_out=0.0))
    for d in (np.inf, 1.0):
        ds = "inf" if np.isinf(d) else f"{d:g}"
        js.append(dict(kind="ch", name=f"A2_leaky_A0.003_d{ds}", n=256, L=128.0, A=0.003, d=d, J=1.7e-3, lam=1e-2,
                       rates=rates_charging(1.7e-3, 1e-2)))
    for A in (0.003, 0.01):
        js.append(dict(kind="ch", name=f"A2_eq_A{A:g}_dinf", n=256, L=128.0, A=A, d=np.inf, R0=24.0,
                       rates=rates_charging(0.0, 0.0), c_out=0.0))
    return js


def jobs_A3():
    """Review (fairness): the supercapacitor's drive strength J was not scanned in A/A2. Charge until the mean rises by
    0.5 (T = 0.5 / J), for no Coulomb, weak Coulomb without/with a close gate, and a stronger Coulomb."""
    js = []
    for J in (1e-4, 5e-4, 2e-3):
        for A, d in ((0.0, np.inf), (0.003, np.inf), (0.003, 1.0), (0.01, np.inf)):
            ds = "inf" if np.isinf(d) else f"{d:g}"
            js.append(dict(kind="ch", name=f"A3_charge_J{J:g}_A{A:g}_d{ds}", n=256, L=128.0, A=A, d=d, J=J, lam=0.0,
                           rates=rates_charging(J, 0.0), c_out=0.0, T=0.5 / J))
    return js


def jobs_GS():
    return [dict(kind="gs", name=f"GS_F{F:.3f}_k{k:.3f}", F=F, k=k)
            for F in (0.010, 0.018, 0.026, 0.034, 0.042, 0.050, 0.058)
            for k in (0.045, 0.050, 0.055, 0.059, 0.063, 0.067)]


def jobs_R(names, seeds=(2, 3, 4), noise=2e-2):
    """Robustness: rerun named specs with 20x larger initial noise and other seeds."""
    base = {s["name"]: s for s in jobs_B() + jobs_A()}
    return [dict(base[nm], name=f"R_{nm}_s{sd}", seed=sd, noise=noise) for nm in names for sd in seeds]


def dispatch(spec):
    return job_gs(spec) if spec["kind"] == "gs" else job_ch(spec)


if __name__ == "__main__":
    which, nw = sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 6
    js = dict(B=jobs_B, A=jobs_A, GS=jobs_GS, A2=jobs_A2, A3=jobs_A3)
    if which == "R":  # python scan.py R n_workers name1,name2,...
        todo = jobs_R(sys.argv[3].split(","))
    else:
        todo = sum((js[w]() for w in (["B", "A", "GS", "A2"] if which == "ALL" else [which])), [])
    # longest single job first, then the cheap Gray-Scott runs, then A, then B from high to low supersaturation
    todo.sort(key=lambda s: (0,) if s.get("n", 0) == 384 else (1,) if s["kind"] == "gs"
              else (2,) if s["name"].startswith("A_") else (3, -s.get("eps", 0)))
    print(f"{len(todo)} iş, {nw} işçi", flush=True)
    with Pool(nw) as pool:
        for i, (name, msg) in enumerate(pool.imap_unordered(dispatch, todo), 1):
            print(f"[{i}/{len(todo)}] {name}: {msg}", flush=True)
    print("TAMAM", flush=True)
