"""Exp17: held-out test of the filament split law (law.py) with the exp12/13 spatial model (distributed drain ballast).
Each job: up sweep only, dV_DD = 0.01 V; split events are extracted from consecutive filament lists (parent width just
before the split, daughters, V_DD, current, offset of the cut from the parent centre).
usage: python run.py [nproc]  -> results/exp17/runs/<name>.json"""
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

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "experiments" / "exp12_filament"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import law  # noqa: E402
import spatial as sp  # noqa: E402

OUT = ROOT / "results" / "exp17" / "runs"


def split_events(hist):
    ev, prev = [], None
    for h in hist:
        cur = [(c - w / 2, c + w / 2) for c, w in zip(h["centres_um"], h["widths_um"])]
        if prev is not None:
            for (s1, e1) in prev:
                kids = [(s2, e2) for (s2, e2) in cur if min(e1, e2) - max(s1, s2) > 0]
                if len(kids) >= 2:
                    ev.append(dict(vdd=h["vdd"], I_uA=h["I_uA"], parent_w=e1 - s1, parent_c=(s1 + e1) / 2,
                                   kids_w=[e - s for s, e in kids], cut_offset=(kids[0][1] + kids[1][0]) / 2 - (s1 + e1) / 2,
                                   n_fil_after=len(cur)))
        prev = cur
    return ev


def run(job):
    t0 = time.time()
    p = job.get("p", law.BASE)
    kw = dict(W=job["W"], dy=0.5, R_sq=1e12, R_L=job["R_L"], vdd_up=(0.6, job["vdd_max"]), dvdd=job.get("dvdd", 0.01 if job["R_L"] < 4e3 else 0.02), seed=job["seed"],
              r_D=job["r_D"], contacts="distributed", r_v=job["r_D"] * job["lam"] ** 2, down=False,
              sigma_scale=job["sigma_scale"], t_relax=job.get("t_relax", 1.0), **p)
    r = sp.simulate(**kw)
    c = law.constants(p)
    ev = split_events(r["hist"])
    res = dict(job=job, prediction=dict(lambda_c=law.lambda_c(job["r_D"], c), w_split=law.w_split(job["lam"], job["r_D"], c),
                                        vdd_split=law.vdd_at_split(job["lam"], job["r_D"], job["R_L"], c)),
               splits=ev, n_max=max(h["n_fil"] for h in r["hist"]), hist=r["hist"], wall=time.time() - t0)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / f"{job['name']}.json").write_text(json.dumps(res), encoding="utf-8")
    first = ev[0] if ev else None
    return (f"{job['name']}: öngörü w_split={res['prediction']['w_split']:.1f} | splits={len(ev)} n_max={res['n_max']} "
            f"ilk: {('w=%.1f vdd=%.2f kids=%s' % (first['parent_w'], first['vdd'], [round(k, 1) for k in first['kids_w']])) if first else '-'}"
            f" ({res['wall']:.0f} s)")


def jobs():
    out = []
    base = dict(W=600.0, R_L=2e3, vdd_max=2.6, seed=0, sigma_scale=1.0, r_D=10.0)
    for lam in (30, 36, 40, 45, 50, 60, 80, 120):
        out.append(dict(base, name=f"A_lam{lam}", lam=float(lam)))
    for lam in (40, 60, 80):
        out.append(dict(base, name=f"B_lowmm_lam{lam}", lam=float(lam), sigma_scale=0.1))
    for r_D, lam in ((5.0, 60), (20.0, 60), (5.0, 80)):
        out.append(dict(base, name=f"C_rD{r_D:g}_lam{lam}", lam=float(lam), r_D=r_D))
    for seed in (1, 2):
        out.append(dict(base, name=f"D_lam60_s{seed}", lam=60.0, seed=seed))
    return out


def jobs_r5k():
    """Clean protocol: R_L = 5 kOhm so that the filament is born much narrower than w_split (35-55 um) and grows slowly
    (dV_DD = 0.02 V -> 4 uA, ~1.6 um per step)."""
    out = []
    base = dict(W=600.0, R_L=5e3, vdd_max=5.6, seed=0, sigma_scale=1.0, r_D=10.0)
    for lam in (40, 45, 50, 60, 80, 120):
        out.append(dict(base, name=f"E_lam{lam}", lam=float(lam)))
    for lam in (40, 50, 60, 80, 120):
        out.append(dict(base, name=f"F_lowmm_lam{lam}", lam=float(lam), sigma_scale=0.1))
    for r_D, lam in ((5.0, 60), (5.0, 80), (20.0, 40), (20.0, 60)):
        out.append(dict(base, name=f"G_rD{r_D:g}_lam{lam}", lam=float(lam), r_D=r_D))
    for seed in (1, 2):
        out.append(dict(base, name=f"H_lam60_s{seed}", lam=60.0, seed=seed))
    return out


U1 = dict(ss=0.03, vs=0.01, fii=3, vov0=-0.1, g_sub=1e-12, von=1.05)
U2 = dict(ss=0.015, vs=0.005, fii=30, vov0=-0.1, g_sub=1e-13, von=1.05)


def jobs_univ():
    """Universality: two device-parameter sets far from the base set; lambda = 1.5 and 2.5 lambda_c (r_D = 10)."""
    out = []
    base = dict(W=600.0, R_L=5e3, vdd_max=5.6, seed=0, sigma_scale=1.0, r_D=10.0)
    for tag, p in (("U1", U1), ("U2", U2)):
        lc = law.lambda_c(10.0, law.constants(p))
        for x in (1.5, 2.5):
            out.append(dict(base, name=f"{tag}_x{x:g}", lam=float(round(x * lc, 1)), p=p))
    return out


NM = dict(ss=0.02, vs=0.005, fii=10, vov0=-0.1, g_sub=1e-15, von=1.05)   # exp13 set A3: no Maxwell point in the window


def jobs_nomax():
    """'No-Maxwell' regime (g_sub = 1e-15): fronts always advance until V_D = V_down, so dV_c ~ 0 and a ballasted
    filament should fragment (comb of narrow filaments) instead of growing to w_split. Control: no ballast."""
    base = dict(W=600.0, R_L=5e3, vdd_max=5.6, seed=0, sigma_scale=1.0, r_D=10.0, p=NM)
    return [dict(base, name="N_ballast_lam60", lam=60.0), dict(base, name="N_noballast", lam=60.0, r_D=1e-9)]


def jobs_slow():
    """Kinetic test: g_sub = 1e-15 gives a body discharge time c_b / g_sub ~ 1 s, equal to the 1 s per V_DD step of the
    standard sweep. Sweep 30x slower (t_relax = 30 s per step)."""
    base = dict(W=600.0, R_L=5e3, vdd_max=3.6, seed=0, sigma_scale=1.0, r_D=10.0, p=NM, t_relax=30.0)
    return [dict(base, name="N_ballast_lam60_slow", lam=60.0)]


if __name__ == "__main__":
    mode = sys.argv[2] if len(sys.argv) > 2 else ""
    js = {"r5k": jobs_r5k, "univ": jobs_univ, "nomax": jobs_nomax, "slow": jobs_slow}.get(mode, jobs)()
    nproc = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    with Pool(nproc) as p:
        for line in p.imap_unordered(run, js):
            print(line, flush=True)
    print("TAMAM", flush=True)
