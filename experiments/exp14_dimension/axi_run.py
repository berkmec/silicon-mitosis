"""Exp14: active droplet in d dimensions with the axisymmetric (z, r) reduction (src/qcell/active_axi.py).
usage: python axi_run.py <jobset>   (jobsets: see JOBS at the bottom; each job writes results/exp14_dimension/axi/<name>.json
and <name>.npz with snapshots). Runs at most 2 jobs in parallel at below-normal priority."""
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

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from qcell import active_sharp as th  # noqa: E402
from qcell.active_axi import AxiCH  # noqa: E402

OUT = ROOT / "results" / "exp14_dimension" / "axi"


def chem(eps, nu_m_over_k, k):
    nu_p, nu_m = eps * k, nu_m_over_k * k
    return dict(nu_p=nu_p, nu_m=nu_m, k_p=k, k_m=k)


def run(job):
    name, d, ch = job["name"], job["d"], job["chem"]
    eps, nmk, k = ch["eps"], ch["nmk"], ch["k"]
    m = AxiCH(d, job["nz"], job["nr"], h=job["h"], dt=job["dt"], S=1.0, **chem(eps, nmk, k))
    c = m.droplet(job["R0"], ecc=job.get("ecc", 0.05))
    c = c + 1e-3 * np.random.default_rng(job.get("seed", 1)).standard_normal(c.shape)
    T, every, snap_every = job["T"], job.get("every", 50.0), job.get("snap_every", 250.0)
    n_every, hist, snaps, t0 = int(round(every / m.dt)), [], {}, time.time()
    t_split, t_hollow, stop_after = None, None, job.get("stop_after_split")   # only the first event is valid here
    while m.t < T - 1e-9:
        c = m.step(c, n_every)
        bod, nb = m.bodies(c), m.bubbles(c)
        hist.append(dict(t=round(m.t, 3), n_axis=sum(b["axis"] for b in bod), n_ring=sum(not b["axis"] for b in bod),
                         n_bubble=nb, bodies=bod, cmin_axis=float(c[:, 0].min()), mass=m.mass(c)))
        if int(round(m.t / every)) % int(round(snap_every / every)) == 0:
            snaps[f"t{int(round(m.t))}"] = c.astype(np.float16)
        if len(bod) >= 2 and t_split is None:
            t_split = m.t
            snaps[f"split_t{int(round(m.t))}"] = c.astype(np.float16)
        if nb > 0 and t_hollow is None:
            t_hollow = m.t
            snaps[f"hollow_t{int(round(m.t))}"] = c.astype(np.float16)
        first = min(x for x in (t_split, t_hollow, np.inf) if x is not None)
        if not bod or (stop_after is not None and np.isfinite(first) and m.t >= first + stop_after):
            break
    res = dict(job=job, t_split=t_split, t_hollow=t_hollow, theory=dict(R_star=th.stable_radius(d, eps, nmk * k, k), R_l2=th.shape_threshold(d, 2, eps, nmk * k, k),
                                     R_tube=th.stable_radius(d - 1.0, eps, nmk * k, k) if d > 2 else th.stable_radius(1.0, eps, nmk * k, k)),
               hist=hist, wall=time.time() - t0)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / f"{name}.json").write_text(json.dumps(res, default=float), encoding="utf-8")
    np.savez_compressed(OUT / f"{name}.npz", **snaps)
    last = hist[-1] if hist else {}
    return f"{name}: t={last.get('t')} split={t_split} hollow={t_hollow} n_axis={last.get('n_axis')} " \
           f"n_ring={last.get('n_ring')} R_eq={[round(b['R_eq'], 2) for b in last.get('bodies', [])]} ({res['wall']:.0f} s)"


def grid(d, R_star, h=0.5):
    """Domain sized from the theory: r_max = R* + 4 l (l = reaction length), z long enough for a 2-3 fold elongation."""
    return dict(h=h, nr=int(np.ceil((1.6 * R_star + 12.0) / h)), nz=2 * int(np.ceil((7.0 * R_star + 20.0) / h / 2)))


CHEM_H = dict(eps=0.20, nmk=2.0, k=0.05)        # small structures, used for the Cartesian hypercube comparison
CHEM_Z = dict(eps=0.20, nmk=0.7, k=0.01)        # Zwicker et al. 2017 Fig. 3A (dnu = -0.001)


def jobs_h():
    out = []
    for d in (2.0, 3.0, 4.0, 5.0, 6.0):
        R = th.stable_radius(d, CHEM_H["eps"], CHEM_H["nmk"] * CHEM_H["k"], CHEM_H["k"])
        out.append(dict(name=f"H_d{d:g}", d=d, chem=CHEM_H, R0=max(0.9 * R, 2.0), T=3000.0, dt=0.25, every=25.0,
                        snap_every=250.0, **grid(d, R)))
    return out


def jobs_z(ds=(2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0), nmks=(0.7,), T=8000.0):
    out = []
    for nmk in nmks:
        for d in ds:
            ch = dict(CHEM_Z, nmk=nmk)
            R = th.stable_radius(d, ch["eps"], nmk * ch["k"], ch["k"])
            out.append(dict(name=f"Z_nmk{nmk:g}_d{d:g}", d=d, chem=ch, R0=5.0, T=T, dt=0.25, every=50.0,
                            snap_every=500.0, **grid(d, R)))
    return out


def jobs_dc():
    """Bracket the critical dimension of the first division for three turnover levels (dnu = -0.004, -0.001, +0.004)."""
    out = []
    plan = {0.4: (2.0, 2.25, 2.5, 3.0), 0.7: (2.1, 2.2, 2.3, 2.4, 3.75, 6.0), 1.2: (2.25, 2.5, 2.75, 3.0, 4.0, 5.0)}
    for nmk, ds in plan.items():
        for d in ds:
            ch = dict(CHEM_Z, nmk=nmk)
            R = th.stable_radius(d, ch["eps"], nmk * ch["k"], ch["k"])
            out.append(dict(name=f"Z_nmk{nmk:g}_d{d:g}", d=d, chem=ch, R0=5.0, T=12000.0, dt=0.25, every=50.0,
                            snap_every=500.0, stop_after_split=500.0, **grid(d, R)))
    return out


def jobs_6d():
    """Can 6D divide at all? Strong internal turnover (dnu = +0.012, +0.022) shifts the window upwards."""
    out = []
    for nmk in (2.0, 3.0):
        ch = dict(CHEM_Z, nmk=nmk)
        out.append(dict(name=f"Z_nmk{nmk:g}_d6", d=6.0, chem=ch, R0=8.0, T=8000.0, dt=0.25, every=50.0, snap_every=500.0,
                        stop_after_split=500.0, h=0.5, nz=240, nr=80))
    return out


def jobs_cmp():
    """3D vs 6D over the same turnover range: 3D at strong turnover, 6D between +0.004 and +0.012, 3D lower edge."""
    out = []
    for d, nmk, R0 in ((3.0, 2.0, 5.0), (3.0, 3.0, 5.0), (3.0, 0.5, 5.0), (6.0, 1.2, 8.0), (6.0, 1.6, 8.0)):
        ch = dict(CHEM_Z, nmk=nmk)
        R = th.stable_radius(d, ch["eps"], nmk * ch["k"], ch["k"])
        g = grid(d, R if np.isfinite(R) else 20.0)
        out.append(dict(name=f"Z_nmk{nmk:g}_d{d:g}", d=d, chem=ch, R0=R0, T=12000.0, dt=0.25, every=50.0,
                        snap_every=500.0, stop_after_split=500.0, **g))
    return out


JOBS = {"H": jobs_h, "Z": jobs_z, "DC": jobs_dc, "6D": jobs_6d, "CMP": jobs_cmp,
        "Zdnu": lambda: jobs_z(nmks=(0.6, 0.8, 1.0)),
        "Z3": lambda: jobs_z(ds=(3.0,))}

if __name__ == "__main__":
    js = JOBS[sys.argv[1]]()
    nproc = int(sys.argv[2]) if len(sys.argv) > 2 else 2
    with Pool(nproc) as p:
        for line in p.imap_unordered(run, js):
            print(line, flush=True)
    print("TAMAM", flush=True)
