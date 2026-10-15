"""Review check for the Gray-Scott ramp claim "existing spots divide into smaller ones".
Settle stage 0 (D x 1), then halve D and track spot lineage for 3000 time units with snapshots every 25:
division = a spot at t overlaps >= 2 spots at t + dt; birth = a spot at t + dt overlaps no spot at t;
death = a spot at t overlaps nothing at t + dt; merge = a spot at t + dt overlaps >= 2 spots at t."""
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
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
from scipy import ndimage  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from qcell.grayscott import GrayScott  # noqa: E402

OUT = ROOT / "results" / "exp11"


def plabel(mask):
    """Connected components on a periodic grid (labels merged across opposite edges)."""
    lab, k = ndimage.label(mask)
    parent = list(range(k + 1))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for a, b in list(zip(lab[0, :], lab[-1, :])) + list(zip(lab[:, 0], lab[:, -1])):
        if a and b:
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[rb] = ra
    roots = np.array([find(i) for i in range(k + 1)])
    _, inv = np.unique(roots, return_inverse=True)  # root 0 (background) maps to 0
    return inv[lab], int(inv.max())


def events(l0, k0, l1, k1):
    pairs = set(zip(l0[(l0 > 0) & (l1 > 0)].tolist(), l1[(l0 > 0) & (l1 > 0)].tolist()))
    ch = {p: 0 for p in range(1, k0 + 1)}
    pa = {c: 0 for c in range(1, k1 + 1)}
    for p, c in pairs:
        ch[p] += 1
        pa[c] += 1
    return dict(div=sum(v >= 2 for v in ch.values()), death=sum(v == 0 for v in ch.values()),
                birth=sum(v == 0 for v in pa.values()), merge=sum(v >= 2 for v in pa.values()))


n, Lbox = 256, 1.25
gs = GrayScott(n=n, h=Lbox / n, F=0.026, k=0.059)
U, V = gs.seed(np.random.default_rng(0), half=10)
dt0 = min(1.0, 0.2 * gs.h**2 / gs.Du)
U, V = gs.run(U, V, int(round(12000 / dt0)), dt=dt0)  # stage 0, same as gs_ramp.py
gs.Du, gs.Dv = gs.Du / 2, gs.Dv / 2
dt1 = min(1.0, 0.2 * gs.h**2 / gs.Du)
every = int(round(25 / dt1))
lab, k = plabel(V > 0.1)
tot = dict(div=0, death=0, birth=0, merge=0)
hist = [dict(t=0.0, n=k)]
for i in range(1, int(round(3000 / 25)) + 1):
    U, V = gs.run(U, V, every, dt=dt1)
    lab1, k1 = plabel(V > 0.1)
    ev = events(lab, k, lab1, k1)
    for key in tot:
        tot[key] += ev[key]
    hist.append(dict(t=i * 25.0, n=k1, **ev))
    lab, k = lab1, k1
res = dict(F=0.026, k=0.059, n=n, L=Lbox, window=3000, snap_dt=25, N_start=hist[0]["n"], N_end=hist[-1]["n"],
           totals=tot, net_change=hist[-1]["n"] - hist[0]["n"], hist=hist)
(OUT / "gs_track.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
print(f"N: {res['N_start']} -> {res['N_end']} (net {res['net_change']:+d}); olaylar: {tot}", flush=True)
print("Denge: bölünme - ölüm + doğum - birleşme =", tot["div"] - tot["death"] + tot["birth"] - tot["merge"], flush=True)
