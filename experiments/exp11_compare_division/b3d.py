"""Exp11 3D check: does my solver reproduce Zwicker et al. (2017) Fig. 1C/3A division in 3D?
Paper (supplement): c+(0)=0, c-(0)=1, k t0 = 1e-2, nu+ t0 = 1.9e-3 (Fig. 1C) / 2e-3 (Fig. 3A), nu- t0 = 7e-3,
cubic 50w box with 100 points for multiple divisions. Here: 64^3, L = 48 w (dx = 0.75 w, coarser than the paper's 0.5),
float32, periodic (paper: no-flux). usage: python b3d.py [T]"""
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
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from qcell.active import ActiveCH, rates_zwicker  # noqa: E402

OUT = ROOT / "results" / "exp11"
T = float(sys.argv[1]) if len(sys.argv) > 1 else 12000.0
EPS_ = float(sys.argv[2]) if len(sys.argv) > 2 else 0.20
DNU_ = float(sys.argv[3]) if len(sys.argv) > 3 else -0.001
TAG = sys.argv[4] if len(sys.argv) > 4 else ""  # "" keeps the original output names (b3d.json, b3d_slices.npz)
m = ActiveCH(n=64, L=48.0, dim=3, S=1.0, dtype=np.float32, **rates_zwicker(EPS_, DNU_))
c = m.disk(5.0, ecc=0.05) + 1e-3 * np.random.default_rng(1).standard_normal((64,) * 3)
hist, snaps, t0 = [], {}, time.time()


def cb(t, c):
    d = m.droplets(c)
    hist.append(dict(t=t, n=d["n"], radii=d["radii"], aspect=d["aspect"], perc=d["perc"], wall=time.time() - t0))
    if int(round(t)) % 1000 == 0:
        snaps[f"t{int(round(t))}"] = c[:, :, 32].astype(np.float16)  # mid-plane slice
        print(f"t={t:6.0f} N={d['n']} R={[round(r, 2) for r in d['radii']]} en-boy={d['aspect']:.2f} "
              f"({time.time() - t0:.0f} s)", flush=True)
        sfx = f"_{TAG}" if TAG else ""
        (OUT / f"b3d{sfx}.json").write_text(json.dumps(dict(hist=hist, eps=EPS_, dnu=DNU_, note=__doc__)),
                                            encoding="utf-8")
        np.savez_compressed(OUT / f"b3d{sfx}_slices.npz", **snaps)


m.run(c, T, 0.2, every=int(round(250 / 0.2)), cb=cb)
print("TAMAM", flush=True)
