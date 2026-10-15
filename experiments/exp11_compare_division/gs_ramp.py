"""Exp11 exploratory: 'dividing while shrinking' in Gray-Scott. Ramping the drive = making the reactions faster
relative to diffusion, i.e. D_u, D_v -> D/2 at each stage (spot size ~ sqrt(D)). Do existing spots divide into
smaller ones, stage after stage? F = 0.026, k = 0.059 (self-replicating band of the exp11 GS scan)."""
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

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from qcell.grayscott import GrayScott  # noqa: E402

OUT = ROOT / "results" / "exp11"
n, Lbox = 256, 1.25
gs = GrayScott(n=n, h=Lbox / n, F=0.026, k=0.059)
U, V = gs.seed(np.random.default_rng(0), half=10)
D0u, D0v = gs.Du, gs.Dv
stages = [(0, 12000), (1, 8000), (2, 8000), (3, 8000)]
hist, snaps = [], {}
for s, T in stages:
    gs.Du, gs.Dv = D0u / 2**s, D0v / 2**s
    dt = min(1.0, 0.2 * gs.h**2 / gs.Du)  # explicit Euler: Du dt / h^2 <= 0.2
    steps = int(round(T / dt))
    t_off = sum(t for _, t in stages[:s])

    def cb(t, U, V, s=s, dt=dt, t_off=t_off):
        sp = gs.spots(V)
        area = sp["area_frac"] * Lbox**2 / sp["n"] if sp["n"] else 0.0
        hist.append(dict(t=t_off + t * dt, stage=s, D_scale=2.0**-s, n=sp["n"], perc=sp["perc"],
                         area_frac=sp["area_frac"], mean_spot_area=area))

    U, V = gs.run(U, V, steps, dt=dt, every=max(1, int(round(500 / dt))), cb=cb)
    snaps[f"stage{s}"] = V.astype(np.float16)
    h = hist[-1]
    print(f"evre {s}: D×{2.0**-s:.3g}, dt={dt:.3g} -> N={h['n']}, ort. leke alanı={h['mean_spot_area']:.3g} "
          f"(√alan={np.sqrt(h['mean_spot_area']):.3g}), şerit={h['perc']}", flush=True)
np.savez_compressed(OUT / "gs_ramp.npz", **snaps)
(OUT / "gs_ramp.json").write_text(json.dumps(dict(F=0.026, k=0.059, n=n, L=Lbox, hist=hist)), encoding="utf-8")
print("TAMAM", flush=True)
