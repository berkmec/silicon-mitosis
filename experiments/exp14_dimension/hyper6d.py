"""Exp14: genuine 6D Cahn-Hilliard on a 16^6 hypercube (16.8 M points, lean float32 solver of cart_run.py) versus the
O(5)-symmetric (z, r) reduction (qcell.active_axi) for the same active droplet: Zwicker chemistry (eps 0.2, nu-/k 0.7,
k 0.01), R0 = 4.5, dx = h = 1, dt = 0.25, t <= 40. Also the reduction at h = 0.5 (resolution effect).
Writes results/exp14_dimension/hyper6d.json."""
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
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import cart_run as C  # noqa: E402
from qcell.active_axi import AxiCH  # noqa: E402

CH = dict(nu_p=0.002, nu_m=0.007, k_p=0.01, k_m=0.01)
R0, T, DT, EVERY = 4.5, 40.0, 0.25, 20
OUT = ROOT / "results" / "exp14_dimension" / "hyper6d.json"
res = {}

# --- O(5)-symmetric reduction, h = 1 (same spacing as the hypercube) and h = 0.5
for h in (1.0, 0.5):
    m = AxiCH(6.0, int(16 / h), int(8 / h), h=h, dt=DT, S=1.0, **CH)
    c = m.droplet(R0, ecc=0.0)
    rows = []
    while m.t < T - 1e-9:
        c = m.step(c, EVERY)
        b = m.bodies(c)
        rows.append(dict(t=m.t, R_eq=b[0]["R_eq"] if b else 0.0, c_centre=float(c[int(8 / h), 0])))
    res[f"axi_h{h}"] = rows
    print(f"axi h={h}: R_eq(t=T)={rows[-1]['R_eq']:.3f}", flush=True)

# --- genuine 6D hypercube
t0 = time.time()
m = C.CartCH((16,) * 6, 1.0, DT, S=1.0, **CH)
c = m.droplet(R0, ecc=0.0)
print(f"6D field {c.nbytes / 1e6:.0f} MB, initial R_eq {C.components(c, 1.0)[0]['R_eq']:.3f}", flush=True)
rows = []
while m.t < T - 1e-9:
    c = m.step(c, EVERY)
    comp = C.components(c, 1.0)
    rows.append(dict(t=m.t, R_eq=comp[0]["R_eq"] if comp else 0.0, n=len(comp), c_centre=float(c[(8,) * 6]),
                     cmin=float(c.min()), cmax=float(c.max()), wall=time.time() - t0))
    print(f"6D t={m.t:5.1f} R_eq={rows[-1]['R_eq']:.3f} n={rows[-1]['n']} c_centre={rows[-1]['c_centre']:.3f} "
          f"({rows[-1]['wall']:.0f} s)", flush=True)
    res["cart6d"] = rows
    OUT.write_text(json.dumps(res, indent=1), encoding="utf-8")
np.save(ROOT / "results" / "exp14_dimension" / "hyper6d_slice.npy", c[:, :, 8, 8, 8, 8].astype(np.float16))
print("TAMAM", flush=True)
