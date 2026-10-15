"""Exp17: resumable 2D run (survives power cuts). The V_DD sweep is done in 1 V segments; after each segment the state
(body potentials), history and snapshots are written to results/exp17/runs2d/<name>_ckpt.npz / .json, and a restart
continues from the last saved V_DD (the mismatch pattern is regenerated from the same seed).
usage: python run2d_resume.py <name> <lam> [periodic 0/1] [seed_dvt V] [vdd_max]"""
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
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import spatial2d as s2  # noqa: E402

OUT = Path(__file__).resolve().parents[2] / "results" / "exp17" / "runs2d"
name, lam = sys.argv[1], float(sys.argv[2])
periodic = bool(int(sys.argv[3])) if len(sys.argv) > 3 else True
seed_dvt = float(sys.argv[4]) if len(sys.argv) > 4 else 0.03
vdd_max = float(sys.argv[5]) if len(sys.argv) > 5 else 16.0
KW = dict(Lx=540.0, Ly=540.0, dx=6.0, R_sq=2.5e10, R_L=1000.0, r_sqD=50.0, lam=lam, dvdd=0.1, seed=0, snap_every=0.5,
          max_it=400, verbose=True, periodic=periodic, seed_dvt=seed_dvt, seed_r=12.0)
ck_npz, ck_json = OUT / f"{name}_ckpt.npz", OUT / f"{name}_ckpt.json"
OUT.mkdir(parents=True, exist_ok=True)
if ck_json.exists() and ck_npz.exists():
    meta = json.loads(ck_json.read_text(encoding="utf-8"))
    z = np.load(ck_npz)
    V, hist, snaps, v_next = z["V_state"], meta["hist"], {k: z[k] for k in z.files if k != "V_state"}, meta["v_next"]
    print(f"devam: V_DD = {v_next:.2f}", flush=True)
else:
    V, hist, snaps, v_next = None, [], {}, 1.2
while v_next <= vdd_max + 1e-9:
    v_end = min(round(v_next + 0.9, 4), vdd_max)
    r = s2.simulate2d_sm(vdd=(v_next, v_end), V0=V, **KW)
    V, v_next = r["V_final"], round(v_end + 0.1, 4)
    hist += r["hist"]
    snaps.update(r["snaps"])
    np.savez_compressed(OUT / f"{name}_ckpt_tmp.npz", V_state=V, **snaps)
    os.replace(OUT / f"{name}_ckpt_tmp.npz", ck_npz)
    ck_json.write_text(json.dumps(dict(params=r["params"], hist=hist, v_next=v_next, done=v_next > vdd_max)), encoding="utf-8")
    print(f"kayıt: V_DD {v_end:.2f}'a kadar", flush=True)
print("TAMAM", flush=True)
