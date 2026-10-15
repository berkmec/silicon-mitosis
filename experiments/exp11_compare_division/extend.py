"""Review check: continue a finished 2D run from its last saved snapshot (float16) with the same parameters.
usage: python extend.py <run_name> <T_extra> [snapshot_every]
Writes runs/<run_name>__ext.json / .npz (time continues from the source run's last time)."""
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
sys.path.insert(0, str(Path(__file__).resolve().parent))
from qcell.active import ActiveCH  # noqa: E402
from scan import far_field_max  # noqa: E402

RUNS = ROOT / "results" / "exp11" / "runs"
name, T_extra = sys.argv[1], float(sys.argv[2])
snap_every = float(sys.argv[3]) if len(sys.argv) > 3 else 2500.0
src = json.loads((RUNS / f"{name}.json").read_text(encoding="utf-8"))
spec = src["spec"]
d = spec.get("d", np.inf)
d = np.inf if isinstance(d, str) else d
m = ActiveCH(n=spec["n"], L=spec["L"], S=1.0, A=spec.get("A", 0.0), d=d, dtype=np.float32, **spec["rates"])
z = np.load(RUNS / f"{name}.npz")
k_last = max(z.files, key=lambda s: int(s[1:]))
t0 = float(k_last[1:])
c = z[k_last].astype(np.float32)
hist, snaps = [], {}


def cb(t, c):
    dd = m.droplets(c)
    tt = t0 + t
    hist.append(dict(t=tt, n=dd["n"], rmax=max(dd["radii"]) if dd["radii"] else 0.0,
                     rmean=float(np.mean(dd["radii"])) if dd["radii"] else 0.0, perc=dd["perc"],
                     aspect=dd["aspect"], cmean=float(c.mean()), cfar=far_field_max(m, c)))
    if abs(t / snap_every - round(t / snap_every)) < 1e-9:
        snaps[f"t{int(round(tt))}"] = c.astype(np.float16)
        print(f"t={tt:.0f} N={dd['n']} R={[round(r, 1) for r in dd['radii']]} perc={dd['perc']}", flush=True)
        (RUNS / f"{name}__ext.json").write_text(json.dumps(dict(spec=spec, source=name, t_start=t0, hist=hist)),
                                                encoding="utf-8")
        np.savez_compressed(RUNS / f"{name}__ext.npz", **snaps)


m.run(c, T_extra, 0.2, every=int(round(250 / 0.2)), cb=cb)
print("TAMAM", flush=True)
