"""Exp14: search (nonlinear sharp-interface theory) for an active chemistry whose 6D droplet fits a 16^6-20^6 hypercube
(R* ~ 4-6.5), is division-prone (sigma_2(R*) > 0), keeps a dense core (c0 well above the spinodal) and has an outer
supply length sqrt(f''(eps)/k) >= 1.5 interface widths. Also reports d = 3, 4 for the same chemistry.
Writes results/exp14_dimension/hyper_scan.json."""
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
import itertools  # noqa: E402
import json  # noqa: E402
import time  # noqa: E402
import warnings  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from qcell import active_sharp_nl as nl  # noqa: E402

OUT = ROOT / "results" / "exp14_dimension" / "hyper_scan.json"


def drop(d, ch):
    R = nl.stable_radius(d, ch)
    if not np.isfinite(R):
        return None
    b = nl.base(R, d, ch)
    return dict(R=R, c0=b["c0"], sig2=nl.sigma(R, d, 2, 0.0, ch, b=b))


def main():
    t0, rows = time.time(), []
    for eps, k, nmk in itertools.product((0.08, 0.10, 0.12, 0.15), (0.02, 0.03, 0.05, 0.08), (0.3, 0.5, 0.7, 1.0)):
        lout = np.sqrt(nl.fpp(eps) / k)
        if lout < 1.5:
            continue
        ch = nl.Chem(eps, nmk * k, k)
        r6 = drop(6.0, ch)
        row = dict(eps=eps, k=k, nmk=nmk, l_out=float(lout), d6=r6)
        if r6 is not None and 3.5 <= r6["R"] <= 7.0:
            row["d3"], row["d4"] = drop(3.0, ch), drop(4.0, ch)
        rows.append(row)
        print(f"eps={eps} k={k} nmk={nmk} l_out={lout:.2f}: d6={r6} d3={row.get('d3')} ({time.time() - t0:.0f}s)",
              flush=True)
        OUT.write_text(json.dumps(rows, default=float), encoding="utf-8")
    print("TAMAM", flush=True)


if __name__ == "__main__":
    main()
