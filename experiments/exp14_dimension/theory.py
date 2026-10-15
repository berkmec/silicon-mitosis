"""Exp14 theory tables: linear (active_sharp) and nonlinear (active_sharp_nl) sharp-interface theory of the active droplet
of qcell.active in d dimensions; Zwicker chemistry eps = 0.2, k = 0.01 and several dnu.
For each (dnu, d): stable droplet radius R*, centre concentration, sigma_l(R*) l = 2..4; tube (m = d-1) radius R_t and
the max over kz of the varicose (l = 0) and sinuous (l = 1) growth rates. Then the critical dimension d_c(dnu) of the tube's
varicose instability (bisection on d). Writes results/exp14_dimension/theory.json."""
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
import warnings  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from qcell import active_sharp as lin  # noqa: E402
from qcell import active_sharp_nl as nl  # noqa: E402

OUT = ROOT / "results" / "exp14_dimension" / "theory.json"
EPS, K = 0.20, 0.01
DNUS = (-0.004, -0.002, -0.001, 0.0, 0.002, 0.004)
DS = (2.0, 2.25, 2.5, 2.75, 3.0, 3.5, 4.0, 5.0, 6.0)
KZ = np.concatenate([np.linspace(0.01, 0.3, 30), np.linspace(0.32, 1.2, 23)])


def tube_nl(d, ch):
    m = d - 1.0
    R = nl.stable_radius(m, ch)
    if not np.isfinite(R):
        return dict(R=None)
    b = nl.base(R, m, ch)
    var = np.array([nl.sigma(R, m, 0, kz, ch, b=b) for kz in KZ])
    sin = np.array([nl.sigma(R, m, 1, kz, ch, b=b) for kz in KZ])
    i = int(np.nanargmax(var))
    return dict(R=R, c0=b["c0"], var_max=float(var[i]), kz_var=float(KZ[i]), sin_max=float(np.nanmax(sin)),
                var=var.tolist(), sin=sin.tolist())


def tube_lin(d, num):
    R = lin.stable_radius(d - 1.0, EPS, num, K)
    if not np.isfinite(R):
        return dict(R=None)
    var, sin = lin.sigma(R, d - 1.0, 0, KZ, EPS, num, K), lin.sigma(R, d - 1.0, 1, KZ, EPS, num, K)
    i = int(np.argmax(var))
    return dict(R=R, var_max=float(var[i]), kz_var=float(KZ[i]), sin_max=float(sin.max()))


def droplet_nl(d, ch):
    R = nl.stable_radius(d, ch)
    if not np.isfinite(R):
        return dict(R=None)
    b = nl.base(R, d, ch)
    return dict(R=R, c0=b["c0"], sig={l: nl.sigma(R, d, l, 0.0, ch, b=b) for l in (2, 3, 4)})


def d_crit(ch, lo=2.0, hi=3.5, tol=0.02):
    f = lambda d: tube_nl(d, ch).get("var_max", np.nan)
    flo, fhi = f(lo), f(hi)
    if not (np.isfinite(flo) and np.isfinite(fhi)):
        return None
    if flo > 0:
        return f"<{lo}"
    if fhi < 0:
        return f">{hi}"
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        if f(mid) > 0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def main():
    t0, res = time.time(), dict(eps=EPS, k=K, kz=KZ.tolist(), rows=[], d_c={})
    for dnu in DNUS:
        num = dnu + K - EPS * K
        ch = nl.Chem(EPS, num, K)
        for d in DS:
            row = dict(dnu=dnu, d=d, lin_R=lin.stable_radius(d, EPS, num, K), lin_tube=tube_lin(d, num),
                       nl_drop=droplet_nl(d, ch), nl_tube=tube_nl(d, ch))
            res["rows"].append(row)
            dr, tb = row["nl_drop"], row["nl_tube"]
            print(f"dnu={dnu:+.3f} d={d}: R*nl={dr.get('R')} sig2={dr.get('sig', {}).get(2)} | tube R={tb.get('R')} "
                  f"var={tb.get('var_max')} sin={tb.get('sin_max')} ({time.time() - t0:.0f}s)", flush=True)
            OUT.write_text(json.dumps(res, default=float), encoding="utf-8")
        res["d_c"][str(dnu)] = d_crit(ch)
        print(f"d_c(dnu={dnu}) = {res['d_c'][str(dnu)]}", flush=True)
        OUT.write_text(json.dumps(res, default=float), encoding="utf-8")
    print("TAMAM", flush=True)


if __name__ == "__main__":
    main()
