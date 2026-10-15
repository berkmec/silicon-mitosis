"""Exp14: passive Rayleigh-Plateau threshold in d dimensions with the axisymmetric Cahn-Hilliard solver.
Tube R0 = 8 (interface width 1), far field at the Gibbs-Thomson value, wall at 1.5 R0 (small reservoir so that the
uniform mode stays stable); growth rate of the varicose mode A(t) = (R(z=0) - R(z=Lz/2)) / 2 fitted for several kz.
Writes results/exp14_dimension/passive_rp.json (sigma(kz R) per d and the interpolated zero crossing vs sqrt(d-2))."""
import json
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from qcell.active_axi import AxiCH  # noqa: E402

LG, R0, H = 1.0 / 12.0, 8.0, 0.5


def iface(c):
    R = np.zeros(c.shape[0])
    for j in range(c.shape[0]):
        idx = np.nonzero((c[j, :-1] >= 0.5) & (c[j, 1:] < 0.5))[0]
        if len(idx):
            i = idx[-1]
            R[j] = (i + 0.5) * H + H * (c[j, i] - 0.5) / (c[j, i] - c[j, i + 1])
    return R


def rate(d, f):
    kz = f * np.sqrt(d - 2.0) / R0
    nz = int(round(2 * np.pi / kz / H))
    nz += nz % 2
    m = AxiCH(d, nz, int(round(1.5 * R0 / H)), h=H, dt=0.5, S=1.0)
    c = m.tube(R0, amp=0.02, c_far=LG * (d - 2.0) / R0)
    ts, A = [], []
    for _ in range(40):
        Rz = iface(c)
        ts.append(m.t)
        A.append(0.5 * (Rz[0] - Rz[nz // 2]))
        c = m.step(c, 100)
    ts, A = np.array(ts), np.array(A)
    sel = (ts >= 300) & (A > 0)
    g = float(np.polyfit(ts[sel], np.log(A[sel]), 1)[0]) if sel.sum() >= 5 else float("nan")
    return float(2 * np.pi / (nz * H) * R0), g


out = {}
for d in (3.0, 4.0, 5.0, 6.0):
    rows = [rate(d, f) for f in (0.6, 0.8, 0.9, 0.95, 1.05, 1.1, 1.2)]
    x, y = np.array(rows).T
    i = np.nonzero((y[:-1] > 0) & (y[1:] <= 0))[0]
    xc = float(x[i[0]] + (x[i[0] + 1] - x[i[0]]) * y[i[0]] / (y[i[0]] - y[i[0] + 1])) if len(i) else float("nan")
    out[str(d)] = dict(kR=x.tolist(), sigma=y.tolist(), kR_c=xc, theory=float(np.sqrt(d - 2.0)))
    print(f"d={d}: measured kR_c = {xc:.3f}, sqrt(d-2) = {np.sqrt(d - 2):.3f}, sigma = {[f'{v:+.1e}' for v in y]}", flush=True)
(ROOT / "results" / "exp14_dimension" / "passive_rp.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
