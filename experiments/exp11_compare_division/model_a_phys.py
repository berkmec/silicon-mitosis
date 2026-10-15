"""Model A (supercapacitor idea), physical: SKY130-like Si 2D electron layer under a gate of capacitance C = eps0 eps_g / t.
For each density n_s find the threshold capacitance C*(n_s) above which the uniform layer is unstable (division) within
q <= k_F: min_q [ K(n; C) + V_gate(q; C) ] = 0.  K includes kinetic + gate-screened HF exchange (+ correlation proxy)."""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from qcell import egas2d as G  # noqa: E402

OUT = ROOT / "results" / "exp11"
OUT.mkdir(parents=True, exist_ok=True)
M, EPS_SI = 0.19, 11.7
EPS0 = G.EPS0


def S_min(n, C_uF, eps_g, with_corr):
    t = EPS0 * eps_g / (C_uF * 1e-2)  # C in uF/cm^2 -> F/m^2 (x 1e-2)
    V = G.V_gate(EPS_SI, eps_g, t)
    K = G.inverse_compressibility(n, g=4, m_eff=M, V_x=V, eps_c=(EPS_SI + eps_g) / 2, with_corr=with_corr)
    kF = np.sqrt(np.pi * n)
    qs = np.linspace(kF / 100, kF, 120)
    return float(np.min(K + np.array([V(q) for q in qs])))


def C_star(n, eps_g, with_corr, lo=0.05, hi=3000.0):
    if S_min(n, hi, eps_g, with_corr) >= 0:
        return None  # stable even at the largest capacitance scanned
    if S_min(n, lo, eps_g, with_corr) < 0:
        return lo
    for _ in range(40):
        mid = np.sqrt(lo * hi)
        if S_min(n, mid, eps_g, with_corr) < 0:
            hi = mid
        else:
            lo = mid
    return float(np.sqrt(lo * hi))


ns_cm2 = np.logspace(10, 13, 25)
res = {"C_SKY130_uF_cm2": 0.8325, "EDL_range_uF_cm2": [1.0, 10.0], "curves": {}}
for eps_g in (3.9, 7.0, 20.0):
    for corr in (False, True):
        key = f"eps_g={eps_g},{'HF+korelasyon' if corr else 'HF'}"
        cs = [C_star(n * 1e4, eps_g, corr) for n in ns_cm2]
        res["curves"][key] = dict(n_cm2=ns_cm2.tolist(), C_star_uF_cm2=cs)
        at1e12 = cs[int(np.argmin(np.abs(ns_cm2 - 1e12)))]
        reach_edl = [float(n) for n, c in zip(ns_cm2, cs) if c is not None and c <= 10.0]
        print(f"{key:28s}: C*(1e12) = {('yok (≤3000)' if at1e12 is None else f'{at1e12:.1f} µF/cm²')};  "
              f"10 µF/cm² ile bölünen en yüksek n = {('yok' if not reach_edl else f'{max(reach_edl):.2g} cm⁻²')};  "
              f"0,83 µF/cm² (SKY130) ile: {('yok' if not [c for c in cs if c is not None and c <= 0.8325] else 'var')}", flush=True)
(OUT / "model_a_phys.json").write_text(json.dumps(res, indent=2), encoding="utf-8")
print("ok")
