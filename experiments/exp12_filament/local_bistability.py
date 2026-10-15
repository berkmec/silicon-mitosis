"""Exp12 step 3: is the floating-body / impact-ionisation feedback of a SKY130 nFET BISTABLE (S-shaped) at 4 K?
Local physics from the official BSIM4 card (ngspice tables spice/tab_vb*.dat: channel current id and isub per um of
width, V_bs 0..1.0 V, V_gs 0.3..1.8, V_ds 0..2.6). 4 K knobs (assumptions): dvt (V_T increase, applied as V_g - dvt),
fii (impact ionisation x fii), body diode j_dio = j_on exp((V_b - V_on)/V_s), frozen-substrate leak g_sub V_b.
Steady states of the body node at fixed (V_g, V_d):  F(V_b) = fii*j_ii(V_g-dvt, V_d, V_b) - j_dio(V_b) - g_sub V_b = 0.
Counts the roots on V_b in [0, 1.2] (V_b > 1.0 uses the V_b = 1.0 table: BSIM4 clamps V_bseff below Phi_s anyway).
Also reports the 'feedback slope' d ln j_ii / d V_b, which must exceed the leak's d ln j / d V_b somewhere for >1 root."""
import itertools
import json
import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator

ROOT = Path(__file__).resolve().parents[2]
SP = Path(__file__).resolve().parent / "spice"
OUT = ROOT / "results" / "exp12"
OUT.mkdir(parents=True, exist_ok=True)

VB = np.round(np.arange(0.0, 1.0001, 0.1), 2)
VG = np.round(np.arange(0.3, 1.8001, 0.05), 3)
VD = np.round(np.arange(0.0, 2.6001, 0.05), 3)


def load_tables():
    ID = np.zeros((len(VB), len(VG), len(VD)))
    II = np.zeros_like(ID)
    for i, vb in enumerate(VB):
        d = np.loadtxt(SP / f"tab_vb{vb:.1f}.dat")
        ID[i] = d[:, 1].reshape(len(VG), len(VD))
        II[i] = d[:, 3].reshape(len(VG), len(VD))
    return ID * 1e0, II  # A per um (tables are for W = 1 um)


ID, II = load_tables()
f_id = RegularGridInterpolator((VB, VG, VD), ID, bounds_error=False, fill_value=None)
f_lii = RegularGridInterpolator((VB, VG, VD), np.log(np.maximum(II, 1e-30)), bounds_error=False, fill_value=None)


def j_ii(vg, vd, vb):
    vb = np.clip(vb, 0.0, 1.0)
    pts = np.column_stack([vb, np.full_like(vb, vg), np.full_like(vb, vd)])
    return np.exp(f_lii(pts))


def j_ch(vg, vd, vb):
    vb = np.clip(vb, 0.0, 1.0)
    pts = np.column_stack([vb, np.full_like(vb, vg), np.full_like(vb, vd)])
    return f_id(pts)


def roots(vg, vd, fii, dvt, von, vs, g_sub, j_on=1e-6):
    vb = np.linspace(0.0, 1.2, 1201)
    F = fii * j_ii(vg - dvt, vd, vb) - j_on * np.exp((vb - von) / vs) - g_sub * vb
    s = np.sign(F)
    idx = np.where(s[:-1] * s[1:] < 0)[0]
    return [float(vb[i] - F[i] * (vb[i + 1] - vb[i]) / (F[i + 1] - F[i])) for i in idx]


if __name__ == "__main__":
    knobs = dict(dvt=(0.0, 0.15, 0.3), fii=(1, 10, 100), von=(0.8, 1.0, 1.1), vs=(0.005, 0.025, 0.05),
                 g_sub=(1e-15, 1e-12, 1e-9))  # g_sub in A/(V um): frozen-substrate leak per um of width
    vg_list, vd_grid = (1.2, 1.8), np.round(np.arange(0.5, 2.6001, 0.05), 3)
    bistable, total, examples = 0, 0, []
    for vg in vg_list:
        for dvt, fii, von, vs, gs in itertools.product(*knobs.values()):
            total += 1
            nroots = [len(roots(vg, vd, fii, dvt, von, vs, gs)) for vd in vd_grid]
            if max(nroots) >= 3:
                bistable += 1
                vd_b = [float(v) for v, k in zip(vd_grid, nroots) if k >= 3]
                examples.append(dict(vg=vg, dvt=dvt, fii=fii, von=von, vs=vs, g_sub=gs, vd_bistable=[min(vd_b), max(vd_b)]))
    # feedback strength: d ln j_ii / d V_b at the nominal point, versus the leak's slope 1/V_s of the diode
    vb = np.linspace(0.0, 1.0, 101)
    slopes = {}
    for vg, vd in ((1.2, 1.2), (1.2, 1.8), (1.8, 1.8), (1.2, 2.4)):
        lj = np.log(j_ii(vg, vd, vb))
        slopes[f"vg{vg}_vd{vd}"] = dict(max_dlnj_dVb=float(np.max(np.gradient(lj, vb))),
                                        jii_ratio_Vb1_over_Vb0=float(np.exp(lj[-1] - lj[0])),
                                        jch_ratio_Vb1_over_Vb0=float(j_ch(vg, vd, np.array([1.0]))[0] /
                                                                     j_ch(vg, vd, np.array([0.0]))[0]))
    out = dict(knob_grid=knobs, vg=vg_list, n_cases=total, n_bistable=bistable, examples=examples[:20],
               feedback=slopes, note=__doc__)
    (OUT / "local_bistability.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    print(f"bistabil düğme kombinasyonu: {bistable}/{total}")
    for e in examples[:8]:
        print("  ", e)
    print(json.dumps(slopes, indent=1))
