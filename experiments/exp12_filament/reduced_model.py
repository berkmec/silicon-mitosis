"""Exp12 step 4: reduced 4 K floating-body model near threshold, anchored to the official SKY130 BSIM4 card.
Body balance per um of channel width (A/um):
    c_b dV_b/dt = kappa d2V_b/dy2 + fii * M(V_ds) * j_d(V_b) - j_on exp((V_b - V_on)/V_s) - g_sub V_b
    j_d(V_b)    = j_T * exp( (V_ov0 + DV_T(V_b)) / (n phi) )  for V_ov < 0 (subthreshold, 4 K slope SS = n phi ln10),
                  continued smoothly above threshold with the BSIM4 table shape (here capped: strong inversion grows ~linearly,
                  weak feedback, see step 3), V_ov0 = V_gs - V_T(0) - dvt.
BSIM4 anchors (300 K tables, spice/tab_vb*.dat): DV_T(V_b) = V_T(0) - V_T(V_b) (body effect, from the 'vth' output) and the
multiplication M(V_ds) = I_sub / I_d near threshold. 4 K inputs (ASSUMPTIONS, scanned): SS (10..30 mV/dec), V_s, V_on, fii,
operating point (V_ov0), g_sub. Output: bistability map (number of steady states of V_b) and hysteresis windows in V_ds."""
import itertools
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import local_bistability as lb  # noqa: E402  (BSIM4 tables + interpolators)

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "exp12"


def anchors():
    """Body effect DV_T(V_b) and near-threshold multiplication M(V_ds) from the BSIM4 tables (W = 1 um, L = 0.15 um)."""
    VB, VG, VD = lb.VB, lb.VG, lb.VD
    ig, idd = int(np.argmin(abs(VG - 0.75))), int(np.argmin(abs(VD - 1.8)))
    vth = []
    for b in VB:
        d = np.loadtxt(lb.SP / f"tab_vb{b:.1f}.dat")
        vth.append(d[:, 9].reshape(len(VG), len(VD))[ig, idd])
    vth = np.array(vth)
    dvt = vth[0] - vth  # >= 0: forward body bias lowers V_T
    vd = VD[VD >= 0.6]
    m = lb.j_ii(0.75, 0.0, np.zeros(1))  # placeholder to init
    M = np.array([lb.j_ii(0.75, v, np.zeros(1))[0] / max(lb.j_ch(0.75, v, np.zeros(1))[0], 1e-30) for v in vd])
    return dict(VB=VB, dvt=dvt, vth0=float(vth[0]), VD=vd, M=M)


AN = anchors()


def dvt_of(vb):
    return np.interp(np.clip(vb, 0, 1.0), AN["VB"], AN["dvt"])


def M_of(vds):
    return np.exp(np.interp(vds, AN["VD"], np.log(np.maximum(AN["M"], 1e-30))))


def j_in(vb, vds, vov0, ss, fii, j_T=1e-7):
    """Impact-ionisation current into the body (A/um). j_T: drain current per um at V_ov = 0 (BSIM4 at threshold ~ 0.1 uA/um)."""
    nphi = ss / np.log(10)
    vov = vov0 + dvt_of(vb)
    jd = np.where(vov < 0, j_T * np.exp(vov / nphi), j_T * (1 + vov / nphi))  # above threshold: weak (linear) growth
    return fii * M_of(vds) * jd


def j_out(vb, von, vs, g_sub, j_on=1e-6):
    return j_on * np.exp((vb - von) / vs) + g_sub * vb


def n_states(vds, vov0, ss, fii, von, vs, g_sub):
    return len(roots_at(vds, vov0, ss, fii, von, vs, g_sub))


def roots_at(vds, vov0, ss, fii, von, vs, g_sub):
    vb = np.linspace(0, 1.3, 2601)
    F = j_in(vb, vds, vov0, ss, fii) - j_out(vb, von, vs, g_sub)
    s = np.sign(F)
    i = np.where(s[:-1] * s[1:] < 0)[0]
    stable = F[i] > F[i + 1]  # F decreasing through the root: stable (dV_b/dt = F / c_b)
    return [float(vb[k]) for k, st_ in zip(i, stable) if st_]


def _follow(vgrid, start, **kw):
    """Follow a stable branch quasi-statically; return (V at the first jump > 0.15 V or None, V_b along the sweep)."""
    path, jump_at = [], None
    prev = None
    for v in vgrid:
        r = roots_at(v, **kw)
        if not r:
            path.append(np.nan)
            continue
        if prev is None:
            prev = min(r) if start == "low" else max(r)
        nxt = min(r, key=lambda x: abs(x - prev))
        if abs(nxt - prev) > 0.15 and jump_at is None:
            jump_at = float(v)
        prev = nxt
        path.append(prev)
    return jump_at, np.array(path)


def sweep(vds_grid, **kw):
    """Up-sweep from the LOW branch (cold device, body near 0) and down-sweep from the HIGH branch.
    V_up: V_ds where the low branch disappears (kink); V_down: where the high branch disappears (None = stays latched)."""
    V_up, up = _follow(vds_grid, "low", **kw)
    V_down, dn = _follow(vds_grid[::-1], "high", **kw)
    return V_up, V_down, up, dn[::-1]


if __name__ == "__main__":
    print(f"BSIM4 dayanakları: V_T(0) = {AN['vth0']:.3f} V; ΔV_T(V_b = 0,4/0,8/1,0) = "
          f"{dvt_of(0.4):.3f}/{dvt_of(0.8):.3f}/{dvt_of(1.0):.3f} V; M(V_ds = 1,2/1,8/2,4) = "
          f"{M_of(1.2):.2e}/{M_of(1.8):.2e}/{M_of(2.4):.2e}", flush=True)
    grid = dict(ss=(0.010, 0.015, 0.020, 0.030, 0.080), vs=(0.005, 0.010, 0.025), fii=(1, 10, 100),
                vov0=(-0.30, -0.20, -0.10, 0.0), g_sub=(1e-15, 1e-13, 1e-11))  # A/(V um): frozen-well leak, physical range
    von = 1.05
    vds_grid = np.round(np.arange(0.6, 2.601, 0.02), 3)
    rows = []
    for ss, vs, fii, vov0, gs in itertools.product(*grid.values()):
        kw = dict(vov0=vov0, ss=ss, fii=fii, von=von, vs=vs, g_sub=gs)
        V_up, V_down, _, _ = sweep(vds_grid, **kw)
        hyst = V_up is not None and (V_down is None or V_down < V_up - 1e-9)
        rows.append(dict(ss=ss, vs=vs, fii=fii, vov0=vov0, g_sub=gs, V_up=V_up, V_down=V_down,
                         kink_seen=V_up is not None, hysteresis=hyst, latched_to_0p6=hyst and V_down is None,
                         width=(V_up - (V_down if V_down is not None else float(vds_grid[0]))) if hyst else None))
    summ = {}
    for ss in grid["ss"]:
        sel = [r for r in rows if r["ss"] == ss]
        summ[f"{ss*1e3:.0f} mV/dek"] = dict(n=len(sel), kink_up=sum(r["kink_seen"] for r in sel),
                                             hysteresis=sum(r["hysteresis"] for r in sel),
                                             latched=sum(r["latched_to_0p6"] for r in sel),
                                             width_median_V=float(np.median([r["width"] for r in sel if r["hysteresis"]]))
                                             if any(r["hysteresis"] for r in sel) else None)
    print("eşik altı eğim -> {kombinasyon, yukarı taramada kink, histerezis, medyan pencere}:")
    for k, v in summ.items():
        print("  ", k, v)
    ex = [r for r in rows if r["hysteresis"]]
    for r in sorted(ex, key=lambda r: -r["width"])[:5]:
        print("   örnek:", r)
    by_ss = summ
    (OUT / "reduced_bistability.json").write_text(json.dumps(dict(anchors=dict(vth0=AN["vth0"],
        dvt_at=[float(dvt_of(x)) for x in (0.2, 0.4, 0.6, 0.8, 1.0)], M_at={str(v): float(M_of(v)) for v in (1.2, 1.8, 2.4)}),
        von=von, grid=grid, by_ss={str(k): v for k, v in by_ss.items()}, rows=rows, note=__doc__), indent=1), encoding="utf-8")
