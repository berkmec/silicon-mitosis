"""Schroedinger FEM: validation (hard wall, separable finite well) and the Si/SiO2 cube well."""
import csv
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from qcell import schro  # noqa: E402

OUT = ROOT / "results" / "exp1"
OUT.mkdir(parents=True, exist_ok=True)
M_IN, M_OUT, V0 = 0.26, 0.5, 3.1
KB_EV = 8.617333262e-5
res = {"validation": [], "wells": []}

# --- V1: hard-wall cube, P2 FEM vs analytic 3,6,6,6,9 ---------------------
L = 3.0
for h in (1.0, 0.6, 0.4, 0.3):
    gm = schro.build_well(L, 0.0, h, h, barrier=False)
    ev, nfree = schro.solve_eigs(gm, np.zeros(gm.mesh.t.shape[1]), np.ones(gm.mesh.t.shape[1]) / M_IN, k=5)
    ex = np.array([3, 6, 6, 6, 9]) * schro.hardwall(L, M_IN, (1, 0, 0)) if False else \
        np.array([3, 6, 6, 6, 9]) * schro.HB * np.pi ** 2 / (M_IN * L ** 2)
    err = np.abs(ev - ex) / ex
    res["validation"].append(dict(test="hardwall", h=h, dofs=nfree, max_rel_err=float(err.max()),
                                  E_fem=ev.tolist(), E_exact=ex.tolist()))
    print(f"V1 hard-wall L={L} h={h}: dofs={nfree} max rel err={err.max():.2e}", flush=True)

# --- V2: separable finite well V(x)+V(y)+V(z), equal masses: exact 3*E1D --
L = 3.0
margin = 0.8
e1 = schro.e1d_finite(L, V0, M_IN)
sep = lambda x: V0 * ((np.abs(x[0]) > L / 2) * 1.0 + (np.abs(x[1]) > L / 2) + (np.abs(x[2]) > L / 2))
for h_int in (0.4, 0.3, 0.2):
    gm = schro.build_well(L, margin, h_int, 1.0, barrier=True)
    nt = gm.mesh.t.shape[1]
    ev, nfree = schro.solve_eigs(gm, np.zeros(nt), np.ones(nt) / M_IN, k=3, potential_fn=sep)
    err = abs(ev[0] - 3 * e1) / (3 * e1)
    res["validation"].append(dict(test="separable_finite", h_int=h_int, dofs=nfree, E_fem=float(ev[0]),
                                  E_exact=float(3 * e1), rel_err=float(err)))
    print(f"V2 separable finite well L={L} h_int={h_int}: dofs={nfree} E={ev[0]:.5f} exact={3*e1:.5f} err={err:.2e}", flush=True)

# --- V3: mesh convergence of the real (Si/SiO2, BenDaniel-Duke) well -------
L = 2.0
prev = None
for h_int in (0.35, 0.25, 0.18, 0.13):
    gm = schro.build_well(L, margin, h_int, max(L / 5, 0.5), barrier=True)
    is_well = gm.region == gm.well_vol
    ev, nfree = schro.solve_eigs(gm, np.where(is_well, 0.0, V0), np.where(is_well, 1 / M_IN, 1 / M_OUT), k=4)
    d = None if prev is None else abs(ev[0] - prev) / prev
    prev = ev[0]
    res["validation"].append(dict(test="real_well_convergence", L=L, h_int=h_int, dofs=nfree, E1=float(ev[0]),
                                  rel_change_vs_prev=d))
    print(f"V3 real well L={L} h_int={h_int}: dofs={nfree} E1={ev[0]:.5f} rel change={d}", flush=True)

# --- main: Si/SiO2 cube well, BenDaniel-Duke -------------------------------
rows = []
for L in (1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0):
    h_int = 0.25 if L <= 4 else 0.3
    gm = schro.build_well(L, margin, h_int, max(L / 5, 0.5), barrier=True)
    nt = gm.mesh.t.shape[1]
    is_well = gm.region == gm.well_vol
    V = np.where(is_well, 0.0, V0)
    invm = np.where(is_well, 1 / M_IN, 1 / M_OUT)
    ev, nfree = schro.solve_eigs(gm, V, invm, k=5)
    e_hw = schro.HB * np.pi ** 2 / (M_IN * L ** 2)
    rows.append(dict(L_nm=L, dofs=nfree, E1_eV=ev[0], E2_eV=ev[1], dE_eV=ev[1] - ev[0],
                     dE_over_kB_K=(ev[1] - ev[0]) / KB_EV,
                     E1_hardwall_eV=3 * e_hw, dE_hardwall_eV=3 * e_hw,
                     ratio_dE_fem_over_hw=(ev[1] - ev[0]) / (3 * e_hw),
                     deg=float(ev[3] - ev[1])))
    print({k: (round(float(v), 5) if isinstance(v, (float, np.floating)) else v) for k, v in rows[-1].items()}, flush=True)
    res["wells"] = [{k: (float(v) if isinstance(v, np.floating) else v) for k, v in r.items()} for r in rows]
    (OUT / "schrodinger.json").write_text(json.dumps(res, indent=2))
res["wells"] = [{k: (float(v) if isinstance(v, np.floating) else v) for k, v in r.items()} for r in rows]
(OUT / "schrodinger.json").write_text(json.dumps(res, indent=2))
with open(OUT / "schrodinger_wells.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)
