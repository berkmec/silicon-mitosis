"""Exp12 step 2: single SKY130 nFET (L = 0.15 um, W = 10 um) with a floating (frozen) body -> kink and hysteresis?
DC sweep of V_ds up then down (ngspice keeps the previous solution as initial guess -> hysteresis if bistable).
Knobs (4 K emulation, see spice_tools.py): fii (impact ionisation x), dvt (V_T shift), von (body diode turn-on), rsub."""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import spice_tools as st  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "exp12"
OUT.mkdir(parents=True, exist_ok=True)
W = 10.0  # um


def netlist(name, vg=1.2, fii=1.0, dvt=0.0, von=1.0, rsub=1e12, n=1.0, fjun=1e-30):
    dm = st.diode_model(von, n=n, i_on=1e-6, w_um=W)
    dname = dm.split()[1]
    return (f"* single floating-body nFET, {name}\n" + st.header(fii, fjun) + dm +
            f"Vd d 0 0\nVg g 0 {vg - dvt}\nVb bs 0 0\n"
            f"X1 d g 0 b sky130_fd_pr__nfet_01v8 l=0.15 w={W} as=0 ad=0 ps=0 pd=0\n"
            f"Rsub b bs {rsub}\nDb b 0 {dname}\n"
            ".control\nsave all @m.x1.msky130_fd_pr__nfet_01v8[isub]\n"
            "dc Vd 0 2.6 0.01\n"
            f"wrdata {st.wsl(st.SP / (name + '_up.dat'))} -i(Vd) v(b) @m.x1.msky130_fd_pr__nfet_01v8[isub]\n"
            "dc Vd 2.6 0 -0.01\n"
            f"wrdata {st.wsl(st.SP / (name + '_dn.dat'))} -i(Vd) v(b) @m.x1.msky130_fd_pr__nfet_01v8[isub]\n"
            "quit\n.endc\n.end\n")


cases = {
    "warm_300K_contacted": dict(fii=1, dvt=0.0, von=0.6, rsub=1e3, fjun=1.0),
    "frozen_bsimdiode_on": dict(fii=1, dvt=0.15, von=1.0, rsub=1e12, fjun=1.0),
    "frozen_von1.0_fii1": dict(fii=1, dvt=0.15, von=1.0, rsub=1e12),
    "frozen_von1.0_fii10": dict(fii=10, dvt=0.15, von=1.0, rsub=1e12),
    "frozen_von0.8_fii1": dict(fii=1, dvt=0.15, von=0.8, rsub=1e12),
    "frozen_von1.0_fii1_n0.3": dict(fii=1, dvt=0.15, von=1.0, rsub=1e12, n=0.3),
    "frozen_von1.0_fii1_rsub1e9": dict(fii=1, dvt=0.15, von=1.0, rsub=1e9),
}
res = {}
for name, kw in cases.items():
    log, bad = st.run(name, netlist(name, **kw))
    up, dn = np.loadtxt(st.SP / f"{name}_up.dat"), np.loadtxt(st.SP / f"{name}_dn.dat")
    vd_u, id_u, vb_u = up[:, 0], up[:, 1], up[:, 3]
    vd_d, id_d, vb_d = dn[:, 0][::-1], dn[:, 1][::-1], dn[:, 3][::-1]
    dI = np.abs(id_u - id_d) / np.maximum(np.abs(id_u), 1e-12)
    hyst = vd_u[dI > 0.05]
    jump = np.argmax(np.diff(vb_u)) if len(vb_u) > 1 else 0
    res[name] = dict(knobs=kw, errors=bad[:3], vb_max_up=float(vb_u.max()),
                     vd_at_largest_vb_jump_up=float(vd_u[jump + 1]), largest_vb_jump_up=float(np.diff(vb_u).max()),
                     hysteresis_vd_range=[float(hyst.min()), float(hyst.max())] if len(hyst) else None,
                     id_at_1p8_up_uA=float(np.interp(1.8, vd_u, id_u) * 1e6),
                     id_at_1p8_dn_uA=float(np.interp(1.8, vd_d, id_d) * 1e6))
    print(name, json.dumps(res[name], ensure_ascii=False), flush=True)
(OUT / "single_kink.json").write_text(json.dumps(res, indent=1, ensure_ascii=False), encoding="utf-8")
