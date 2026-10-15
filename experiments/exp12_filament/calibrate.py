"""Exp12 step 1: calibrate the compact model (src/qcell/fet4k.py) to the official SKY130 BSIM4 card in ngspice at 300 K.
Data: spice/iv_vb{0,0.2,0.4,0.6}.dat from spice/iv_grid.cir (L = 0.15 um, W = 1 um, Vgs 0.6..1.8, Vds 0..2.6).
Fit 1: vt0, mu, lam, nss to Id (relative error, Id > 1 uA/um).  Fit 2: effective beta0 and an alpha scale to Isub
(Isub > 1e-12 A/um), because my V_dsat is not BSIM4's.  Output: results/exp12/calibration.json + fig_calibration.png."""
import json
import sys
from pathlib import Path

import matplotlib
import numpy as np
from scipy.optimize import least_squares

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from qcell.fet4k import FetParams, id_per_w, iii_per_w  # noqa: E402

SP = Path(__file__).resolve().parent / "spice"
OUT = ROOT / "results" / "exp12"
OUT.mkdir(parents=True, exist_ok=True)
W = 1e-6  # m (data are for W = 1 um)

rows = []
for vb in ("0", "0.2", "0.4", "0.6"):
    d = np.loadtxt(SP / f"iv_vb{vb}.dat")
    nvd = 53
    vgs = np.repeat(np.round(np.arange(0.6, 1.8001, 0.1), 3), nvd)
    rows.append(np.column_stack([vgs, d[:, 0], np.full(len(d), float(vb)), d[:, 1] / W, d[:, 3] / W]))
data = np.vstack(rows)  # vgs, vds, vbs, Id [A/m], Isub [A/m]
vgs, vds, vbs, idd, isub = data.T

sel = idd > 1.0  # > 1 uA/um
p0 = FetParams()


def res_id(x):
    p = p0.at(vt0=x[0], mu=x[1], lam=x[2], nss=x[3])
    return np.log(id_per_w(p, vgs[sel], vds[sel], vbs[sel]) / idd[sel])


f1 = least_squares(res_id, x0=[0.6, 0.03, 0.05, 1.5], bounds=([0.2, 0.002, 0.0, 0.8], [1.0, 0.1, 1.0, 5.0]))
p1 = p0.at(vt0=f1.x[0], mu=f1.x[1], lam=f1.x[2], nss=f1.x[3])
e_id = np.exp(np.abs(res_id(f1.x))) - 1

sel2 = (isub > 1e-6) & (vgs >= 0.8)  # > 1 pA/um
raw = iii_per_w(p1, vgs[sel2], vds[sel2], vbs[sel2])
e_raw = np.log(raw / isub[sel2])


def res_ii(x):
    p = p1.at(beta0=x[0], alpha0=p1.alpha0 * x[1], alpha1=p1.alpha1 * x[1])
    return np.log(iii_per_w(p, vgs[sel2], vds[sel2], vbs[sel2]) / isub[sel2])


f2 = least_squares(res_ii, x0=[13.86, 1.0], bounds=([3.0, 0.01], [40.0, 100.0]))
p2 = p1.at(beta0=f2.x[0], alpha0=p1.alpha0 * f2.x[1], alpha1=p1.alpha1 * f2.x[1])
e_ii = np.exp(np.abs(res_ii(f2.x))) - 1

out = dict(
    note="ngspice 44.2 + official SKY130 nfet_01v8 tt card (BSIM4 4.5), 300 K, L=0.15 um, W=1 um",
    fit_id=dict(vt0=f1.x[0], mu=f1.x[1], lam=f1.x[2], nss=f1.x[3], n_points=int(sel.sum()),
                rel_err_median=float(np.median(e_id)), rel_err_p90=float(np.percentile(e_id, 90)),
                rel_err_max=float(e_id.max())),
    isub_with_official_alpha_beta=dict(n_points=int(sel2.sum()), log_ratio_median=float(np.median(e_raw)),
                                       log_ratio_p10=float(np.percentile(e_raw, 10)),
                                       log_ratio_p90=float(np.percentile(e_raw, 90))),
    fit_isub=dict(beta0_eff=f2.x[0], alpha_scale=f2.x[1], rel_err_median=float(np.median(e_ii)),
                  rel_err_p90=float(np.percentile(e_ii, 90)), rel_err_max=float(e_ii.max())),
    params=dict(p2.__dict__),
)
(OUT / "calibration.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
print(json.dumps({k: v for k, v in out.items() if k != "params"}, indent=1))

fig, axs = plt.subplots(1, 2, figsize=(12, 4.6))
for vb, ls in ((0.0, "-"), (0.4, "--")):
    for vg in (0.9, 1.2, 1.8):
        m = (np.isclose(vbs, vb)) & (np.isclose(vgs, vg))
        axs[0].plot(vds[m], idd[m] * 1e-6 * 1e6, "k" + ls, lw=1)
        axs[0].plot(vds[m], id_per_w(p2, vgs[m], vds[m], vbs[m]) * 1e-6 * 1e6, "r" + ls, lw=1)
        axs[1].semilogy(vds[m], np.maximum(isub[m] * 1e-6, 1e-16), "k" + ls, lw=1)
        axs[1].semilogy(vds[m], np.maximum(iii_per_w(p2, vgs[m], vds[m], vbs[m]) * 1e-6, 1e-16), "r" + ls, lw=1)
axs[0].set_xlabel("V_ds (V)")
axs[0].set_ylabel("I_d (µA/µm)")
axs[0].set_title("siyah: ngspice+SKY130 BSIM4, kırmızı: kompakt model\nV_gs = 0,9 / 1,2 / 1,8 V; düz V_bs=0, kesik V_bs=0,4 V", fontsize=9)
axs[1].set_ylim(1e-14, 1e-4)
axs[1].set_xlabel("V_ds (V)")
axs[1].set_ylabel("I_sub (A/µm)")
axs[1].set_title(f"darbe iyonlaşması: β₀,eff = {f2.x[0]:.2f} V (resmî 13,86), α ölçeği {f2.x[1]:.2f}", fontsize=9)
fig.tight_layout()
fig.savefig(OUT / "fig_calibration.png", dpi=110)
