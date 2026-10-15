"""Exp16 step 2e (conservation-law check on existing data): exp11 'A2_charge' runs, one round droplet growing under
uniform charge injection J (s = J, no leakage) in an L x L periodic box, with weak gate-screened Coulomb repulsion A.
Mass balance: d/dt [pi R^2 (c_in - c_far) + c_far L^2] = J L^2  ->  dR/dt ~ J L^2 / (2 pi dc R) if c_far stays fixed.
Sparse regression of dR/dt (central differences, Delta t = 250) on {1, R, 1/R, R^2, 1/R^2}; trained on the first 60 %
of each run, validated on the last 40 % (extrapolation to larger R). Also the integral form R^2(t) = R0^2 + s t.
Output: results/exp16_equations/s2e_charging.json, fig_s2e_charging.png."""
import common as C  # noqa: I001

import glob
import json
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

C.low_priority()
NAMES = ["1", "R", "1/R", "R^2", "1/R^2"]


def lib(R):
    return np.column_stack([np.ones_like(R), R, 1 / R, R**2, 1 / R**2])


TAUS = np.logspace(-5, -0.5, 20)
out = {}
fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
for f in sorted(glob.glob(str(C.RESULTS / "exp11" / "runs" / "A2_charge_*.json"))):
    d = json.loads(Path(f).read_text(encoding="utf-8"))
    sp, h = d["spec"], d["hist"]
    n = np.array([x["n"] for x in h])
    if np.any(n != 1):
        out[sp["name"]] = dict(skipped="n != 1 at some time")
        continue
    t = np.array([x["t"] for x in h])
    R = np.array([x["rmean"] for x in h])
    asp = np.array([x.get("aspect", 1.0) for x in h])
    cfar = np.array([x["cfar"] for x in h])
    J, L = sp["J"], sp["L"]
    dRdt = np.gradient(R, t)
    k = slice(1, -1)
    Rk, Dk = R[k], dRdt[k]
    ntr = int(0.6 * len(Rk))
    path = C.threshold_path(lib(Rk[:ntr]), Dk[:ntr], lib(Rk[ntr:]), Dk[ntr:], TAUS)
    sel = C.select_model_r2(path, r2_tol=1e-3)
    cf = sel["coef"][0]
    s_fit = np.polyfit(t, R**2, 1)
    dc = 1.0 - cfar
    pred_coef = J * L**2 / (2 * np.pi * np.median(dc))
    # integral variable z = R^2: conservation predicts dz/dt = const = J L^2/(pi dc). A near-constant target makes R^2
    # meaningless, so the selection uses validation rel-RMSE (rule A: sparsest within 5 % of the best).
    z = R**2
    dz = np.gradient(z, t)[k]
    zk = z[k]
    libz = lambda zz: np.column_stack([np.ones_like(zz), np.sqrt(zz), zz, 1 / np.sqrt(zz), zz**2])  # noqa: E731
    NZ = ["1", "sqrt z", "z", "1/sqrt z", "z^2"]
    pz = C.threshold_path(libz(zk[:ntr]), dz[:ntr], libz(zk[ntr:]), dz[ntr:], TAUS)
    sz = C.select_model(pz, tol=0.05)
    cz = sz["coef"][0]
    zres = dict(support=[NZ[i] for i in np.flatnonzero(cz)], coef={NZ[i]: float(cz[i]) for i in np.flatnonzero(cz)},
                val_rel_rmse=sz["rmse_val"][0], tau=sz["tau"],
                tau_range_const_only=[min(tt), max(tt)] if (tt := [r["tau"] for r in pz
                                                                    if list(np.flatnonzero(r["coef"][0])) == [0]])
                else None)
    out[sp["name"]] = dict(A=sp["A"], d=sp["d"], J=J, L=L, aspect_max=float(np.nanmax(asp)),
                           R_range=[float(R.min()), float(R.max())], cfar_range=[float(cfar.min()), float(cfar.max())],
                           support=[NAMES[i] for i in np.flatnonzero(cf)],
                           coef={NAMES[i]: float(cf[i]) for i in np.flatnonzero(cf)}, r2_val_extrap=sel["r2_val"][0],
                           tau=sel["tau"], predicted_1_over_R_coef=float(pred_coef),
                           ratio_found_over_pred=float(cf[2] / pred_coef) if cf[2] != 0 else None,
                           R2_vs_t_slope=float(s_fit[0]), R2_slope_pred=float(J * L**2 / (np.pi * np.median(dc))),
                           R2_linear_r2=C.r2(R**2, np.polyval(s_fit, t)), dR_val_rel_rmse=sel["rmse_val"][0],
                           integral_form_z=zres)
    print("   z=R^2:", zres, flush=True)
    print(sp["name"], out[sp["name"]]["support"], {k_: round(v, 5) for k_, v in out[sp["name"]]["coef"].items()},
          "pred 1/R", round(pred_coef, 5), "R2val", round(sel["r2_val"][0], 4),
          "R^2 slope", round(s_fit[0], 5), "pred", round(J * L**2 / (np.pi * np.median(dc)), 5), flush=True)
    axs[0].plot(t, R**2, ".", ms=2, label=sp["name"])
    axs[1].plot(Rk, Dk, ".", ms=3)
    o = np.argsort(Rk)
    axs[1].plot(Rk[o], (lib(Rk) @ cf)[o], "k-", lw=0.6)
C.save_json("s2e_charging.json", out)
axs[0].set_xlabel("t")
axs[0].set_ylabel("R²")
axs[0].set_title("R² doğrusal büyüme (kütle korunumu)", fontsize=9)
axs[0].legend(fontsize=5)
axs[1].set_xlabel("R")
axs[1].set_ylabel("dR/dt")
axs[1].set_title("noktalar: veri; çizgi: seçilen seyrek model", fontsize=9)
fig.tight_layout()
fig.savefig(C.OUT / "fig_s2e_charging.png", dpi=110)
