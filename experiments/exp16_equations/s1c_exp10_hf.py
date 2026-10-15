"""Exp16 step 1c (ground truth): 2D electron gas inverse compressibility K(n)/K_kin from exp10 (phi_physics.json).

Variant B (Si, g = 4 valleys x spin, Hartree-Fock, NO gate screening, eps_bar) has the closed form
    K/K_kin = 1 - (sqrt(g)/pi) r_s = 1 - c n^(-1/2),   r_s = 1/(sqrt(pi n) a_B*)   (kinetic const + exchange ~ -n^(-1/2)).
Library: powers z^p, z = n / 1e12 cm^-2, p in {-2,-1.5,-1,-0.5,0,0.5,1,1.5,2}; selection on an EXTRAPOLATION split
(train n >= 1e10, validate n < 1e10). Bootstrap + threshold path. The other variants (A: DMC correlation, C: gate-screened
exchange, D/E: correlation proxies) are run through the same pipeline as a contrast (no closed form expected).
Output: results/exp16_equations/s1c_exp10.json, fig_s1c_exp10.png."""
import common as C  # noqa: I001

import json

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

C.low_priority()
d = json.loads((C.RESULTS / "exp10" / "phi_physics.json").read_text(encoding="utf-8"))
P = np.array([-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2])
NAMES = [("1" if p == 0 else f"z^{p:g}") for p in P]
TAUS = np.logspace(-7, -0.5, 30)
out = dict(names=NAMES, cases={})
fig, axs = plt.subplots(1, 2, figsize=(12, 4.3))
for name, cs in d["cases"].items():
    n = np.array(cs["n_cm2"])
    y = np.array(cs["K_over_Kkin"])
    rs = np.array(cs["rs"])
    z = n / 1e12
    Th = z[:, None] ** P[None, :]
    tr, va = n >= 1e10 - 1, n < 1e10 - 1
    path = C.threshold_path(Th[tr], y[tr], Th[va], y[va], TAUS)
    sel = C.select_model(path)
    coef = sel["coef"][0]
    bs = C.bootstrap(Th, y, sel["tau"], n_boot=300)
    # exponent scan without a library: 1 - K/Kkin = c z^p (log-log, all points)
    lp = np.polyfit(np.log(z), np.log(np.abs(1 - y)), 1)
    alpha_rs = float(np.median((1 - y) / rs))  # (1 - K/Kkin)/r_s ; constant only for pure HF without screening
    res = dict(n_points=len(n), n_range=[float(n.min()), float(n.max())], tau=sel["tau"],
               support=[NAMES[i] for i in np.flatnonzero(coef)],
               coef={NAMES[i]: float(coef[i]) for i in np.flatnonzero(coef)},
               r2_val_extrap=sel["r2_val"][0], rel_rmse_val=sel["rmse_val"][0],
               bootstrap_p_incl={NAMES[i]: float(bs["p_incl"][0, i]) for i in range(len(P))},
               bootstrap_mean={NAMES[i]: float(bs["mean"][0, i]) for i in range(len(P))},
               bootstrap_std={NAMES[i]: float(bs["std"][0, i]) for i in range(len(P))},
               loglog_exponent_of_1_minus_K=float(lp[0]), loglog_prefactor=float(np.exp(lp[1])),
               median_one_minus_K_over_rs=alpha_rs,
               spread_one_minus_K_over_rs=float(np.ptp((1 - y) / rs) / abs(alpha_rs)),
               path=[dict(tau=r["tau"], n_terms=r["n_terms"][0], r2_val=r["r2_val"][0]) for r in path])
    out["cases"][name] = res
    print(name, "| support", res["support"], "| coef", {k: round(v, 5) for k, v in res["coef"].items()},
          "| R2val", round(res["r2_val_extrap"], 6), "| loglog p", round(res["loglog_exponent_of_1_minus_K"], 4),
          "| (1-K)/rs", round(alpha_rs, 5), flush=True)
    axs[0].loglog(n, np.abs(1 - y), "o", ms=3, label=name)
    axs[1].semilogx(n, y - Th @ coef, "-", lw=1, label=name)
# analytic reference for B
cB = [k for k in d["cases"] if k.startswith("B_")][0]
rs_1e12 = float(d["operating_point_1e12"]["rs"])
g = 4
out["analytic_B"] = dict(alpha=np.sqrt(g) / np.pi, coef_z_m05=-np.sqrt(g) / np.pi * rs_1e12, const=1.0,
                         formula="K/K_kin = 1 - (sqrt(g)/pi) r_s = 1 - (2/pi) r_s(1e12) z^-1/2, g = 4")
fb = out["cases"][cB]["coef"]
out["analytic_B"]["rel_err_coef_z_m05"] = (fb.get("z^-0.5", np.nan) - out["analytic_B"]["coef_z_m05"]) / abs(
    out["analytic_B"]["coef_z_m05"])
out["analytic_B"]["rel_err_const"] = fb.get("1", np.nan) - 1.0
print("analytic B:", out["analytic_B"], flush=True)
C.save_json("s1c_exp10.json", out)
axs[0].set_xlabel("n (cm⁻²)")
axs[0].set_ylabel("|1 − K/K_kin|")
axs[0].set_title("B (ekransız HF) tam n^(−1/2); diğerleri sapıyor", fontsize=9)
axs[0].legend(fontsize=6)
axs[1].set_xlabel("n (cm⁻²)")
axs[1].set_ylabel("artık (veri − seyrek model)")
axs[1].set_title("seyrek kuvvet yasası modelinin artıkları", fontsize=9)
axs[1].legend(fontsize=6)
fig.tight_layout()
fig.savefig(C.OUT / "fig_s1c_exp10.png", dpi=110)
