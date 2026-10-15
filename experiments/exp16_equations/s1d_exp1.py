"""Exp16 step 1d (ground truth, exp1 FEM data).

(a) C proportional to a: exp1 stores capacitances in units of eps0*a at a = 1 only, and transport.json builds C(L) by
    scaling (C = c * eps0 * eps_r * a). A fitted exponent is therefore 1 by construction -> reported as NOT a test.
(b) Concentric-sphere FEM convergence (validate_sphere.json): known law C(h) = C_exact + A h^p, p ~ 2 for P1 elements.
    Sparse polynomial library {1, h, h^2, h^3} (LOO-CV) and a 3-parameter fit (hold out the finest mesh).
(c) Finite-barrier quantum dot levels E1(L) (schrodinger_wells.csv, Si/SiO2 cube, V0 = 3.1 eV, m_w = 0.26, m_b = 0.5):
    hard wall E1 = 3 pi^2 hbar^2/(2 m_w L^2). Candidate law: E1 = A/(L + L0)^p, with the BenDaniel-Duke penetration
    expectation L0 = 2 m_b/(m_w kappa), kappa = sqrt(2 m_b V0)/hbar (deep-well limit).
Output: results/exp16_equations/s1d_exp1.json, fig_s1d_exp1.png."""
import common as C  # noqa: I001

import csv
import json

import matplotlib
import numpy as np
from scipy.optimize import least_squares

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

C.low_priority()
E1DIR = C.RESULTS / "exp1"
out = {}

# (a) circular check
tr = json.loads((E1DIR / "transport.json").read_text(encoding="utf-8"))
rows = tr["eps"]["3.9"]["table"]
L = np.array([r["L_nm"] for r in rows])
Cf = np.array([r["C_aF"] for r in rows])
out["a_C_vs_size"] = dict(loglog_slope=float(np.polyfit(np.log(L), np.log(Cf), 1)[0]),
                          note="C(L) in transport.json = c_sigma * eps0 * eps_r * a by construction; exp1 FEM itself "
                               "was solved only at a = 1 (units eps0*a). Slope 1 is circular, not a discovery/test.")

# (b) sphere convergence
vs = json.loads((E1DIR / "validate_sphere.json").read_text(encoding="utf-8"))
h = np.array([r["h"] for r in vs["rows"]])
cf = np.array([r["C_fem"] for r in vs["rows"]])
Cex = vs["C_exact"]
Th = np.column_stack([np.ones_like(h), h, h**2, h**3])
names_b = ["1", "h", "h^2", "h^3"]
loo_b = {}
for tau in (1e-4, 1e-3, 1e-2, 3e-2):
    pr = C.loo_cv(Th, cf, tau)
    loo_b[tau] = float(np.sqrt(np.mean((pr - cf) ** 2)))
tau_b = min(loo_b, key=loo_b.get)
cb = C.stlsq_scaled(Th, cf, tau_b)[0]


def model_b(p, hh):
    return p[0] + p[1] * hh ** p[2]


fit_all = least_squares(lambda p: model_b(p, h) - cf, x0=[Cex, 10.0, 2.0])
fit_ho = least_squares(lambda p: model_b(p, h[:-1]) - cf[:-1], x0=[Cex, 10.0, 2.0])
J = fit_all.jac
cov = np.linalg.inv(J.T @ J) * np.sum(fit_all.fun**2) / max(len(h) - 3, 1)
fit_noc = least_squares(lambda p: model_b(p, h[1:]) - cf[1:], x0=[Cex, 10.0, 2.0])  # drop pre-asymptotic h = 0.3
p_known = np.polyfit(np.log(h), np.log(cf / Cex - 1), 1)  # C_exact given: error ~ h^p
out["b_sphere_convergence"] = dict(
    C_exact=Cex, sparse=dict(tau=tau_b, loo_rmse=loo_b, support=[names_b[i] for i in np.flatnonzero(cb)],
                             coef={names_b[i]: float(cb[i]) for i in np.flatnonzero(cb)},
                             C_inf_rel_err=float((cb[0] - Cex) / Cex)),
    power_fit_all=dict(C_inf=fit_all.x[0], A=fit_all.x[1], p=fit_all.x[2], p_std=float(np.sqrt(cov[2, 2])),
                       C_inf_rel_err=float((fit_all.x[0] - Cex) / Cex)),
    power_fit_without_h0p3=dict(C_inf=fit_noc.x[0], A=fit_noc.x[1], p=fit_noc.x[2],
                                C_inf_rel_err=float((fit_noc.x[0] - Cex) / Cex)),
    exponent_with_C_exact_given=dict(p=float(p_known[0]), prefactor=float(np.exp(p_known[1]))),
    power_fit_holdout_finest=dict(C_inf=fit_ho.x[0], p=fit_ho.x[2],
                                  pred_finest=float(model_b(fit_ho.x, h[-1])), true_finest=float(cf[-1]),
                                  rel_err_pred=float((model_b(fit_ho.x, h[-1]) - cf[-1]) / cf[-1])))
print("(b)", out["b_sphere_convergence"], flush=True)

# (c) E1(L)
with open(E1DIR / "schrodinger_wells.csv", encoding="utf-8") as f:
    rr = list(csv.DictReader(f))
Lw = np.array([float(r["L_nm"]) for r in rr])
E1 = np.array([float(r["E1_eV"]) for r in rr])
E2 = np.array([float(r["E2_eV"]) for r in rr])
E1hw = np.array([float(r["E1_hardwall_eV"]) for r in rr])
A_hw = float(np.median(E1hw * Lw**2))
hb2_2m0 = 0.0380998  # eV nm^2
mw, mb, V0 = 0.26, 0.5, 3.1
kappa = np.sqrt(mb * V0 / hb2_2m0)
L0_bdd = 2 * mb / (mw * kappa)


def model_c(p, LL):
    return p[0] / (LL + p[1]) ** p[2]


fc = least_squares(lambda p: model_c(p, Lw) - E1, x0=[A_hw, 0.5, 2.0], x_scale=[1, 0.1, 0.1])
Jc = fc.jac
covc = np.linalg.inv(Jc.T @ Jc) * np.sum(fc.fun**2) / max(len(Lw) - 3, 1)
loo = []
for i in range(len(Lw)):
    m = np.arange(len(Lw)) != i
    fi = least_squares(lambda p: model_c(p, Lw[m]) - E1[m], x0=fc.x)
    loo.append(float((model_c(fi.x, Lw[i]) - E1[i]) / E1[i]))
# p fixed to 2 and A fixed to the hard-wall value: one free parameter L0 per size (is it constant?)
L0_each = Lw * (np.sqrt(E1hw / E1) - 1)
# one-parameter law E1 = A_hw/(L + L0)^2 (physics-constrained), LOO; and the zero-parameter BenDaniel-Duke prediction
f1 = least_squares(lambda p: A_hw / (Lw + p[0]) ** 2 / E1 - 1, x0=[0.5])
loo1 = []
for i in range(len(Lw)):
    m = np.arange(len(Lw)) != i
    fi = least_squares(lambda p: A_hw / (Lw[m] + p[0]) ** 2 / E1[m] - 1, x0=[0.5])
    loo1.append(float(A_hw / (Lw[i] + fi.x[0]) ** 2 / E1[i] - 1))
bdd_err = (A_hw / (Lw + L0_bdd) ** 2 / E1 - 1).tolist()
# sparse library in 1/L for E1 (shows the collinearity problem honestly)
Thc = np.column_stack([Lw**-1.0, Lw**-2.0, Lw**-3.0, Lw**-4.0])
names_c = ["1/L", "1/L^2", "1/L^3", "1/L^4"]
loo_c = {tau: float(np.sqrt(np.mean((C.loo_cv(Thc, E1, tau) - E1) ** 2))) for tau in (1e-4, 1e-3, 1e-2, 5e-2)}
tau_c = min(loo_c, key=loo_c.get)
cc = C.stlsq_scaled(Thc, E1, tau_c)[0]
fc2 = least_squares(lambda p: model_c(p, Lw) - (E2 - E1), x0=[A_hw, 0.5, 2.0], x_scale=[1, 0.1, 0.1])
out["c_levels"] = dict(
    A_hardwall_eVnm2=A_hw, kappa_per_nm=kappa, L0_BenDanielDuke_nm=L0_bdd,
    fit_E1=dict(A=fc.x[0], L0=fc.x[1], p=fc.x[2], std=np.sqrt(np.diag(covc)).tolist(),
                A_rel_to_hardwall=float(fc.x[0] / A_hw), loo_rel_err=loo, max_abs_loo=float(np.max(np.abs(loo)))),
    L0_from_each_size_with_A_hw_p2=dict(zip(map(str, Lw), L0_each.tolist())),
    one_param_fit=dict(L0=float(f1.x[0]), loo_rel_err=loo1, max_abs_loo=float(np.max(np.abs(loo1)))),
    zero_param_BDD_rel_err=dict(zip(map(str, Lw), bdd_err)),
    fit_dE=dict(A=fc2.x[0], L0=fc2.x[1], p=fc2.x[2]),
    sparse_1overL=dict(tau=tau_c, loo_rmse=loo_c, support=[names_c[i] for i in np.flatnonzero(cc)],
                       coef={names_c[i]: float(cc[i]) for i in np.flatnonzero(cc)}))
print("(c)", out["c_levels"], flush=True)
C.save_json("s1d_exp1.json", out)

fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
hh = np.linspace(0.0, 0.32, 100)
axs[0].plot(h, cf, "ko", label="FEM (validate_sphere)")
axs[0].plot(hh, model_b(fit_all.x, hh), "r-", label=f"C∞ + A h^p, p={fit_all.x[2]:.2f}")
axs[0].axhline(Cex, color="b", ls=":", label="analitik")
axs[0].set_xlabel("h")
axs[0].set_ylabel("C")
axs[0].legend(fontsize=8)
LL = np.linspace(0.8, 9, 200)
axs[1].loglog(Lw, E1, "ko", label="FEM E1")
axs[1].loglog(LL, A_hw / LL**2, "b:", label="sert duvar A/L²")
axs[1].loglog(LL, model_c(fc.x, LL), "r-", label=f"serbest A/(L+L0)^p: L0={fc.x[1]:.3f}, p={fc.x[2]:.2f}")
axs[1].loglog(LL, A_hw / (LL + L0_bdd) ** 2, "g--", label="BenDaniel–Duke A_sd/(L+0,603)² (parametresiz)")
axs[1].set_xlabel("L (nm)")
axs[1].set_ylabel("E1 (eV)")
axs[1].set_title(f"BenDaniel–Duke beklentisi L0 = {L0_bdd:.3f} nm", fontsize=9)
axs[1].legend(fontsize=8)
fig.tight_layout()
fig.savefig(C.OUT / "fig_s1d_exp1.png", dpi=110)
