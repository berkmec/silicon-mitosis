"""Exp16 step 2e-2: droplet growth law with the far-field supersaturation as a second state variable (exp11 A2_charge).
dR/dt depends on c_far(t) (it rises 0.011 -> 0.078 and decays to 0.033), so dR/dt = f(R) alone is ill-posed (s2e).
Expectation (2D diffusion-limited growth + Gibbs-Thomson, Delta c ~ 1): dR/dt ~ a c_far / R - b / R^2 with
b/a = gamma = 1/12 (c_eq(R) = gamma / R for f = c^2 (1-c)^2 / 2, kappa = 1/4); a carries a slowly varying log factor.
Pooled over the 8 runs, library {1, R, 1/R, 1/R^2, c, c/R, c R} (c = c_far), leave-one-run-out validation (rel-RMSE).
Output: results/exp16_equations/s2e2_growth.json"""
import common as C  # noqa: I001

import glob
import json
from pathlib import Path

import numpy as np

C.low_priority()
NAMES = ["1", "R", "1/R", "1/R^2", "c", "c/R", "c·R"]


def lib(R, c):
    return np.column_stack([np.ones_like(R), R, 1 / R, 1 / R**2, c, c / R, c * R])


RR, CC, DD, GG = [], [], [], []
for gi, f in enumerate(sorted(glob.glob(str(C.RESULTS / "exp11" / "runs" / "A2_charge_*.json")))):
    h = json.loads(Path(f).read_text(encoding="utf-8"))["hist"]
    if any(x["n"] != 1 for x in h):
        continue
    t = np.array([x["t"] for x in h])
    R = np.array([x["rmean"] for x in h])
    c = np.array([x["cfar"] for x in h])
    d = np.gradient(R, t)
    RR.append(R[1:-1])
    CC.append(c[1:-1])
    DD.append(d[1:-1])
    GG.append(np.full(len(R) - 2, gi))
R, c, D, G = map(np.concatenate, (RR, CC, DD, GG))
Th = lib(R, c)
path = []
for tau in np.logspace(-4, -0.3, 16):
    pred = np.zeros_like(D)
    for g in np.unique(G):
        m = G != g
        pred[~m] = Th[~m] @ C.stlsq_scaled(Th[m], D[m], tau)[0]
    cf = C.stlsq_scaled(Th, D, tau)[0]
    path.append(dict(tau=float(tau), lorocv_rel_rmse=C.rel_rmse(D, pred), lorocv_r2=C.r2(D, pred),
                     support=[NAMES[i] for i in np.flatnonzero(cf)], coef={NAMES[i]: float(cf[i])
                                                                          for i in np.flatnonzero(cf)}))
best = min(p["lorocv_rel_rmse"] for p in path)
sel = min((p for p in path if p["lorocv_rel_rmse"] <= 1.05 * best), key=lambda p: (len(p["support"]), -p["tau"]))
# the 2-term physical model forced (OLS) and its leave-one-run-out error
idx = [NAMES.index("c/R"), NAMES.index("1/R^2")]
pred2 = np.zeros_like(D)
for g in np.unique(G):
    m = G != g
    cg, *_ = np.linalg.lstsq(Th[m][:, idx], D[m], rcond=None)
    pred2[~m] = Th[~m][:, idx] @ cg
c2, *_ = np.linalg.lstsq(Th[:, idx], D, rcond=None)
out = dict(n_points=len(D), n_runs=len(np.unique(G)), path=path, selected=sel,
           forced_two_term=dict(a_c_over_R=float(c2[0]), b_inv_R2=float(-c2[1]), b_over_a=float(-c2[1] / c2[0]),
                                gamma_expected=1 / 12, lorocv_rel_rmse=C.rel_rmse(D, pred2), lorocv_r2=C.r2(D, pred2)))
C.save_json("s2e2_growth.json", out)
print("selected:", sel)
print("forced:", out["forced_two_term"])
