"""Exp16 step 2a: equations for the exp11 2D active droplets (Zwicker-type chemically active emulsion, qcell.active).

Part A (existing data, results/exp11/runs/B_e*.json, Delta t = 250): stationary radius R*(eps, nu_-) of single round
droplets; sparse regression on a flux-balance-motivated library, leave-one-out CV; dissolution test; relaxation rate.
Flux-balance expectation (2D, one droplet, c_in ~ 1): d(pi R^2)/dt = 2 pi R j(R, eps) - nu_- pi R^2  ->
dR/dt = j - nu_- R / 2, so R* = 2 j / nu_-  with nu_- = dnu + k (1 - eps), k = 0.01.
Part B (NEW dense runs, same solver/units, flagged as new data): R(t) every 5 t0 for 6 (eps, dnu) x 2 initial radii,
L = 64 w, n = 128 (dx = 0.5 as in exp11), dt = 0.2, float32. dR/dt from PySINDy SmoothedFiniteDifference; sparse
regression of dR/dt on {1, R, 1/R, R^2, 1/R^2} x {1, eps, nu_-}; leave-one-condition-out validation; test of the sharp
flux-balance prediction coefficient(nu_- R) = -1/(2 dc).
Output: results/exp16_equations/s2a_droplets.json, fig_s2a_droplets.png."""
import common as C  # noqa: I001

import glob
import json
import time

import matplotlib
import numpy as np
import pysindy as ps

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from qcell.active import ActiveCH, rates_zwicker  # noqa: E402

C.low_priority()
t_start = time.time()
K = 0.01
out = {}

# =============================================================== Part A: existing histories
rows = []
CLS = json.loads((C.RESULTS / "exp11" / "classification.json").read_text(encoding="utf-8"))
for f in sorted(glob.glob(str(C.RESULTS / "exp11" / "runs" / "B_e*.json"))):
    d = json.loads(open(f, encoding="utf-8").read())
    sp, h = d["spec"], d["hist"]
    eps, dnu = sp["eps"], sp["dnu"]
    num = sp["rates"]["nu_m"]
    t = np.array([x["t"] for x in h])
    n = np.array([x["n"] for x in h])
    rm = np.array([x.get("rmean", 0.0) for x in h])
    asp = np.array([x.get("aspect", np.nan) for x in h])
    dissolved = n[-1] == 0
    cls = CLS[sp["name"]]["class"]
    # round stationary droplets only: exp11 class 'tek, kararlı', or 'bölündü' runs that END as stationary round
    # droplets; 'solucan' runs end as folded worms with R_eq ~ 40-50 and aspect ~ 1 and are excluded.
    stat = (cls == "tek, kararlı" or cls == "bölündü") and (not dissolved) and np.ptp(rm[-3:]) < 0.05 \
        and (np.isnan(asp[-1]) or asp[-1] < 1.15) and rm[-1] < 10
    rows.append(dict(name=sp["name"], eps=eps, dnu=dnu, nu_m=num, n_end=int(n[-1]), dissolved=bool(dissolved),
                     exp11_class=cls,
                     stationary=bool(stat), R_star=float(rm[-1]) if stat else None,
                     R_series=rm.tolist(), t_series=t.tolist()))
st = [r for r in rows if r["stationary"]]
eps_a = np.array([r["eps"] for r in st])
num_a = np.array([r["nu_m"] for r in st])
Rs = np.array([r["R_star"] for r in st])
LIB = {"1": lambda e, v: np.ones_like(e), "eps": lambda e, v: e, "nu-": lambda e, v: v,
       "eps/nu-": lambda e, v: e / v, "1/nu-": lambda e, v: 1 / v, "eps^2/nu-": lambda e, v: e**2 / v,
       "sqrt(eps/nu-)": lambda e, v: np.sqrt(e / v), "eps/nu-^2": lambda e, v: e / v**2}
names = list(LIB)
Th = np.column_stack([LIB[k](eps_a, num_a) for k in names])
cv = {}
for tau in np.logspace(-4, -0.5, 15):
    pr = C.loo_cv(Th, Rs, tau)
    cv[float(tau)] = dict(loo_rel_rmse=C.rel_rmse(Rs, pr), loo_r2=C.r2(Rs, pr),
                          support=[names[i] for i in np.flatnonzero(C.stlsq_scaled(Th, Rs, tau)[0])])
best = min(v["loo_rel_rmse"] for v in cv.values())
tau_a = min((t for t, v in cv.items() if v["loo_rel_rmse"] <= 1.1 * best), key=lambda t: (len(cv[t]["support"]), -t))
coef_a = C.stlsq_scaled(Th, Rs, tau_a)[0]
bs_a = C.bootstrap(Th, Rs, tau_a, n_boot=500)
# one-parameter flux-balance law R* = a eps / nu_- (forced) and its LOO error
a1 = float(np.sum(Rs * eps_a / num_a) / np.sum((eps_a / num_a) ** 2))
loo1 = np.array([np.sum(np.delete(Rs * eps_a / num_a, i)) / np.sum(np.delete((eps_a / num_a) ** 2, i)) * eps_a[i]
                 / num_a[i] for i in range(len(Rs))])
two_j = Rs * num_a  # = 2 j(eps) under flux balance


def pred_sparse(e, v):
    return np.column_stack([LIB[k](e, v) for k in names]) @ coef_a


diss = [r for r in rows if r["dissolved"]]
R_pred_diss = [float(pred_sparse(np.array([r["eps"]]), np.array([r["nu_m"]]))[0]) for r in diss]
worm = [r for r in rows if not r["stationary"] and not r["dissolved"]]
R_pred_worm = [float(pred_sparse(np.array([r["eps"]]), np.array([r["nu_m"]]))[0]) for r in worm]
# relaxation rate from the existing transients: lambda_i = -ln((R_{i+1}-R*)/(R_i-R*))/dt  for |R-R*| > 0.15
lam = []
for r in st:
    Rser, tser = np.array(r["R_series"]), np.array(r["t_series"])
    dR = Rser - r["R_star"]
    for i in range(len(dR) - 1):
        if abs(dR[i]) > 0.15 and abs(dR[i + 1]) > 0.05 and np.sign(dR[i]) == np.sign(dR[i + 1]):
            lam.append(dict(name=r["name"], eps=r["eps"], nu_m=r["nu_m"],
                            lam=float(-np.log(dR[i + 1] / dR[i]) / (tser[i + 1] - tser[i])),
                            R_from=float(Rser[i])))
lam_ratio = [x["lam"] / (x["nu_m"] / 2) for x in lam]
out["A_existing"] = dict(
    n_stationary=len(st), n_dissolved=len(diss), n_other=len(worm),
    stationary=[dict(name=r["name"], eps=r["eps"], dnu=r["dnu"], nu_m=r["nu_m"], R_star=r["R_star"]) for r in st],
    library=names, cv_path=cv, tau=tau_a, support=[names[i] for i in np.flatnonzero(coef_a)],
    coef={names[i]: float(coef_a[i]) for i in np.flatnonzero(coef_a)},
    bootstrap_p_incl={names[i]: float(bs_a["p_incl"][0, i]) for i in range(len(names))},
    bootstrap_std={names[i]: float(bs_a["std"][0, i]) for i in range(len(names))},
    forced_flux_balance=dict(a=a1, loo_rel_rmse=C.rel_rmse(Rs, loo1), loo_r2=C.r2(Rs, loo1),
                             two_j_over_eps=(two_j / eps_a).tolist()),
    dissolved=[dict(name=r["name"], R_pred=p) for r, p in zip(diss, R_pred_diss)],
    min_stationary_R=float(Rs.min()),
    not_round_or_not_stationary=[dict(name=r["name"], R_pred=p) for r, p in zip(worm, R_pred_worm)],
    relaxation=dict(points=lam, ratio_lambda_over_half_nu_m=dict(
        median=float(np.median(lam_ratio)) if lam else None, n=len(lam),
        p10=float(np.percentile(lam_ratio, 10)) if lam else None, p90=float(np.percentile(lam_ratio, 90)) if lam else None)))
print("A:", {k: out["A_existing"][k] for k in ("n_stationary", "support", "coef", "tau")},
      "forced", out["A_existing"]["forced_flux_balance"]["a"], out["A_existing"]["forced_flux_balance"]["loo_rel_rmse"],
      "relax", out["A_existing"]["relaxation"]["ratio_lambda_over_half_nu_m"], flush=True)

# =============================================================== Part B: new dense runs
CONDS = [(0.14, -0.002), (0.14, 0.0), (0.14, 0.002), (0.14, 0.004), (0.11, 0.0), (0.17, 0.0)]
R0S = [3.0, 8.0]
DT, T, REC = 0.2, 1500.0, 5.0
CACHE = C.OUT / "s2a_dense_runs.json"
runs = json.loads(CACHE.read_text(encoding="utf-8")) if CACHE.exists() else []
for eps, dnu in (CONDS if not runs else []):
    for R0 in R0S:
        t0 = time.time()
        rt = rates_zwicker(eps, dnu, k=K)
        m = ActiveCH(n=128, L=64.0, S=1.0, dtype=np.float32, **rt)
        c = m.disk(R0)
        c = c + 1e-3 * np.random.default_rng(int(R0 * 10 + eps * 1000)).standard_normal(c.shape)
        ts, Rsm, Rth = [], [], []

        def cb(t, cc, m=m, eps=eps, ts=ts, Rsm=Rsm, Rth=Rth):
            ex = np.clip((cc - eps) / (1 - eps), 0, None).sum() * m.dx**2   # smooth excess-area radius
            ts.append(t)
            Rsm.append(float(np.sqrt(ex / np.pi)))
            dd = m.droplets(cc)
            Rth.append(max(dd["radii"]) if dd["radii"] else 0.0)
            if Rth[-1] == 0.0 and len(Rth) > 3 and Rth[-2] == 0.0:
                return False

        m.run(c, T, DT, every=int(round(REC / DT)), cb=cb)
        runs.append(dict(eps=eps, dnu=dnu, nu_m=rt["nu_m"], R0=R0, t=ts, R=Rsm, R_thr=Rth, wall=time.time() - t0))
        print(f"run eps={eps} dnu={dnu} R0={R0}: R_end={Rsm[-1]:.3f} thr={Rth[-1]:.3f} ({time.time() - t0:.0f}s)",
              flush=True)

if not CACHE.exists():
    C.save_json("s2a_dense_runs.json", runs)
sfd = ps.SmoothedFiniteDifference(smoother_kws={"window_length": 9, "polyorder": 3})
RB, EB, VB_, DB, GB = [], [], [], [], []
for gi, r in enumerate(runs):
    t, R = np.array(r["t"]), np.array(r["R"])
    alive = np.array(r["R_thr"]) > 1.0
    if alive.sum() < 12:
        continue
    t, R = t[alive], R[alive]
    dR = sfd._differentiate(R[:, None], t)[:, 0]
    keep = slice(4, -4)
    RB.append(R[keep])
    DB.append(dR[keep])
    EB.append(np.full(len(R[keep]), r["eps"]))
    VB_.append(np.full(len(R[keep]), r["nu_m"]))
    GB.append(np.full(len(R[keep]), CONDS.index((r["eps"], r["dnu"]))))
RB, DB, EB, VB_, GB = map(np.concatenate, (RB, DB, EB, VB_, GB))
base = {"1": lambda R: np.ones_like(R), "R": lambda R: R, "1/R": lambda R: 1 / R, "R^2": lambda R: R**2,
        "1/R^2": lambda R: 1 / R**2}
par = {"": lambda e, v: np.ones_like(e), "eps·": lambda e, v: e, "nu-·": lambda e, v: v}
NB = [p + b for p in par for b in base]
ThB = np.column_stack([par[p](EB, VB_) * base[b](RB) for p in par for b in base])
# leave-one-condition-out validation along a threshold path
taus = np.logspace(-4, -0.3, 18)
lcoo = []
for tau in taus:
    pred = np.zeros_like(DB)
    for g in np.unique(GB):
        m_ = GB != g
        Cg = C.stlsq_scaled(ThB[m_], DB[m_], tau)[0]
        pred[~m_] = ThB[~m_] @ Cg
    Cf = C.stlsq_scaled(ThB, DB, tau)[0]
    lcoo.append(dict(tau=float(tau), r2_lcoo=C.r2(DB, pred), rel_rmse=C.rel_rmse(DB, pred),
                     support=[NB[i] for i in np.flatnonzero(Cf)], coef=Cf))
bestB = min(x["rel_rmse"] for x in lcoo)
selB = min((x for x in lcoo if x["rel_rmse"] <= 1.1 * bestB), key=lambda x: (len(x["support"]), -x["tau"]))
bsB = C.bootstrap(ThB, DB, selB["tau"], n_boot=200)
# the sharp test: flux-balance loss term -nu_- R / 2 forced into a minimal model dR/dt = a eps/R + b/R + c nu_- R + d
mini = ["eps·1/R", "1/R", "nu-·R", "1", "R"]
idx = [NB.index(x) for x in mini]
cm, *_ = np.linalg.lstsq(ThB[:, idx], DB, rcond=None)
out["B_dense_new_runs"] = dict(
    conditions=CONDS, R0=R0S, dt=DT, T=T, rec=REC, L=64.0, n=128,
    runs=[dict(eps=r["eps"], dnu=r["dnu"], R0=r["R0"], R_end=r["R"][-1], R_thr_end=r["R_thr"][-1], wall=r["wall"])
          for r in runs],
    n_points=len(DB), library=NB,
    path=[dict(tau=x["tau"], r2_lcoo=x["r2_lcoo"], n_terms=len(x["support"]), support=x["support"]) for x in lcoo],
    selected=dict(tau=selB["tau"], support=selB["support"], r2_lcoo=selB["r2_lcoo"],
                  coef={NB[i]: float(selB["coef"][i]) for i in np.flatnonzero(selB["coef"])},
                  bootstrap_p_incl={NB[i]: float(bsB["p_incl"][0, i]) for i in range(len(NB))},
                  bootstrap_std={NB[i]: float(bsB["std"][0, i]) for i in range(len(NB))}),
    minimal_ols=dict(zip(mini, map(float, cm))))
print("B selected:", out["B_dense_new_runs"]["selected"], flush=True)
print("B minimal OLS:", out["B_dense_new_runs"]["minimal_ols"], flush=True)
C.save_json("s2a_droplets.json", dict(out, wall_s=time.time() - t_start))

fig, axs = plt.subplots(1, 3, figsize=(15, 4.3))
xx = eps_a / num_a
axs[0].plot(xx, Rs, "ko", label="exp11 durağan damlalar")
for r, p in zip(diss, R_pred_diss):
    axs[0].plot(r["eps"] / r["nu_m"], 0.3, "rx")
xs = np.linspace(0, xx.max() * 1.1, 50)
axs[0].plot(xs, a1 * xs, "b-", lw=1, label=f"R* = {a1:.3f} ε/ν₋")
axs[0].set_xlabel("ε/ν₋")
axs[0].set_ylabel("R* (w)")
axs[0].set_title("durağan yarıçap (x: çözünen koşular)", fontsize=9)
axs[0].legend(fontsize=8)
for r in runs:
    axs[1].plot(r["t"], r["R"], lw=1, label=f"ε={r['eps']},Δν={r['dnu']},R0={r['R0']:g}")
axs[1].set_xlabel("t (t₀)")
axs[1].set_ylabel("R (pürüzsüz fazlalık-alan yarıçapı)")
axs[1].legend(fontsize=5, ncol=2)
axs[1].set_title("yeni yoğun örneklenmiş koşular (Part B)", fontsize=9)
coefsel = selB["coef"]
for g in np.unique(GB):
    m_ = GB == g
    o = np.argsort(RB[m_])
    axs[2].plot(RB[m_][o], DB[m_][o], ".", ms=2)
    axs[2].plot(RB[m_][o], (ThB[m_] @ coefsel)[o], "k-", lw=0.7)
axs[2].axhline(0, color="gray", lw=0.5)
axs[2].set_xlabel("R")
axs[2].set_ylabel("dR/dt")
axs[2].set_title("dR/dt: noktalar veri, çizgi seçilen seyrek model", fontsize=9)
fig.tight_layout()
fig.savefig(C.OUT / "fig_s2a_droplets.png", dpi=110)
print(f"done {time.time() - t_start:.0f}s")
