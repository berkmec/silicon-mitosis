"""Exp16 step 2d: scaling laws from small datasets (honest: few points, wide uncertainty).

(a) exp11 3D active droplets (b3d_dnu*.json, eps = 0.2, k = 0.01): parent radius just before division and daughter
    radius at birth vs nu_- = dnu + k (1 - eps). Flux-balance limits: R ~ nu_-^(-1/2) (R << l) ... nu_-^(-1) (R >> l),
    l = sqrt(D/k) = 10 w. Power-law fit on 4 points, LOO, bootstrap over points.
(b) exp8 cluster fission threshold g_c(N), d = inf / 4 / 2: g_c taken from the zero crossing of the energy gain
    E_single - E_split(g) (linear interpolation in log g; removes the x1.45 grid quantisation), power law fits.
    Liquid-drop (Rayleigh) expectation for unscreened 2D clusters: g_c ~ 1/N.
(c) exp11 gs_ramp: Gray-Scott plateau spot number and mean spot area vs diffusion scale D (dimensional analysis:
    spot area ~ D, N ~ 1/D at fixed box).
Output: results/exp16_equations/s2d_scaling.json, fig_s2d_scaling.png."""
import common as C  # noqa: I001

import json

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

C.low_priority()
E11 = C.RESULTS / "exp11"
out = {}
rng = np.random.default_rng(0)


def powerfit(x, y):
    """log-log LS fit y = A x^p; LOO relative errors and a bootstrap over points for p."""
    lx, ly = np.log(x), np.log(y)
    p, lA = np.polyfit(lx, ly, 1)
    loo = []
    for i in range(len(x)):
        m = np.arange(len(x)) != i
        pi, lAi = np.polyfit(lx[m], ly[m], 1)
        loo.append(float(np.exp(lAi + pi * lx[i]) / y[i] - 1))
    bs = []
    for _ in range(2000):
        idx = rng.integers(0, len(x), len(x))
        if len(np.unique(x[idx])) >= 2:
            bs.append(np.polyfit(lx[idx], ly[idx], 1)[0])
    resid = ly - (lA + p * lx)
    se = np.sqrt(np.sum(resid**2) / max(len(x) - 2, 1) / np.sum((lx - lx.mean()) ** 2)) if len(x) > 2 else np.nan
    return dict(p=float(p), A=float(np.exp(lA)), p_se=float(se), p_boot_p5_p95=[float(np.percentile(bs, 5)),
                                                                                 float(np.percentile(bs, 95))],
                loo_rel_err=loo, max_abs_loo=float(np.max(np.abs(loo))), n=len(x))


# (a) 3D division
rows = []
for fn, dnu in (("b3d_dnum0p002", -0.002), ("b3d_dnu0", 0.0), ("b3d_dnu0p002", 0.002), ("b3d_dnu0p004", 0.004),
                ("b3d_dnum0p004", -0.004)):
    d = json.loads((E11 / f"{fn}.json").read_text(encoding="utf-8"))
    h = d["hist"]
    i2 = next((i for i, x in enumerate(h) if x["n"] == 2), None)
    num = dnu + 0.01 * (1 - 0.2)
    if i2 is None:
        rows.append(dict(dnu=dnu, nu_m=num, divided=False, R_final=h[-1]["radii"][0]))
        continue
    rows.append(dict(dnu=dnu, nu_m=num, divided=True, t_div=h[i2]["t"], R_parent=h[i2 - 1]["radii"][0],
                     R_daughter=float(np.mean(h[i2]["radii"])), vol_ratio=float(np.mean(h[i2]["radii"]) ** 3 * 2
                                                                                / h[i2 - 1]["radii"][0] ** 3)))
dv = [r for r in rows if r["divided"]]
num = np.array([r["nu_m"] for r in dv])
out["a_3d_division"] = dict(rows=rows,
                            daughter=powerfit(num, np.array([r["R_daughter"] for r in dv])),
                            parent=powerfit(num, np.array([r["R_parent"] for r in dv])),
                            expectation="flux balance: exponent -1/2 for R << l = 10 w, -1 for R >> l")
print("(a)", {k: (v["p"], v["p_boot_p5_p95"], v["max_abs_loo"]) for k, v in out["a_3d_division"].items()
              if isinstance(v, dict) and "p" in v}, flush=True)

# (b) exp8 fission thresholds
E8 = C.RESULTS / "exp8"
gc = {}
for fn, dlab in (("fission_scan.json", "inf"), ("fission_scan_d4.0.json", "4.0"), ("fission_scan_d2.0.json", "2.0")):
    sc = json.loads((E8 / fn).read_text(encoding="utf-8"))["scan"]
    for N in sorted(set(x["N"] for x in sc)):
        pts = sorted([(x["g"], x["gain"]) for x in sc if x["N"] == N], key=lambda p: p[0])
        g, gain = np.array(pts).T
        k = next((i for i in range(len(g)) if gain[i] > 0), None)
        if k is None or k == 0:
            continue
        lg0, lg1 = np.log(g[k - 1]), np.log(g[k])
        lgc = lg0 + (0 - gain[k - 1]) * (lg1 - lg0) / (gain[k] - gain[k - 1])
        gc.setdefault(dlab, []).append(dict(N=int(N), g_c_interp=float(np.exp(lgc)), g_lo=float(g[k - 1]),
                                            g_hi=float(g[k])))
out["b_exp8_fission"] = {}
for dlab, lst in gc.items():
    N = np.array([x["N"] for x in lst], float)
    g = np.array([x["g_c_interp"] for x in lst])
    gh = np.array([x["g_hi"] for x in lst])
    out["b_exp8_fission"][dlab] = dict(points=lst, fit_interp=powerfit(N, g), fit_grid_first_fission=powerfit(N, gh))
    print("(b) d =", dlab, "p_interp", round(out["b_exp8_fission"][dlab]["fit_interp"]["p"], 3),
          out["b_exp8_fission"][dlab]["fit_interp"]["p_boot_p5_p95"], "p_grid",
          round(out["b_exp8_fission"][dlab]["fit_grid_first_fission"]["p"], 3), flush=True)

# (c) gs_ramp
r = json.loads((E11 / "gs_ramp.json").read_text(encoding="utf-8"))
st = {}
for h in r["hist"]:
    st.setdefault(h["D_scale"], []).append(h)
Ds, Np, Ap = [], [], []
for D, hs in sorted(st.items(), reverse=True):
    tail = hs[-6:]
    Ds.append(D)
    Np.append(float(np.mean([x["n"] for x in tail])))
    Ap.append(float(np.mean([x["mean_spot_area"] for x in tail])))
Ds, Np, Ap = map(np.array, (Ds, Np, Ap))
out["c_gs_ramp"] = dict(D_scale=Ds.tolist(), N_plateau=Np.tolist(), mean_spot_area=Ap.tolist(),
                        N_vs_D=powerfit(Ds, Np), N_vs_D_without_first_stage=powerfit(Ds[1:], Np[1:]),
                        area_vs_D=powerfit(Ds, Ap), expectation="area ~ D^1, N ~ D^-1 (diffusion length ~ sqrt(D))")
print("(c)", {k: (v["p"], v["p_boot_p5_p95"]) for k, v in out["c_gs_ramp"].items() if isinstance(v, dict)}, flush=True)
C.save_json("s2d_scaling.json", out)

fig, axs = plt.subplots(1, 3, figsize=(15, 4.2))
for key, mk in (("R_daughter", "o"), ("R_parent", "s")):
    axs[0].loglog(num, [r[key] for r in dv], mk, label=key)
xx = np.linspace(num.min() * 0.9, num.max() * 1.1, 20)
for key, f in (("daughter", out["a_3d_division"]["daughter"]), ("parent", out["a_3d_division"]["parent"])):
    axs[0].loglog(xx, f["A"] * xx ** f["p"], "-", lw=0.8, label=f"{key}: p={f['p']:.2f}")
axs[0].set_xlabel("ν₋ = Δν + k(1−ε)")
axs[0].set_ylabel("yarıçap (w)")
axs[0].set_title("3B bölünme (4 nokta)", fontsize=9)
axs[0].legend(fontsize=7)
for dlab, v in out["b_exp8_fission"].items():
    N = np.array([x["N"] for x in v["points"]])
    axs[1].loglog(N, [x["g_c_interp"] for x in v["points"]], "o-", ms=4, lw=0.8,
                  label=f"d={dlab}: p={v['fit_interp']['p']:.2f}")
axs[1].set_xlabel("N")
axs[1].set_ylabel("g_c (kazanç sıfır geçişi)")
axs[1].set_title("exp8 küme fisyon eşiği", fontsize=9)
axs[1].legend(fontsize=7)
axs[2].loglog(Ds, Np, "ko-", label=f"N: p={out['c_gs_ramp']['N_vs_D']['p']:.2f}")
axs[2].loglog(Ds, Ap * 1e4, "bs-", label=f"alan×1e4: p={out['c_gs_ramp']['area_vs_D']['p']:.2f}")
axs[2].set_xlabel("D ölçeği")
axs[2].set_title("Gray–Scott: leke sayısı/alanı – D", fontsize=9)
axs[2].legend(fontsize=7)
fig.tight_layout()
fig.savefig(C.OUT / "fig_s2d_scaling.png", dpi=110)
