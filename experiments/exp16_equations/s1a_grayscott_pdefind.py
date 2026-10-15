"""Exp16 step 1a (ground truth): can PySINDy PDE-FIND rediscover the Gray-Scott equations from simulated fields?

Truth (qcell.grayscott, Pearson 1993):  u_t = Du (u_11 + u_22) - u v^2 + F (1 - u),  v_t = Dv (v_11 + v_22) + u v^2 - (F + k) v
F = 0.026, k = 0.059, Du = 2e-5, Dv = 1e-5, h = 2.5/256, dt = 1 (forward Euler), n = 128, periodic.
Library: PDELibrary, polynomials of (u, v) up to degree 3, all 1st/2nd spatial derivatives (incl. u_12) and all products
(110 candidate terms per equation). Training run and a held-out validation run (different noise seed, different time
window). The sparsity threshold is chosen on the validation run only (sparsest model within 5 % of the best validation
error); the truth is used only for scoring afterwards. Variants: time derivative from forward differences (consistent with
the integrator), from 2nd-order central differences at snapshot spacing s = 1, 2, 5, noisy data (strong form) and the weak
form (WeakPDELibrary) with noise. Forecast test: the identified PDE is integrated from the validation state for 200 steps.
Output: results/exp16_equations/s1a_grayscott.json, fig_s1a_grayscott.png."""
import common as C  # noqa: I001  (must be first: thread limits)

import time

import matplotlib
import numpy as np
import pysindy as ps

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from qcell.grayscott import GrayScott  # noqa: E402

C.low_priority()
F, K, DU, DV, N, H = 0.026, 0.059, 2e-5, 1e-5, 128, 2.5 / 256
TRUE = [{"1": F, "u": -F, "u v^2": -1.0, "u_11": DU, "u_22": DU},
        {"v": -(F + K), "u v^2": 1.0, "v_11": DV, "v_22": DV}]
TAUS = np.logspace(-4, -0.2, 16)
ROWS_PER_SLICE = 1000
t_start = time.time()


def simulate(seed, t_warm, n_rec):
    gs = GrayScott(n=N, F=F, k=K, Du=DU, Dv=DV, h=H)
    U, V = gs.seed(np.random.default_rng(seed), half=10, noise=0.01)
    U, V = gs.run(U, V, t_warm)
    X = np.empty((N, N, n_rec, 2))
    X[:, :, 0, 0], X[:, :, 0, 1] = U, V
    for i in range(1, n_rec):
        U, V = gs.run(U, V, 1)
        X[:, :, i, 0], X[:, :, i, 1] = U, V
    return gs, X


gs, X_tr = simulate(1, 1500, 81)
_, X_va = simulate(2, 2500, 81)
cover = {"train_area_frac_v>0.1": float((X_tr[:, :, 0, 1] > 0.1).mean()),
         "val_area_frac_v>0.1": float((X_va[:, :, 0, 1] > 0.1).mean()),
         "train_spots": gs.spots(X_tr[:, :, 0, 1])["n"], "val_spots": gs.spots(X_va[:, :, 0, 1])["n"]}
print("data", cover, f"{time.time() - t_start:.1f}s", flush=True)

xg = np.arange(N) * H
grid = np.stack(np.meshgrid(xg, xg, indexing="ij"), axis=-1)


def make_lib():
    return ps.PDELibrary(function_library=ps.PolynomialLibrary(degree=3, include_bias=False), derivative_order=2,
                         spatial_grid=grid, include_bias=True, include_interaction=True,
                         differentiation_method=ps.FiniteDifference,
                         diff_kwargs={"periodic": True, "is_uniform": True})


LIB = make_lib().fit(X_tr[:, :, :2])
NAMES = LIB.get_feature_names(["u", "v"])
for j in range(2):
    assert all(n in NAMES for n in TRUE[j]), (TRUE[j], NAMES[:25])
print("library size", len(NAMES), flush=True)


def theta_full(X):
    out = LIB.transform(X)
    return np.asarray(out[0] if isinstance(out, list) else out)


def strong_rows(X, s, deriv, seed):
    """Subsampled (Theta, u_t) rows. X already noisy if needed; s = snapshot spacing in steps (dt = s)."""
    Xs = X[:, :, ::s]
    T = Xs.shape[2]
    if deriv == "fwd":
        Yd = (Xs[:, :, 1:] - Xs[:, :, :-1]) / s
        tid = np.arange(T - 1)
    else:  # PySINDy 2nd-order central finite differences along time, interior points only
        Yd = ps.FiniteDifference(order=2, d=1, axis=2, is_uniform=True)._differentiate(Xs, float(s))[:, :, 1:-1]
        tid = np.arange(1, T - 1)
    rng = np.random.default_rng(seed)
    TH, YY = [], []
    for c0 in range(0, len(tid), 8):
        cid = tid[c0:c0 + 8]
        th = theta_full(Xs[:, :, cid]).reshape(N * N, len(cid), -1)
        yd = Yd[:, :, c0:c0 + 8].reshape(N * N, len(cid), 2)
        for i in range(len(cid)):
            r = rng.choice(N * N, ROWS_PER_SLICE, replace=False)
            TH.append(th[r, i])
            YY.append(yd[r, i])
    return np.vstack(TH), np.vstack(YY)


def weak_rows(X, s, seed, K_dom=1000, hw=8):
    Xs = X[:, :, ::s]
    T = Xs.shape[2]
    xyt = np.zeros((N, N, T, 3))
    xyt[..., 0], xyt[..., 1] = grid[:, :, None, 0], grid[:, :, None, 1]
    xyt[..., 2] = (np.arange(T) * s)[None, None, :]
    np.random.seed(seed)  # WeakPDELibrary draws its subdomain centres with np.random
    wl = ps.WeakPDELibrary(function_library=ps.PolynomialLibrary(degree=3, include_bias=False), derivative_order=2,
                           spatiotemporal_grid=xyt, include_bias=True, include_interaction=True, K=K_dom,
                           H_xt=[hw * H, hw * H, hw * s], p=4)
    th = wl.fit_transform(Xs)
    th = np.asarray(th[0] if isinstance(th, list) else th)
    yd = np.asarray(wl.convert_u_dot_integral(Xs))
    assert wl.get_feature_names(["u", "v"]) == NAMES
    return th, yd


def add_noise(X, frac, seed):
    if frac == 0:
        return X
    rng = np.random.default_rng(seed)
    sd = X.reshape(-1, 2).std(axis=0)
    return X + frac * sd * rng.standard_normal(X.shape)


def score(Csel):
    out = []
    for j in range(2):
        found = {NAMES[i]: float(Csel[j, i]) for i in np.flatnonzero(Csel[j])}
        tr = TRUE[j]
        rel = {n: (found[n] - v) / abs(v) for n, v in tr.items() if n in found}
        out.append(dict(found=found, missing=[n for n in tr if n not in found],
                        false_terms=[n for n in found if n not in tr], rel_err=rel,
                        max_abs_rel_err=max([abs(v) for v in rel.values()], default=None),
                        exact_support=(set(found) == set(tr))))
    return out


def forecast(Csel, steps=200):
    """Integrate the identified PDE (forward Euler, dt = 1, same discrete library) from the validation state and
    compare with the true Gray-Scott trajectory (relative L2 error of v)."""
    x = X_va[:, :, :1].copy()
    U, V = X_va[:, :, 0, 0].copy(), X_va[:, :, 0, 1].copy()
    errs = []
    for st in range(1, steps + 1):
        th = theta_full(x).reshape(N * N, -1)
        x = x + (th @ Csel.T).reshape(N, N, 1, 2)
        U, V = gs.run(U, V, 1)
        if not np.all(np.isfinite(x)) or np.abs(x).max() > 1e3:
            errs.append(np.inf)
            break
        if st % 50 == 0:
            errs.append(float(np.linalg.norm(x[:, :, 0, 1] - V) / np.linalg.norm(V)))
    return errs, x[:, :, 0, 1], V


VARIANTS = [
    dict(name="ileri fark, s=1 (integratörle tutarlı)", form="strong", s=1, deriv="fwd", noise=0.0),
    dict(name="merkezi FD, s=1", form="strong", s=1, deriv="central", noise=0.0),
    dict(name="merkezi FD, s=2", form="strong", s=2, deriv="central", noise=0.0),
    dict(name="merkezi FD, s=5", form="strong", s=5, deriv="central", noise=0.0),
    dict(name="merkezi FD, s=1, %1 gürültü", form="strong", s=1, deriv="central", noise=0.01),
    dict(name="zayıf form, s=1, %1 gürültü", form="weak", s=1, noise=0.01),
    dict(name="zayıf form, s=2, %5 gürültü", form="weak", s=2, noise=0.05),
]
results = []
paths_for_fig = {}
for vi, v in enumerate(VARIANTS):
    t0 = time.time()
    Xtr, Xva = add_noise(X_tr, v["noise"], 10 + vi), add_noise(X_va, v["noise"], 20 + vi)
    if v["form"] == "strong":
        Ttr, Ytr = strong_rows(Xtr, v["s"], v["deriv"], 100 + vi)
        Tva, Yva = strong_rows(Xva, v["s"], v["deriv"], 200 + vi)
    else:
        Ttr, Ytr = weak_rows(Xtr, v["s"], 100 + vi)
        Tva, Yva = weak_rows(Xva, v["s"], 200 + vi)
    path = C.threshold_path(Ttr, Ytr, Tva, Yva, TAUS)
    res = dict(v, n_rows_train=len(Ttr), n_rows_val=len(Tva), rules={})
    for rule, fsel in (("A_rmse5pct", lambda j: C.select_model(path, j)),
                       ("B_dR2_1e-3", lambda j: C.select_model_r2(path, j, r2_tol=1e-3))):
        sel = [fsel(j) for j in range(2)]
        Csel = np.vstack([sel[0]["coef"][0], sel[1]["coef"][1]])
        sc = score(Csel)
        errs, vf, vt = forecast(Csel)
        res["rules"][rule] = dict(tau=[sel[0]["tau"], sel[1]["tau"]], r2_val=[sel[0]["r2_val"][0], sel[1]["r2_val"][1]],
                                  n_terms=[sel[0]["n_terms"][0], sel[1]["n_terms"][1]],
                                  equations=[C.eq_string(NAMES, Csel[j], 5) for j in range(2)], score=sc,
                                  forecast_relL2_v_every50=errs)
        if rule.startswith("B"):
            paths_for_fig[v["name"]] = (path, vf, vt)
        print(f"[{vi}] {v['name']} | rule {rule}: terms {res['rules'][rule]['n_terms']} "
              f"exact {[s['exact_support'] for s in sc]} maxrel {[s['max_abs_rel_err'] for s in sc]} "
              f"R2val {res['rules'][rule]['r2_val']} forecast {errs}", flush=True)
        if max(res["rules"][rule]["n_terms"]) <= 12:
            for j in range(2):
                print("   ", ["u_t", "v_t"][j], "=", res["rules"][rule]["equations"][j], flush=True)
    # threshold sensitivity: range of tau that returns exactly the true support (for each equation)
    res["path"] = []
    for r in path:
        ex = [set(NAMES[i] for i in np.flatnonzero(r["coef"][j])) == set(TRUE[j]) for j in range(2)]
        res["path"].append(dict(tau=r["tau"], n_terms=r["n_terms"], r2_val=r["r2_val"], exact=ex))
    for j in range(2):
        ok = [p["tau"] for p in res["path"] if p["exact"][j]]
        res.setdefault("tau_range_exact_support", []).append([min(ok), max(ok)] if ok else None)
    res["wall_s"] = time.time() - t0
    results.append(res)
    print(f"    tau range with exact support: {res['tau_range_exact_support']} ({res['wall_s']:.0f}s)", flush=True)

C.save_json("s1a_grayscott.json", dict(truth=TRUE, params=dict(F=F, k=K, Du=DU, Dv=DV, n=N, h=H, dt=1.0),
                                       data=cover, library_size=len(NAMES), library_names=NAMES,
                                       selection_rules={"A_rmse5pct": "sparsest with val rel-RMSE <= 1.05 x min",
                                                        "B_dR2_1e-3": "sparsest with val R2 >= max R2 - 1e-3"},
                                       variants=results, wall_total_s=time.time() - t_start))

fig, axs = plt.subplots(1, 3, figsize=(15, 4.4))
for name, (path, _, _) in paths_for_fig.items():
    axs[0].plot([r["n_terms"][0] for r in path], [r["r2_val"][0] for r in path], "o-", ms=3, lw=1, label=name)
axs[0].axvline(5, color="k", ls=":", lw=0.8)
axs[0].set_xscale("log")
axs[0].set_ylim(0.5, 1.005)
axs[0].set_xlabel("u_t denklemindeki terim sayısı")
axs[0].set_ylabel("doğrulama R² (ayrı koşu)")
axs[0].set_title("Pareto: eşik taraması (doğru model 5 terim)", fontsize=9)
axs[0].legend(fontsize=6)
names_plot = list(TRUE[0]) + [n for n in TRUE[1] if n not in TRUE[0]]
labels = [r["name"] for r in results]
for vi, r in enumerate(results):
    scB = r["rules"]["B_dR2_1e-3"]["score"]
    errs_u = [abs(scB[0]["rel_err"].get(n, np.nan)) if n in TRUE[0] else np.nan for n in names_plot]
    errs_v = [abs(scB[1]["rel_err"].get(n, np.nan)) if n in TRUE[1] else np.nan for n in names_plot]
    e = np.where(np.isnan(errs_u), errs_v, errs_u)
    axs[1].semilogy(np.arange(len(names_plot)) + 0.1 * vi, np.maximum(e, 1e-16), "o", ms=4, label=labels[vi])
axs[1].set_xticks(range(len(names_plot)))
axs[1].set_xticklabels(names_plot, rotation=45, fontsize=7)
axs[1].set_ylabel("|bağıl katsayı hatası| (gerçek terimler)")
axs[1].set_title("katsayı hataları, kural B (u ve v denklemleri birlikte)", fontsize=9)
_, vf, vt = paths_for_fig[VARIANTS[1]["name"]]
axs[2].imshow(np.hstack([vt, vf]), cmap="viridis")
axs[2].set_title("200 adım tahmin: sol gerçek v, sağ bulunan PDE (merkezi FD s=1, kural B)", fontsize=9)
axs[2].axis("off")
fig.tight_layout()
fig.savefig(C.OUT / "fig_s1a_grayscott.png", dpi=110)
print(f"done {time.time() - t_start:.0f}s")
