"""Exp16 step 1a-2: noise robustness of the Gray-Scott rediscovery with a REDUCED library (no polynomial x derivative
products: 1 + 9 polynomials of (u, v) up to degree 3 + 10 spatial derivatives = 20 terms), strong vs weak form.
Same simulated data as s1a (seeds 1/2, warm-up 1500/2500, 81 snapshots, every step). Noise: 1 % and 5 % of the field
std. Selection on the held-out run with rule B (sparsest with val R^2 >= max - 1e-3).
Output: results/exp16_equations/s1a2_weak_reduced.json"""
import common as C  # noqa: I001

import time

import numpy as np
import pysindy as ps

from qcell.grayscott import GrayScott  # noqa: E402

C.low_priority()
F, K, DU, DV, N, H = 0.026, 0.059, 2e-5, 1e-5, 128, 2.5 / 256
TRUE = [{"1": F, "u": -F, "u v^2": -1.0, "u_11": DU, "u_22": DU},
        {"v": -(F + K), "u v^2": 1.0, "v_11": DV, "v_22": DV}]
TAUS = np.logspace(-4, -0.2, 16)
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
    return X


X_tr, X_va = simulate(1, 1500, 81), simulate(2, 2500, 81)
xg = np.arange(N) * H
grid = np.stack(np.meshgrid(xg, xg, indexing="ij"), axis=-1)
xyt = np.zeros((N, N, 81, 3))
xyt[..., 0], xyt[..., 1] = grid[:, :, None, 0], grid[:, :, None, 1]
xyt[..., 2] = np.arange(81)[None, None, :]


def weak(X, seed, K_dom=2000, hw=12):
    np.random.seed(seed)
    wl = ps.WeakPDELibrary(function_library=ps.PolynomialLibrary(degree=3, include_bias=False), derivative_order=2,
                           spatiotemporal_grid=xyt, include_bias=True, include_interaction=False, K=K_dom,
                           H_xt=[hw * H, hw * H, hw], p=4)
    th = wl.fit_transform(X)
    th = np.asarray(th[0] if isinstance(th, list) else th)
    return th, np.asarray(wl.convert_u_dot_integral(X)), wl.get_feature_names(["u", "v"])


def strong(X, seed, rows=1000):
    lib = ps.PDELibrary(function_library=ps.PolynomialLibrary(degree=3, include_bias=False), derivative_order=2,
                        spatial_grid=grid, include_bias=True, include_interaction=False,
                        differentiation_method=ps.FiniteDifference, diff_kwargs={"periodic": True, "is_uniform": True})
    lib.fit(X[:, :, :2])
    th = np.asarray(lib.transform(X))[:, :, 1:-1]
    yd = ps.FiniteDifference(order=2, d=1, axis=2, is_uniform=True)._differentiate(X, 1.0)[:, :, 1:-1]
    rng = np.random.default_rng(seed)
    T = th.shape[2]
    th, yd = th.reshape(N * N * T, -1), yd.reshape(N * N * T, 2)
    r = rng.choice(len(th), rows * T, replace=False)
    return th[r], yd[r], lib.get_feature_names(["u", "v"])


out = dict(truth=TRUE, variants=[])
for form in ("weak", "strong"):
    for noise in (0.0, 0.01, 0.05):
        t0 = time.time()
        rng = np.random.default_rng(int(noise * 1000) + 7)
        Xtr = X_tr + noise * X_tr.reshape(-1, 2).std(axis=0) * rng.standard_normal(X_tr.shape)
        Xva = X_va + noise * X_va.reshape(-1, 2).std(axis=0) * rng.standard_normal(X_va.shape)
        f = weak if form == "weak" else strong
        Ttr, Ytr, names = f(Xtr, 11)
        Tva, Yva, _ = f(Xva, 12)
        path = C.threshold_path(Ttr, Ytr, Tva, Yva, TAUS)
        res = dict(form=form, noise=noise, n_terms_library=len(names), eq=[])
        for j in range(2):
            sel = C.select_model_r2(path, j, r2_tol=1e-3)
            c = sel["coef"][j]
            found = {names[i]: float(c[i]) for i in np.flatnonzero(c)}
            rel = {n: (found[n] - v) / abs(v) for n, v in TRUE[j].items() if n in found}
            exact_taus = [r["tau"] for r in path if set(names[i] for i in np.flatnonzero(r["coef"][j])) == set(TRUE[j])]
            res["eq"].append(dict(tau=sel["tau"], r2_val=sel["r2_val"][j], n_terms=sel["n_terms"][j], found=found,
                                  exact_support=set(found) == set(TRUE[j]),
                                  false_terms=[n for n in found if n not in TRUE[j]],
                                  missing=[n for n in TRUE[j] if n not in found], rel_err=rel,
                                  max_abs_rel_err=max([abs(v) for v in rel.values()], default=None),
                                  tau_range_exact=[min(exact_taus), max(exact_taus)] if exact_taus else None,
                                  equation=C.eq_string(names, c, 5)))
        res["wall_s"] = time.time() - t0
        out["variants"].append(res)
        print(f"{form} noise={noise}: " + " | ".join(
            f"{['u', 'v'][j]}: n={e['n_terms']} exact={e['exact_support']} maxrel={e['max_abs_rel_err']} "
            f"R2={e['r2_val']:.5f} tau_ok={e['tau_range_exact']}" for j, e in enumerate(res["eq"])) +
            f" ({res['wall_s']:.0f}s)", flush=True)
        for e in res["eq"]:
            if e["n_terms"] <= 10:
                print("    ", e["equation"], flush=True)
C.save_json("s1a2_weak_reduced.json", out)
print(f"done {time.time() - t_start:.0f}s")
