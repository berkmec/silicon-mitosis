"""Exp16 step 2b: growth law of the Gray-Scott spot count N(t) (exp11: runs/GS_*.json, pilotGS.json, gs_track.json).

Candidates (each 3 parameters r, K, N0, integrated with solve_ivp and fitted on sqrt(N), a variance-stabilising scale
for counts):  logistic r N (1 - N/K);  Gompertz r N ln(K/N);  front-limited (2D rim growth) r sqrt(N) (1 - N/K);
saturating linear r (1 - N/K).  Compared by BIC and 5-fold CV (interleaved time points). Independent check with PySINDy
(SmoothedFiniteDifference + STLSQ, library {N, N^2, sqrt N, N^1.5} on N/N_max) and the doubling time ln2/r.
Output: results/exp16_equations/s2b_gs_count.json, fig_s2b_gs_count.png."""
import common as C  # noqa: I001

import json

import matplotlib
import numpy as np
import pysindy as ps
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

C.low_priority()
E11 = C.RESULTS / "exp11"
series = {}
for fk in ("F0.010_k0.050", "F0.018_k0.055", "F0.026_k0.059", "F0.034_k0.063", "F0.042_k0.063"):
    h = json.loads((E11 / "runs" / f"GS_{fk}.json").read_text(encoding="utf-8"))["hist"]
    series[f"GS_{fk} (128², L=1.25)"] = (np.array([x["t"] for x in h], float), np.array([x["n"] for x in h], float))
pg = np.array(json.loads((E11 / "pilotGS.json").read_text(encoding="utf-8")), float)
series["pilotGS F0.026_k0.059 (256², L=2.5)"] = (pg[:, 0], pg[:, 1])
tr = json.loads((E11 / "gs_track.json").read_text(encoding="utf-8"))
series["gs_track (ölçek değişimi sonrası, 256²)"] = (np.array([x["t"] for x in tr["hist"]], float),
                                                     np.array([x["n"] for x in tr["hist"]], float))
MODELS = {
    "lojistik": lambda N, r, K: r * N * (1 - N / K),
    "Gompertz": lambda N, r, K: r * N * np.log(K / np.maximum(N, 1e-9)),
    "cephe (sqrt N)": lambda N, r, K: r * np.sqrt(np.maximum(N, 0)) * (1 - N / K),
    "doğrusal doyum": lambda N, r, K: r * (1 - N / K),
}


def integrate(f, p, t):
    r, K, N0 = p
    s = solve_ivp(lambda tt, y: f(y, r, K), (t[0], t[-1]), [N0], t_eval=t, rtol=1e-8, atol=1e-10)
    if not s.success or s.y.shape[1] != len(t):
        return np.full_like(t, np.nan)
    return s.y[0]


def fit(f, t, N, mask):
    Kg, rg = N.max(), 3.0 / (t[-1] - t[0])
    best = None
    for r0 in (rg, 3 * rg, 0.3 * rg):
        try:
            res = least_squares(lambda p: np.nan_to_num(np.sqrt(np.maximum(integrate(f, p, t), 0))[mask]
                                                        - np.sqrt(N[mask]), nan=1e3),
                                x0=[r0, Kg, max(N[0], 1.0)], bounds=([1e-7, 1.0, 0.1], [1.0, 10 * Kg, 10 * Kg]))
            if best is None or res.cost < best.cost:
                best = res
        except ValueError:
            pass
    return best.x


out = {}
fig, axs = plt.subplots(1, 3, figsize=(15, 4.3))
for si, (name, (t, N)) in enumerate(series.items()):
    t = t - t[0]
    res = {}
    for mname, f in MODELS.items():
        p = fit(f, t, N, np.ones_like(t, bool))
        pred = integrate(f, p, t)
        rss = float(np.sum((np.sqrt(np.maximum(pred, 0)) - np.sqrt(N)) ** 2))
        bic = len(t) * np.log(rss / len(t)) + 3 * np.log(len(t))
        cvp = np.zeros_like(N)
        for k in range(5):
            m = np.arange(len(t)) % 5 != k
            m[0] = True
            pk = fit(f, t, N, m)
            cvp[~m] = integrate(f, pk, t)[~m]
        res[mname] = dict(r=float(p[0]), K=float(p[1]), N0=float(p[2]), bic=float(bic), r2=C.r2(N, pred),
                          cv_rmse_N=float(np.sqrt(np.mean((cvp[1:] - N[1:]) ** 2))),
                          doubling_time=float(np.log(2) / p[0]) if mname in ("lojistik", "Gompertz") else None)
    bmin = min(v["bic"] for v in res.values())
    for v in res.values():
        v["dBIC"] = v["bic"] - bmin
    # PySINDy cross-check on N / N_max
    Nn = N / N.max()
    lib = ps.CustomLibrary(library_functions=[lambda x: x, lambda x: x**2, lambda x: np.sqrt(np.abs(x)),
                                              lambda x: np.abs(x) ** 1.5],
                           function_names=[lambda x: "N", lambda x: "N^2", lambda x: "sqrt(N)", lambda x: "N^1.5"])
    sindy_path = []
    for thr in (0.01, 0.05, 0.1, 0.2, 0.5):
        mdl = ps.SINDy(feature_library=lib, optimizer=ps.STLSQ(threshold=thr, alpha=1e-6),
                       differentiation_method=ps.SmoothedFiniteDifference(smoother_kws={"window_length": 7,
                                                                                         "polyorder": 2}))
        mdl.fit(Nn[:, None], t=t / (t[-1] - t[0]))
        sindy_path.append(dict(threshold=thr, equation=mdl.equations(precision=3)[0]))
    out[name] = dict(n_points=len(t), N_max=float(N.max()), models=res, sindy=sindy_path,
                     best_by_bic=min(res, key=lambda k: res[k]["bic"]))
    print(name, "best:", out[name]["best_by_bic"],
          {k: (round(v["dBIC"], 1), round(v["r2"], 3), round(v["cv_rmse_N"], 2)) for k, v in res.items()}, flush=True)
    print("   sindy:", [s["equation"] for s in sindy_path], flush=True)
    if si in (2, 5, 6):
        ax = axs[[2, 5, 6].index(si)]
        ax.plot(t, N, "ko", ms=3, label="veri")
        for mname, f in MODELS.items():
            v = res[mname]
            ax.plot(t, integrate(f, (v["r"], v["K"], v["N0"]), t), lw=1, label=f"{mname} ΔBIC={v['dBIC']:.1f}")
        ax.set_title(name, fontsize=8)
        ax.set_xlabel("t")
        ax.set_ylabel("N")
        ax.legend(fontsize=6)
C.save_json("s2b_gs_count.json", out)
fig.tight_layout()
fig.savefig(C.OUT / "fig_s2b_gs_count.png", dpi=110)
