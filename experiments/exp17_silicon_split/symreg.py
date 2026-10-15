"""Exp17: equation from data for the filament split width (clean protocol runs E/F/G/H, first split of a filament that
was born narrower than the split width).
(1) Model selection by leave-one-out (LOO) relative error among candidate forms:
    P0  physics, no free parameter:   w = 2 lam ln(lam / (lam - lam_c)), lam_c = sqrt(2 dV_c / (r_D j_on))
    P0n physics, self-consistent j_high(V_D) (law.w_split_nl), no free parameter
    P1  physics form, lam_c = k / sqrt(r_D) with one fitted k (compare k with sqrt(2 dV_c / j_on))
    Q3  power law w = A lam^a r_D^b
    Q2  w = A + B / (lam sqrt(r_D))      (a smooth 2-parameter alternative)
(2) Unbiased symbolic regression with gplearn (genetic programming; run in env/venv_sindy) on (lam, r_D) -> w.
Writes results/exp17/symreg.json."""
import glob
import json
import sys
from pathlib import Path

import numpy as np
from scipy import optimize

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "exp17"


def load():
    rows = []
    for f in sorted(glob.glob(str(OUT / "runs" / "[EFGH]_*.json"))):
        r = json.loads(Path(f).read_text(encoding="utf-8"))
        h = r["hist"]
        birth = next((x for x in h if x["n_fil"] > 0), None)
        if not r["splits"] or birth is None:
            rows.append(dict(name=Path(f).stem, lam=r["job"]["lam"], r_D=r["job"]["r_D"], w=None,
                             w_birth=max(birth["widths_um"]) if birth else None, mm=r["job"]["sigma_scale"]))
            continue
        s = r["splits"][0]
        rows.append(dict(name=Path(f).stem, lam=r["job"]["lam"], r_D=r["job"]["r_D"], w=s["parent_w"],
                         w_birth=max(birth["widths_um"]), mm=r["job"]["sigma_scale"], vdd=s["vdd"],
                         pred_lin=r["prediction"]["w_split"]))
    return rows


def loo(fit, predict, X, y):
    err = []
    for i in range(len(y)):
        m = np.ones(len(y), bool)
        m[i] = False
        p = fit(X[m], y[m])
        err.append((predict(p, X[i:i + 1])[0] - y[i]) / y[i])
    return float(np.sqrt(np.mean(np.square(err)))), err


def main():
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import law
    c = law.constants()
    rows = load()
    ok = [r for r in rows if r["w"] is not None and r["w_birth"] < 0.8 * r["w"] and r["mm"] == 1.0]
    X = np.array([[r["lam"], r["r_D"]] for r in ok])
    y = np.array([r["w"] for r in ok])
    k_th = np.sqrt(2 * c["dV_c"] / c["j_on"])                       # um * sqrt(Ohm/um)

    def phys(k, X):
        lc = k / np.sqrt(X[:, 1])
        return np.where(X[:, 0] > lc, 2 * X[:, 0] * np.log(X[:, 0] / np.maximum(X[:, 0] - lc, 1e-9)), 1e6)

    curve = law.j_high_curve(c=c)
    res = dict(n=len(ok), data=ok, k_theory=float(k_th), models={})
    p0 = phys(k_th, X)
    res["models"]["P0"] = dict(params=0, rel_rmse=float(np.sqrt(np.mean(((p0 - y) / y) ** 2))), pred=p0.tolist())
    p0n = np.array([law.w_split_nl(l, rd, c, curve=curve) for l, rd in X])
    res["models"]["P0n"] = dict(params=0, rel_rmse=float(np.sqrt(np.mean(((p0n - y) / y) ** 2))), pred=p0n.tolist())
    fitP1 = lambda Xa, ya: optimize.minimize_scalar(lambda k: np.sum(((phys(k, Xa) - ya) / ya) ** 2),
                                                    bounds=(0.3 * k_th, 3 * k_th), method="bounded").x
    kfit = fitP1(X, y)
    res["models"]["P1"] = dict(params=1, k_fit=float(kfit), k_over_theory=float(kfit / k_th),
                               loo=loo(fitP1, lambda k, Xa: phys(k, Xa), X, y)[0])
    fitQ3 = lambda Xa, ya: np.linalg.lstsq(np.c_[np.ones(len(ya)), np.log(Xa)], np.log(ya), rcond=None)[0]
    predQ3 = lambda b, Xa: np.exp(np.c_[np.ones(len(Xa)), np.log(Xa)] @ b)
    b3 = fitQ3(X, y)
    res["models"]["Q3"] = dict(params=3, coef=b3.tolist(), loo=loo(fitQ3, predQ3, X, y)[0])
    fitQ2 = lambda Xa, ya: np.linalg.lstsq(np.c_[np.ones(len(ya)), 1 / (Xa[:, 0] * np.sqrt(Xa[:, 1]))], ya, rcond=None)[0]
    predQ2 = lambda b, Xa: np.c_[np.ones(len(Xa)), 1 / (Xa[:, 0] * np.sqrt(Xa[:, 1]))] @ b
    res["models"]["Q2"] = dict(params=2, coef=fitQ2(X, y).tolist(), loo=loo(fitQ2, predQ2, X, y)[0])
    (OUT / "symreg.json").write_text(json.dumps(res, indent=1, default=float), encoding="utf-8")
    for k, v in res["models"].items():
        print(k, {kk: (round(vv, 4) if isinstance(vv, float) else vv) for kk, vv in v.items() if kk != "pred"})
    print("data:", [(r["name"], r["w"], round(r["pred_lin"], 1)) for r in ok])
    np.save(OUT / "symreg_xy.npy", np.c_[X, y])


if __name__ == "__main__":
    main()
