"""Exp16 shared helpers: process priority/threads, paths, scaled STLSQ (PySINDy) threshold paths, bootstrap, metrics.

Import this module FIRST in every exp16 script (it limits BLAS threads before numpy is imported).
Sparse regression is always PySINDy's STLSQ (https://github.com/dynamicslab/pysindy). Columns of the library and the
target are scaled to unit RMS before the fit, so the threshold tau is dimensionless and independent of sample size:
a term is dropped when |xi_j| * rms(theta_j) < tau * rms(y)."""
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "2")

import ctypes  # noqa: E402
import json  # noqa: E402
import sys  # noqa: E402
import warnings  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results"
OUT = RESULTS / "exp16_equations"
OUT.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(ROOT / "src"))
warnings.filterwarnings("ignore", category=UserWarning)


def low_priority():
    """BELOW_NORMAL_PRIORITY_CLASS for this process (keeps the desktop and other jobs responsive)."""
    if sys.platform != "win32":
        return
    from ctypes import wintypes

    k = ctypes.windll.kernel32
    k.GetCurrentProcess.restype = wintypes.HANDLE
    k.SetPriorityClass.argtypes = [wintypes.HANDLE, wintypes.DWORD]
    k.SetPriorityClass(k.GetCurrentProcess(), 0x4000)


def _clean(o):
    if isinstance(o, dict):
        return {str(k): _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_clean(v) for v in o]
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, np.ndarray):
        return _clean(o.tolist())
    if isinstance(o, float) and not np.isfinite(o):
        return None
    return o


def save_json(name, obj):
    p = OUT / name
    p.write_text(json.dumps(_clean(obj), indent=1, ensure_ascii=False), encoding="utf-8")
    return p


def r2(y, yhat):
    y, yhat = np.asarray(y, float).ravel(), np.asarray(yhat, float).ravel()
    ss = np.sum((y - y.mean()) ** 2)
    return float(1.0 - np.sum((y - yhat) ** 2) / ss) if ss > 0 else float("nan")


def rel_rmse(y, yhat):
    y, yhat = np.asarray(y, float).ravel(), np.asarray(yhat, float).ravel()
    return float(np.sqrt(np.mean((y - yhat) ** 2)) / np.sqrt(np.mean(y**2)))


def stlsq_scaled(Theta, y, tau, alpha=1e-10, max_iter=50):
    """PySINDy STLSQ on RMS-scaled columns. Theta (N, p), y (N,) or (N, m). Returns coefficients (m, p) in ORIGINAL
    units (zeros for dropped terms)."""
    import pysindy as ps

    Theta = np.asarray(Theta, float)
    Y = np.asarray(y, float)
    Y = Y[:, None] if Y.ndim == 1 else Y
    cs = np.sqrt(np.mean(Theta**2, axis=0))
    cs[cs == 0] = 1.0
    ys = np.sqrt(np.mean(Y**2, axis=0))
    ys[ys == 0] = 1.0
    opt = ps.STLSQ(threshold=tau, alpha=alpha, max_iter=max_iter, normalize_columns=False, unbias=True)
    opt.fit(Theta / cs, Y / ys)
    return opt.coef_ * ys[:, None] / cs[None, :]


def threshold_path(Theta_tr, y_tr, Theta_va, y_va, taus, alpha=1e-10):
    """Fit STLSQ for every tau; record support size, train/validation R^2 (per target) and coefficients."""
    rows = []
    for tau in taus:
        C = stlsq_scaled(Theta_tr, y_tr, tau, alpha=alpha)
        Ytr = np.asarray(y_tr, float).reshape(len(Theta_tr), -1)
        Yva = np.asarray(y_va, float).reshape(len(Theta_va), -1)
        ptr, pva = Theta_tr @ C.T, Theta_va @ C.T
        rows.append(dict(tau=float(tau), n_terms=[int(np.sum(c != 0)) for c in C],
                         r2_train=[r2(Ytr[:, j], ptr[:, j]) for j in range(C.shape[0])],
                         r2_val=[r2(Yva[:, j], pva[:, j]) for j in range(C.shape[0])],
                         rmse_val=[rel_rmse(Yva[:, j], pva[:, j]) for j in range(C.shape[0])],
                         coef=C))
    return rows


def select_model(path, j=0, tol=0.05, min_terms=1):
    """Sparsest model whose validation rel. RMSE (target j) is within (1+tol) of the best one on the path.
    The true model is never consulted (selection uses held-out data only)."""
    best = min(r["rmse_val"][j] for r in path if r["n_terms"][j] >= min_terms)
    ok = [r for r in path if r["n_terms"][j] >= min_terms and r["rmse_val"][j] <= (1 + tol) * best]
    return min(ok, key=lambda r: (r["n_terms"][j], r["rmse_val"][j]))


def select_model_r2(path, j=0, r2_tol=1e-3, min_terms=1):
    """Parsimony rule B: sparsest model whose validation R^2 (target j) is at most r2_tol below the best one.
    Robust to systematic derivative errors that extra terms can partly absorb (rule A then over-selects)."""
    best = max(r["r2_val"][j] for r in path if r["n_terms"][j] >= min_terms)
    ok = [r for r in path if r["n_terms"][j] >= min_terms and r["r2_val"][j] >= best - r2_tol]
    return min(ok, key=lambda r: (r["n_terms"][j], -r["r2_val"][j]))


def bootstrap(Theta, y, tau, n_boot=200, seed=0, alpha=1e-10):
    """Row bootstrap at fixed tau: inclusion probability and mean/std of every coefficient (target 0 if 1-D)."""
    rng = np.random.default_rng(seed)
    N = len(Theta)
    Cs = []
    for _ in range(n_boot):
        idx = rng.integers(0, N, N)
        Cs.append(stlsq_scaled(Theta[idx], np.asarray(y)[idx], tau, alpha=alpha))
    Cs = np.array(Cs)  # (B, m, p)
    return dict(p_incl=(Cs != 0).mean(axis=0), mean=Cs.mean(axis=0), std=Cs.std(axis=0))


def eq_string(names, coef, prec=4):
    terms = [f"{c:+.{prec}g}·{n}" if n != "1" else f"{c:+.{prec}g}" for n, c in zip(names, coef) if c != 0]
    return " ".join(terms) if terms else "0"


def loo_cv(Theta, y, tau, alpha=1e-10):
    """Leave-one-out predictions with re-selection of the support inside every fold (honest CV for tiny data)."""
    N = len(Theta)
    pred = np.zeros(N)
    for i in range(N):
        m = np.ones(N, bool)
        m[i] = False
        C = stlsq_scaled(Theta[m], y[m], tau, alpha=alpha)
        pred[i] = (Theta[i] @ C.T).item()
    return pred
