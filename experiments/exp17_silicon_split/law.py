"""Exp17: split law of a current filament in a wide SKY130 nFET at 4 K with a distributed drain ballast (exp13 geometry).

Rail: drop u = V_c - V_D obeys u'' = u / lambda^2 - r_D j_d (lambda^2 = r_v / r_D). For a filament of width w carrying j_on
per um, the drop is largest at the centre and the centre-to-edge difference is
    Delta_u(w) = (r_v j_on / 2) (1 - exp(-w / (2 lambda)))^2 .
The fronts (edges) sit at the Maxwell voltage V_M; the centre loses its high state when its V_D falls below V_down
(saddle-node of the local bistability). Split when Delta_u = V_M - V_down = dV_c:
    lambda_c = sqrt(2 dV_c / (r_D j_on)),   w_split = 2 lambda ln(lambda / (lambda - lambda_c))   (lambda > lambda_c)
    (minimum 2 lambda_c for lambda >> lambda_c; no split for lambda <= lambda_c).
Total filament width (exp16 law, here derived): w_tot = I / j_on. The constants come from the reduced 4 K model of exp12
(constants()); nothing is fitted."""
import json
import sys
from pathlib import Path

import numpy as np
from scipy import integrate, optimize

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "experiments" / "exp12_filament"))
import reduced_model as rm  # noqa: E402

BASE = dict(ss=0.02, vs=0.005, fii=10, vov0=-0.1, g_sub=1e-13, von=1.05)


def _jd(V, p):
    nphi = p["ss"] / np.log(10)
    vov = p["vov0"] + rm.dvt_of(V)
    return np.where(vov < 0, 1e-7 * np.exp(vov / nphi), 1e-7 * (1 + vov / nphi))


def _F(V, vds, p):
    return p["fii"] * rm.M_of(vds) * _jd(V, p) - 1e-6 * np.exp((V - p["von"]) / p["vs"]) - p["g_sub"] * V


def constants(p=BASE):
    """V_down (lowest V_ds with a high state), V_M (equal-area front condition), j_on (high-state current at V_M)."""
    Vg = np.linspace(0.0, 1.3, 26001)

    def roots(vds):
        f = _F(Vg, vds, p)
        idx = np.nonzero(np.sign(f[:-1]) != np.sign(f[1:]))[0]
        return [optimize.brentq(_F, Vg[i], Vg[i + 1], args=(vds, p)) for i in idx]

    has_high = lambda v: bool(roots(v)) and max(roots(v)) > 0.5
    has_low = lambda v: bool(roots(v)) and min(roots(v)) < 0.5
    a, b = 0.3, 3.0
    for _ in range(50):
        m = 0.5 * (a + b)
        a, b = (a, m) if has_high(m) else (m, b)
    v_down = b
    a, b = v_down, 3.0                          # V_up: largest V_ds with a low state (end of the bistable window)
    for _ in range(50):
        m = 0.5 * (a + b)
        a, b = (m, b) if has_low(m) else (a, m)
    v_up = a
    area = lambda v: integrate.quad(lambda V: _F(V, v, p), min(roots(v)), max(roots(v)), limit=200)[0]
    lo, hi = v_down + 1e-5, v_up - 1e-5
    if not (area(lo) < 0 < area(hi)):
        # no equal-area point inside the window: fronts always advance (high state invades) until V_D reaches V_down
        return dict(V_down=v_down, V_up=v_up, V_M=float("nan"), dV_c=float("nan"),
                    j_on=float(_jd(max(roots(lo)), p)), regime="no-Maxwell")
    v_m = optimize.brentq(area, lo, hi)
    return dict(V_down=v_down, V_up=v_up, V_M=v_m, dV_c=v_m - v_down, j_on=float(_jd(max(roots(v_m)), p)),
                regime="Maxwell")


def lambda_c(r_D, c):
    return float(np.sqrt(2 * c["dV_c"] / (r_D * c["j_on"])))


def w_split(lam, r_D, c):
    lc = lambda_c(r_D, c)
    return float(2 * lam * np.log(lam / (lam - lc))) if lam > lc else float("inf")


def du_array(w, P, lam, r_D, c):
    """Centre-to-edge drop of a filament in a periodic array (period P; P -> inf gives the isolated result):
    Delta_u = r_v j_on (cosh a - 1) / (cosh a + sinh a coth b), a = w / 2 lam, b = (P - w) / 2 lam."""
    a, b = w / (2 * lam), (P - w) / (2 * lam)
    return float(r_D * lam**2 * c["j_on"] * (np.cosh(a) - 1) / (np.cosh(a) + np.sinh(a) / np.tanh(b)))


def w_split_array(lam, r_D, P, c):
    g = lambda w: du_array(w, P, lam, r_D, c) - c["dV_c"]
    ws = np.linspace(1e-3, P * (1 - 1e-6), 4000)
    v = np.array([g(w) for w in ws])
    i = np.nonzero((v[:-1] < 0) & (v[1:] >= 0))[0]
    return float(optimize.brentq(g, ws[i[0]], ws[i[0] + 1])) if len(i) else float("inf")


def du_config(intervals, parent, W, lam, r_D, c, dy=0.5):
    """Exact (discrete rail, finite device with reflecting ends) centre-to-edge drop of the filament `parent` for a given
    configuration of filaments `intervals` [(start, end) in um], all carrying j_on."""
    import spatial as sp  # exp12 rail solver
    n = int(round(W / dy))
    y = (np.arange(n) + 0.5) * dy
    jd = np.zeros(n)
    for s, e in intervals:
        jd[(y >= s) & (y < e)] = c["j_on"]
    V = sp.rail_voltage(jd, 1.0, r_D, dy, "distributed", r_D * lam**2)
    s, e = parent
    edge = 0.5 * (np.interp(s, y, V) + np.interp(e, y, V))
    return float(edge - np.interp(0.5 * (s + e), y, V))


def vdd_at_split(lam, r_D, R_L, c):
    """V_DD = V_c + R_L I at the split; V_c = V_M + u(edge), u(edge) = r_v j_on (1 - exp(-2a)) / 2, a = w / 2 lambda."""
    w = w_split(lam, r_D, c)
    if not np.isfinite(w):
        return float("inf")
    r_v = r_D * lam**2
    a = w / (2 * lam)
    return float(c["V_M"] + r_v * c["j_on"] * (1 - np.exp(-2 * a)) / 2 + R_L * c["j_on"] * w)


def j_high_curve(p=BASE, c=None, n=400):
    """High-state current density j_high(V_D) for V_D from V_down (fold) to 1.2 V (reduced model, no mismatch)."""
    c = constants(p) if c is None else c
    Vg = np.linspace(0.0, 1.3, 26001)
    vd = np.linspace(c["V_down"] + 1e-6, 1.2, n)
    jh = []
    for v in vd:
        f = _F(Vg, v, p)
        idx = np.nonzero(np.sign(f[:-1]) != np.sign(f[1:]))[0]
        r = [optimize.brentq(_F, Vg[i], Vg[i + 1], args=(v, p)) for i in idx]
        jh.append(float(_jd(max(r), p)))
    return vd, np.array(jh)


def w_split_nl(lam, r_D, c, curve=None, wmax=2000.0):
    """Split width with the self-consistent current density j_high(V_D) inside the filament (zero outside). Inside:
    u'' = u / lam^2 - r_D j_high(V_c - u); outside u = u_e exp(-(y - w/2) / lam); edge V_D = V_M; split when the centre
    reaches V_down. Shooting from the edge (u = u_e, u' = -u_e / lam) to the centre (u' = 0)."""
    vd, jh = curve if curve is not None else j_high_curve(c=c)
    jfun = lambda v: np.interp(v, vd, jh, left=jh[0], right=jh[-1])

    def centre(w):
        def shoot(ue, full=False):
            vc = c["V_M"] + ue
            f = lambda y, s: [s[1], s[0] / lam**2 - r_D * jfun(vc - s[0])]
            sol = integrate.solve_ivp(f, [w / 2, 0.0], [ue, -ue / lam], rtol=1e-10, atol=1e-14)
            return sol.y[1, -1] if not full else sol.y[0, -1]
        lo, hi = 1e-9, r_D * lam**2 * c["j_on"] * 2
        ue = optimize.brentq(shoot, lo, hi, xtol=1e-14)
        return c["V_M"] + ue - shoot(ue, True)          # V_D at the centre

    g = lambda w: centre(w) - c["V_down"]
    ws = np.geomspace(5.0, wmax, 60)
    vals = [g(w) for w in ws]
    for i in range(len(ws) - 1):
        if vals[i] > 0 >= vals[i + 1]:
            return float(optimize.brentq(g, ws[i], ws[i + 1], xtol=1e-3))
    return float("inf")


if __name__ == "__main__":
    c = constants()
    out = dict(constants=c, lambda_c={str(r): lambda_c(r, c) for r in (5.0, 10.0, 20.0)}, predictions={})
    for r_D in (5.0, 10.0, 20.0):
        for lam in (20, 30, 36, 40, 45, 50, 60, 80, 120):
            out["predictions"][f"rD{r_D:g}_lam{lam}"] = dict(w_split=w_split(lam, r_D, c),
                                                            vdd_RL2k=vdd_at_split(lam, r_D, 2e3, c),
                                                            vdd_RL5k=vdd_at_split(lam, r_D, 5e3, c))
    print(json.dumps(out, indent=1))
    (ROOT / "results" / "exp17").mkdir(parents=True, exist_ok=True)
    dst = ROOT / "results" / "exp17" / "predictions.json"
    if dst.exists():   # the pre-registered record: never overwrite it
        dst = dst.with_name("predictions_recomputed.json")
    dst.write_text(json.dumps(out, indent=1), encoding="utf-8")
    print("written:", dst)
