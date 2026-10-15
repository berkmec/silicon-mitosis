"""Exp12 step 5: does the near-threshold floating-body latch form CURRENT FILAMENTS along the channel width?
1D along the width y (cells of dy um), reduced 4 K model of reduced_model.py (BSIM4-anchored), per um of width:
    c_b dV_i/dt = kappa (V_{i+1} - 2 V_i + V_{i-1}) / dy^2 + j_in(V_i; V_ds, V_ov0 + d_i) - j_out(V_i)
    V_ds = V_DD - R_L * sum_i j_d(V_i) dy          (load line: global coupling)
d_i: V_T mismatch per cell from the SKY130 Pelgrom coefficient (vth0_slope = 3.356e-3 V um, units assumed V um):
sigma = 3.356e-3 / sqrt(dy * L). kappa = L_body / R_sq (lateral body conductance; frozen p-well sheet resistance R_sq is
UNKNOWN at 4 K and is scanned). Quasi-static sweep of V_DD up then down with backward-Euler relaxation at each step.
Edges: no-flux. Output: domains of high body potential (V_b > 0.5 V) = current filaments: number, widths, hysteresis."""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve

sys.path.insert(0, str(Path(__file__).resolve().parent))
import reduced_model as rm  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "exp12"
L_CH, L_BODY = 0.15, 0.5  # um


def j_d(vb, vov0, ss, j_T=1e-7):
    nphi = ss / np.log(10)
    vov = vov0 + rm.dvt_of(vb)
    return np.where(vov < 0, j_T * np.exp(vov / nphi), j_T * (1 + vov / nphi))


def rail_voltage(jd, vc, r_D, dy, contacts="both", r_v=None):
    """Drain-rail node voltages at the cell centres (exp13 local inhibitor). Cell i draws jd_i*dy (A) from its node;
    neighbouring nodes are joined by r_D*dy (Ohm; r_D in Ohm per um of rail).
    contacts = "both"/"left": contact(s) at the end(s) held at vc through half a segment.
    contacts = "distributed": every node is tied to the common node vc through its own ballast r_v/dy (r_v in Ohm um),
    no end contacts; the inhibition (screening) length is lambda = sqrt(r_v / r_D) and is position independent.
    KCL: G V = -jd dy + contact terms. r_D = 0 -> V = vc everywhere (exp12 behaviour)."""
    n = len(jd)
    if r_D <= 0:
        return np.full(n, vc)
    g, gc = 1.0 / (r_D * dy), 2.0 / (r_D * dy)
    main = np.full(n, 2.0 * g)
    b = -jd * dy
    if contacts == "distributed":
        gv = dy / r_v
        main = main + gv
        main[0] -= g
        main[-1] -= g
        b = b + gv * vc
        A = diags([-g * np.ones(n - 1), main, -g * np.ones(n - 1)], [-1, 0, 1])
        return spsolve(A.tocsc(), b)
    if n == 1:
        return np.array([vc - jd[0] * dy / (2 * gc if contacts == "both" else gc)])
    main[0] = g + gc
    b[0] += gc * vc
    if contacts == "both":
        main[-1] = g + gc
        b[-1] += gc * vc
    else:
        main[-1] = g
    A = diags([-g * np.ones(n - 1), main, -g * np.ones(n - 1)], [-1, 0, 1])
    return spsolve(A.tocsc(), b)


def filaments(high):
    edges = np.diff(np.concatenate([[0], high.astype(int), [0]]))
    return list(zip(np.where(edges == 1)[0].tolist(), np.where(edges == -1)[0].tolist()))


def lineage(prev, cur):
    """Overlap bookkeeping between consecutive steps: births, deaths, splits (one -> >=2), merges (>=2 -> one)."""
    ov = [[max(0, min(e1, e2) - max(s1, s2)) > 0 for (s2, e2) in cur] for (s1, e1) in prev]
    births = sum(1 for j in range(len(cur)) if not any(ov[i][j] for i in range(len(prev))))
    deaths = sum(1 for i in range(len(prev)) if not any(ov[i]))
    splits = sum(1 for i in range(len(prev)) if sum(ov[i]) >= 2)
    merges = sum(1 for j in range(len(cur)) if sum(ov[i][j] for i in range(len(prev))) >= 2)
    return dict(births=births, deaths=deaths, splits=splits, merges=merges)


def simulate(W=100.0, dy=0.5, R_sq=1e9, R_L=5e3, ss=0.02, vs=0.005, fii=10, vov0=-0.1, g_sub=1e-13, von=1.05,
             vdd_up=(0.6, 3.0), dvdd=0.02, mismatch=True, seed=0, c_b=1e-15, t_relax=1.0, n_sub=12, store_at=(),
             r_D=0.0, contacts="both", r_v=None, down=True, sigma_scale=1.0):
    """down=False: up sweep only; sigma_scale scales the Pelgrom mismatch (exp17)."""
    n = int(round(W / dy))
    rng = np.random.default_rng(seed)
    sig = sigma_scale * 3.356e-3 / np.sqrt(dy * L_CH) if mismatch else 0.0
    d = rng.normal(0.0, sig, n)
    kappa = L_BODY / R_sq
    if n == 1:
        lap = diags([np.zeros(1)], [0]).tocsr()  # a single cell has no lateral flux
    else:
        lap = diags([np.ones(n - 1), -2 * np.ones(n), np.ones(n - 1)], [-1, 0, 1]).tolil()
        lap[0, 0] = lap[-1, -1] = -1.0  # no-flux edges
        lap = (kappa / dy**2) * lap.tocsr()
    nphi = ss / np.log(10)

    def rhs_and_jac(V, vds):
        vov = vov0 + d + rm.dvt_of(V)
        jd = np.where(vov < 0, 1e-7 * np.exp(vov / nphi), 1e-7 * (1 + vov / nphi))
        M = rm.M_of(vds)
        ddvt = np.gradient(rm.AN["dvt"], rm.AN["VB"])
        dvt_p = np.interp(np.clip(V, 0, 1.0), rm.AN["VB"], ddvt)
        djd = np.where(vov < 0, jd / nphi, 1e-7 / nphi) * dvt_p
        jo = 1e-6 * np.exp((V - von) / vs)
        F = lap @ V + fii * M * jd - jo - g_sub * V
        dF = fii * M * djd - jo / vs - g_sub
        return F, dF, jd

    V = np.zeros(n)
    VD = None
    hist, profiles = [], {}
    prev_fil = []
    events = {"up": dict(births=0, deaths=0, splits=0, merges=0), "down": dict(births=0, deaths=0, splits=0, merges=0)}
    vdd_list = list(np.round(np.arange(vdd_up[0], vdd_up[1] + 1e-9, dvdd), 4))
    sweep = [("up", v) for v in vdd_list] + ([("down", v) for v in vdd_list[::-1]] if down else [])
    for direction, vdd in sweep:
        dt = t_relax / n_sub
        for _ in range(n_sub):  # backward Euler, Newton per step (V_ds frozen within a Newton iteration)
            for _it in range(30):
                jd_now = j_d(V, vov0 + d, ss)
                vc = max(vdd - R_L * np.sum(jd_now) * dy, 0.0)
                vd_new = rail_voltage(jd_now, vc, r_D, dy, contacts, r_v)
                VD = vd_new if (VD is None or r_D <= 0) else 0.5 * VD + 0.5 * vd_new  # under-relaxed rail update
                vds = np.maximum(VD, 0.0) if r_D > 0 else vc
                F, dF, _ = rhs_and_jac(V, vds)
                # pseudo-time step toward the stable steady state. In the runaway (positive-feedback) region dF > 0 and
                # c_b/dt - dF < 0 would reverse the step; limit the effective step to the growth time: A_ii >= dF > 0.
                A = diags(np.maximum(c_b / dt, 2.0 * dF) - dF) - lap
                dV = spsolve(A.tocsc(), F)
                dV = np.clip(dV, -0.05, 0.05)
                V = np.clip(V + dV, 0.0, 1.3)
                if np.max(np.abs(dV)) < 1e-7:
                    break
        jd = j_d(V, vov0 + d, ss)
        I = float(np.sum(jd) * dy)
        vds = max(vdd - R_L * I, 0.0)
        high = V > 0.5
        fil = filaments(high)
        widths = [float((e - s) * dy) for s, e in fil]
        ev = lineage(prev_fil, fil)
        for k_ in ev:
            events[direction][k_] += ev[k_]
        prev_fil = fil
        vd_prof = VD if (VD is not None and r_D > 0) else np.full(n, vds)
        hist.append(dict(dir=direction, vdd=float(vdd), vds=float(vds), I_uA=I * 1e6, n_fil=len(widths),
                         widths_um=widths, centres_um=[float((s + e) / 2 * dy) for s, e in fil],
                         frac_high=float(high.mean()), vd_min=float(np.min(vd_prof)), vd_max=float(np.max(vd_prof)),
                         **ev))
        if any(direction == d_ and abs(vdd - v_) < 1e-9 for d_, v_ in store_at):
            profiles[f"{direction}_{vdd:.2f}"] = V.tolist()
            profiles[f"{direction}_{vdd:.2f}_VD"] = vd_prof.tolist()
    return dict(params=dict(W=W, dy=dy, R_sq=R_sq, R_L=R_L, ss=ss, vs=vs, fii=fii, vov0=vov0, g_sub=g_sub, von=von,
                            sigma_vt_cell_mV=sig * 1e3, kappa=kappa, mismatch=mismatch, r_D=r_D, contacts=contacts,
                            r_v=r_v, lambda_um=(float(np.sqrt(r_v / r_D)) if (r_v and r_D) else None)),
                hist=hist, V_final=V.tolist(), profiles=profiles, events=events)


if __name__ == "__main__":
    res = {}
    for R_sq in (1e6, 1e8, 1e10, 1e12):
        for mm in (True, False):
            key = f"Rsq{R_sq:.0e}_{'mismatch' if mm else 'uniform'}"
            r = simulate(R_sq=R_sq, mismatch=mm)
            res[key] = r
            up = [h for h in r["hist"] if h["dir"] == "up"]
            dn = [h for h in r["hist"] if h["dir"] == "down"]
            nmax = max(h["n_fil"] for h in r["hist"])
            first = next((h for h in up if h["n_fil"] > 0), None)
            last_dn = next((h for h in dn[::-1] if h["n_fil"] > 0), None)
            print(f"{key:24s} maks filaman={nmax}  ilk ateşleme V_DD={first['vdd'] if first else None}  "
                  f"aşağıda son V_DD={last_dn['vdd'] if last_dn else None}  "
                  f"V_DD=2,0'da (yukarı) N={next(h['n_fil'] for h in up if abs(h['vdd']-2.0)<1e-9)} "
                  f"genişlik={next(h['widths_um'] for h in up if abs(h['vdd']-2.0)<1e-9)}", flush=True)
    (OUT / "spatial.json").write_text(json.dumps(res), encoding="utf-8")
