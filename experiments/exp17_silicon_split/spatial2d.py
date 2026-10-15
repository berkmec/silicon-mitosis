"""Exp17: 2D version of the exp12/13 filament model - an array of nFET finger segments (cells dx x dx; one finger per
row, segment width dx, L = 0.15 um) sharing a frozen p-well body (sheet R_sq) and a 2D drain sheet (r_sq_D, Ohm/sq) that
is tied to the common drain node V_c through a distributed vertical ballast r_A (Ohm um^2):
    c_b dV/dt = kappa lap V + F(V; V_D),   F = fii M(V_D) j_d(V) - j_diode(V) - g_sub V      (exp12 reduced 4 K model)
    -lap u + u / lambda^2 = r_sq_D j_d / dx,  u = V_c - V_D,  lambda^2 = r_A / r_sq_D        (2D screened inhibitor)
    V_c = V_DD - R_L I,  I = sum j_d dx                                                         (load line)
No-flux edges. Quasi-static up sweep of V_DD; pseudo-time Newton as in spatial.py. Output per step: spots (connected
V > 0.5 regions), their areas / equivalent radii / holes (ring = enclosed low region), V_D range; snapshots."""
import numpy as np
from scipy import ndimage
from scipy.sparse import diags, identity, kron
from scipy.sparse.linalg import splu, spsolve

import law  # noqa: F401  (puts exp12 on the path)
import reduced_model as rm  # noqa: E402

L_CH, L_BODY = 0.15, 0.5


def lap2d(nx, ny, h, periodic=False):
    def lap1(n):
        m = diags([np.ones(n - 1), -2 * np.ones(n), np.ones(n - 1)], [-1, 0, 1]).tolil()
        if periodic:
            m[0, -1] = m[-1, 0] = 1.0
        else:
            m[0, 0] = m[-1, -1] = -1.0
        return m.tocsr()
    return (kron(lap1(nx), identity(ny)) + kron(identity(nx), lap1(ny))).tocsr() / h**2


def spots(high, dx):
    lab, n = ndimage.label(high)
    out = []
    for i in range(1, n + 1):
        m = lab == i
        area = float(m.sum()) * dx**2
        filled = ndimage.binary_fill_holes(m)
        holes = int(ndimage.label(filled & ~m)[1])
        com = ndimage.center_of_mass(m)
        out.append(dict(area=area, R_eq=float(np.sqrt(area / np.pi)), holes=holes,
                        cx=float(com[0] * dx), cy=float(com[1] * dx)))
    return out


def simulate2d_sm(Lx=400.0, Ly=400.0, dx=4.0, R_sq=5e10, R_L=1e3, r_sqD=10.0, lam=60.0, vdd=(0.6, 2.0), dvdd=0.01,
                  seed=0, sigma_scale=1.0, c_b=1e-15, dt=0.125, max_it=600, tol=1e-6, snap_every=0.05, ss=0.02, vs=0.005,
                  fii=10, vov0=-0.1, g_sub=1e-13, von=1.05, verbose=False, periodic=False, seed_dvt=0.0, seed_r=0.0,
                  V0=None):
    """V0: initial body potentials (resume from a checkpoint); the final state is returned as V_final.
    Same model as simulate2d, but the load-line feedback (every V_D depends on the total current through V_c) is put
    into the pseudo-time Newton matrix exactly as a rank-1 term (Sherman-Morrison), V_D is recomputed exactly from V at
    every iteration (no under-relaxation) and each V_DD step is iterated to convergence (max |dV| < tol, flagged).
    periodic: periodic body and drain sheet (no walls). seed_dvt, seed_r: designed nucleation site at the centre
    (V_T lowered by seed_dvt inside radius seed_r)."""
    nx, ny = int(round(Lx / dx)), int(round(Ly / dx))
    N = nx * ny
    rng = np.random.default_rng(seed)
    d = rng.normal(0.0, sigma_scale * 3.356e-3 / np.sqrt(dx * L_CH), N)
    if seed_dvt and seed_r:
        X, Y = np.meshgrid((np.arange(nx) + 0.5) * dx - Lx / 2, (np.arange(ny) + 0.5) * dx - Ly / 2, indexing="ij")
        d += np.where((X**2 + Y**2).ravel() < seed_r**2, seed_dvt, 0.0)       # vov = vov0 + d: + raises overdrive
    body = (L_BODY / R_sq) * lap2d(nx, ny, dx, periodic)
    drain = splu((-lap2d(nx, ny, dx, periodic) + identity(N) / lam**2).tocsc())
    nphi = ss / np.log(10)
    ddvt = np.gradient(rm.AN["dvt"], rm.AN["VB"])
    vgrid = np.linspace(0.0, 3.0, 3001)
    dM = np.gradient(rm.M_of(vgrid), vgrid)

    def state(V, v):
        vov = vov0 + d + rm.dvt_of(V)
        jd = np.where(vov < 0, 1e-7 * np.exp(vov / nphi), 1e-7 * (1 + vov / nphi))
        I = float(jd.sum() * dx)
        vc = v - R_L * I
        VD = vc - drain.solve(r_sqD * jd / dx)
        return jd, vov, I, vc, VD

    V = np.zeros(N) if V0 is None else np.asarray(V0, float).ravel().copy()
    hist, snaps, nxt = [], {}, vdd[0]
    for v in np.round(np.arange(vdd[0], vdd[1] + 1e-9, dvdd), 4):
        conv = False
        for it in range(max_it):
            jd, vov, I, vc, VD = state(V, v)
            vds = np.maximum(VD, 0.0)
            M = rm.M_of(vds)
            dvt_p = np.interp(np.clip(V, 0, 1.0), rm.AN["VB"], ddvt)
            djd = np.where(vov < 0, jd / nphi, 1e-7 / nphi) * dvt_p
            jo = 1e-6 * np.exp((V - von) / vs)
            F = body @ V + fii * M * jd - jo - g_sub * V
            dF = fii * M * djd - jo / vs - g_sub
            A = (diags(np.maximum(c_b / dt, 2.0 * dF) - dF) - body).tocsc()
            a = fii * np.interp(vds, vgrid, dM) * jd            # dF_i / dV_D,i
            b = -R_L * dx * djd                                   # dV_D,i / dV_j through V_c (rank 1)
            lu = splu(A)
            x1, x2 = lu.solve(F), lu.solve(a)
            dV = x1 + x2 * (b @ x1) / (1.0 - b @ x2)              # (A - a b^T)^-1 F
            dV = np.clip(dV, -0.05, 0.05)
            V = np.clip(V + dV, 0.0, 1.3)
            if np.max(np.abs(dV)) < tol:
                conv = True
                break
        jd, _, I, vc, VD = state(V, v)
        hi = (V > 0.5).reshape(nx, ny)
        sp_ = spots(hi, dx)
        hist.append(dict(vdd=float(v), I_uA=I * 1e6, vc=float(vc), n=len(sp_), spots=sp_, conv=conv, it=it + 1,
                         vd_min=float(VD.min()), vd_max=float(VD.max())))
        if v >= nxt - 1e-9:
            snaps[f"V_{v:.2f}"] = V.reshape(nx, ny).astype(np.float16)
            snaps[f"VD_{v:.2f}"] = VD.reshape(nx, ny).astype(np.float32)
            nxt = v + snap_every
        if verbose:
            print(f"V_DD={v:.2f} I={I*1e3:.2f} mA vc={vc:.3f} conv={conv} it={it+1} "
                  f"spots={[(round(s['R_eq'],1), s['holes']) for s in sp_]} VD {VD.min():.4f}-{VD.max():.4f}", flush=True)
    return dict(params=dict(Lx=Lx, Ly=Ly, dx=dx, R_sq=R_sq, R_L=R_L, r_sqD=r_sqD, lam=lam, seed=seed,
                            sigma_scale=sigma_scale, solver="sherman-morrison"), hist=hist, snaps=snaps, V_final=V)


def simulate2d(Lx=400.0, Ly=400.0, dx=4.0, R_sq=5e10, R_L=1e3, r_sqD=10.0, lam=60.0, vdd=(0.6, 2.0), dvdd=0.01,
               seed=0, sigma_scale=1.0, c_b=1e-15, t_relax=1.0, n_sub=8, snap_every=0.05, ss=0.02, vs=0.005, fii=10,
               vov0=-0.1, g_sub=1e-13, von=1.05, verbose=False):
    nx, ny = int(round(Lx / dx)), int(round(Ly / dx))
    N = nx * ny
    rng = np.random.default_rng(seed)
    d = rng.normal(0.0, sigma_scale * 3.356e-3 / np.sqrt(dx * L_CH), N)
    kappa = L_BODY / R_sq
    Lh = lap2d(nx, ny, dx)
    body = kappa * Lh
    drain = splu((-Lh + identity(N) / lam**2).tocsc())
    nphi = ss / np.log(10)
    ddvt = np.gradient(rm.AN["dvt"], rm.AN["VB"])

    def jd_of(V):
        vov = vov0 + d + rm.dvt_of(V)
        return np.where(vov < 0, 1e-7 * np.exp(vov / nphi), 1e-7 * (1 + vov / nphi)), vov

    V = np.zeros(N)
    VD = None
    hist, snaps, nxt = [], {}, vdd[0]
    for v in np.round(np.arange(vdd[0], vdd[1] + 1e-9, dvdd), 4):
        dt = t_relax / n_sub
        for _ in range(n_sub):
            for _it in range(30):
                jd, vov = jd_of(V)
                I = float(jd.sum() * dx)
                vc = max(v - R_L * I, 0.0)
                vd_new = vc - drain.solve(r_sqD * jd / dx)
                VD = vd_new if VD is None else 0.5 * VD + 0.5 * vd_new
                vds = np.maximum(VD, 0.0)
                M = rm.M_of(vds)
                dvt_p = np.interp(np.clip(V, 0, 1.0), rm.AN["VB"], ddvt)
                djd = np.where(vov < 0, jd / nphi, 1e-7 / nphi) * dvt_p
                jo = 1e-6 * np.exp((V - von) / vs)
                F = body @ V + fii * M * jd - jo - g_sub * V
                dF = fii * M * djd - jo / vs - g_sub
                A = diags(np.maximum(c_b / dt, 2.0 * dF) - dF) - body
                dV = np.clip(spsolve(A.tocsc(), F), -0.05, 0.05)
                V = np.clip(V + dV, 0.0, 1.3)
                if np.max(np.abs(dV)) < 1e-7:
                    break
        jd, _ = jd_of(V)
        I = float(jd.sum() * dx)
        hi = (V > 0.5).reshape(nx, ny)
        sp_ = spots(hi, dx)
        hist.append(dict(vdd=float(v), I_uA=I * 1e6, vc=float(max(v - R_L * I, 0.0)), n=len(sp_), spots=sp_,
                         vd_min=float(VD.min()), vd_max=float(VD.max())))
        if v >= nxt - 1e-9:
            snaps[f"V_{v:.2f}"] = V.reshape(nx, ny).astype(np.float16)
            snaps[f"VD_{v:.2f}"] = VD.reshape(nx, ny).astype(np.float32)
            nxt = v + snap_every
        if verbose:
            print(f"V_DD={v:.2f} I={I*1e6:.1f} uA spots={[(round(s['R_eq'],1), s['holes']) for s in sp_]}", flush=True)
    return dict(params=dict(Lx=Lx, Ly=Ly, dx=dx, R_sq=R_sq, R_L=R_L, r_sqD=r_sqD, lam=lam, seed=seed,
                            sigma_scale=sigma_scale, kappa=kappa), hist=hist, snaps=snaps)
