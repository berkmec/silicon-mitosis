"""Exp15 analiz: her kosu icin yapi sayisi N(t), birlesme/bolunme olaylari, enine uyum, dik isinma, kink.

Cikti: results/exp15_phase_space/analysis.json, analysis_series.npz
"""
import os
os.environ["OMP_NUM_THREADS"] = "1"
import sys
if sys.platform == "win32":  # below-normal priority, Windows only
    import ctypes
    from ctypes import wintypes

    k32 = ctypes.windll.kernel32
    k32.GetCurrentProcess.restype = wintypes.HANDLE
    k32.SetPriorityClass.argtypes = [wintypes.HANDLE, wintypes.DWORD]
    k32.SetPriorityClass(k32.GetCurrentProcess(), 0x4000)
import sys
import json
import glob
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import structures as S

RES = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..",
                                   "results", "exp15_phase_space"))
RUNS = os.path.join(RES, "runs")
FAMILY_REF = {"L40": "1d_L40_seed", "L20": "1d_L20_seed", "L100": "1d_L40_seed"}
FRACS = (0.3, 0.5)


def load(name):
    p = os.path.join(RUNS, name + ".npz")
    if not os.path.exists(p):
        return None
    z = np.load(p)
    meta = json.load(open(os.path.join(RUNS, name + ".json"), encoding="utf-8"))
    return z, meta


def family(name):
    for f in ("L100", "L40", "L20"):
        if f"_{f}" in name:
            return f
    return None


TAVG = 7   # anlik goruntu (= 7/omega_p ~ 1.1 plazma periyodu): Langmuir salinim titresimini bastirmak icin


def tavg(phi):
    """Zaman ekseninde merkezli kayan ortalama (TAVG kare); yapilar yavas, Langmuir salinimi hizli."""
    from scipy.ndimage import uniform_filter1d
    return uniform_filter1d(phi.astype(float), size=TAVG, axis=0, mode="nearest")


def phi_ref(name):
    r = load(name)
    if r is None:
        return None
    z, _ = r
    t, phi = z["t"], tavg(z["phi"])
    sel = (t >= 10) & (t <= 40)
    return float(max(S.smooth(p).max() for p in phi[sel]))


def analyze(name, refs):
    r = load(name)
    if r is None:
        return None, None
    z, meta = r
    t, phi, WE, KE = z["t"], z["phi"], z["WE"], z["KE"]
    L, N = z["L"], z["N"]
    h = L / N
    d = len(N)
    fam = family(name)
    pref = refs[fam]
    min_size = max(2, int(round(np.prod(1.0 / h))))
    out = dict(name=name, d=d, L=L.tolist(), N=N.tolist(), B=meta["config"]["B"], npart=meta["npart"],
               energy_drift_max=meta["energy_drift_max"], energy_drift_end=meta["energy_drift_end"],
               runtime_s=meta["runtime_s"], phi_ref=pref, min_size_cells=min_size)
    ser = dict(t=t, WE=WE)
    # doygunluk
    early = t <= 60
    i_sat = int(np.argmax(np.where(early, WE, -1)))
    out["t_sat"] = float(t[i_sat])
    out["WE_sat"] = float(WE[i_sat])
    out["phi_max_sat"] = float(phi[i_sat].max())
    # esiksiz hayatta kalma olcutleri
    pmax = np.array([S.smooth(p.astype(float)).max() for p in phi])
    ser["phimax"] = pmax
    out["WE_over_WEsat_at"] = {str(int(tt)): float(WE[np.argmin(np.abs(t - tt))] / WE[i_sat])
                               for tt in (30, 50, 75, 100, 150, 200, 300, 600) if tt <= t[-1] + 0.5}
    out["phimax_over_ref_at"] = {str(int(tt)): float(pmax[np.argmin(np.abs(t - tt))] / pref)
                                 for tt in (30, 50, 75, 100, 150, 200, 300, 600) if tt <= t[-1] + 0.5}
    out["KE_x_end_over_0"] = float(KE[-1, 0] / KE[0, 0])
    # yapi sayisi ve olaylar (zaman-ortalamali phi uzerinde)
    phi_t = tavg(phi)
    for fr in FRACS:
        cnt, ev = S.track(phi_t, fr * pref, min_size)
        key = f"th{int(fr * 100)}"
        ser[f"N_{key}"] = cnt
        for kk, vv in ev.items():
            ser[f"{kk}_{key}"] = vv
        after = t >= out["t_sat"] + 5
        n_form = int(cnt[np.argmin(np.abs(t - (out["t_sat"] + 5)))])
        tw = t[1:]
        sel_ev = tw >= out["t_sat"] + 5
        out[key] = dict(theta=fr * pref, N_formed=n_form, N_max_after=int(cnt[after].max()),
                        t_N_max=float(t[after][np.argmax(cnt[after])]), N_end=int(cnt[-1]),
                        splits=int(ev["splits"][sel_ev].sum()), merges=int(ev["merges"][sel_ev].sum()),
                        splits_persistent=int(ev["splits_persistent"][sel_ev].sum()),
                        births=int(ev["births"][sel_ev].sum()), deaths=int(ev["deaths"][sel_ev].sum()),
                        N_at={str(int(tt)): int(cnt[np.argmin(np.abs(t - tt))]) for tt in
                              (20, 30, 50, 75, 100, 150, 200, 300, 500, 1000) if tt <= t[-1] + 0.5})
    # enine uyum ve dik enerji
    if d > 1:
        coh = np.array([S.coherence(p.astype(float), h) for p in phi])
        ser["C"] = coh[:, 0]; ser["Wpar"] = coh[:, 1]; ser["Wperp"] = coh[:, 2]
        out["C_at"] = {str(int(tt)): float(coh[np.argmin(np.abs(t - tt)), 0]) for tt in (20, 30, 50, 75, 100, 150, 200)
                       if tt <= t[-1] + 0.5}
        Kp = KE[:, 1:].sum(axis=1)
        ser["Kperp"] = Kp
        out["Kperp_gain_rel_end"] = float(Kp[-1] / Kp[0] - 1.0)
        out["Kperp_gain_abs_end"] = float(Kp[-1] - Kp[0])
        out["Kperp_gain_frac_of_beam_KE"] = float((Kp[-1] - Kp[0]) / (KE[0].sum()))
        # dik alan enerjisinin buyumesi: Wperp'in doygunluk sonrasi ilk 2 basamakta ustel uyumu
        Wp = coh[:, 2]

        def steepest(t0, t1, win):
            sel = (t >= t0) & (t <= t1) & (Wp > 0)
            if sel.sum() < 6:
                return None
            tt, lw = t[sel], np.log(Wp[sel])
            best = None
            for i0 in range(len(tt)):
                i1 = np.searchsorted(tt, tt[i0] + win)
                if i1 - i0 < 5 or tt[i1 - 1] - tt[i0] < 0.8 * win:
                    continue
                g = np.polyfit(tt[i0:i1], lw[i0:i1], 1)[0] / 2.0
                if best is None or g > best[0]:
                    best = (float(g), float(tt[i0]), float(tt[i1 - 1]))
            return best

        ipk = int(np.argmax(Wp))
        out["t_Wperp_peak"] = float(t[ipk])
        out["Wperp_peak_over_W"] = float(Wp[ipk] / (coh[ipk, 1] + Wp[ipk]))
        kt = meta["config"].get("kick_t")
        t0p = out["t_sat"] + 3.0 if kt is None else max(out["t_sat"] + 3.0, kt + 2.0)
        if kt is not None:
            ipk = int(np.argmax(np.where(t >= kt, Wp, -1)))
            out["t_Wperp_peak"] = float(t[ipk])
            out["Wperp_peak_over_W"] = float(Wp[ipk] / (coh[ipk, 1] + Wp[ipk]))
            out["kick_t"] = kt
        e1 = steepest(0.0, out["t_sat"], 8.0) if kt is None else None
        e2 = steepest(t0p, max(t0p, out["t_Wperp_peak"]), 10.0)
        if kt is not None:
            # saglam uyum: [kick+5, W_perp'in itme sonrasi tepesinin %10'unu ilk astigi an]
            post = t >= kt + 5.0
            Wmax = Wp[t >= kt].max()
            idx_hi = np.where(post & (Wp > 0.1 * Wmax))[0]
            t_hi = float(t[idx_hi[0]]) if idx_hi.size else float(t[-1])
            sel = post & (t <= t_hi) & (Wp > 0)
            if sel.sum() >= 5:
                g = float(np.polyfit(t[sel], np.log(Wp[sel]), 1)[0] / 2.0)
                out["gamma_perp_kick_fit"] = dict(gamma=g, window=[kt + 5.0, t_hi],
                                                  growth_factor_W=float(Wp[sel][-1] / Wp[sel][0]),
                                                  note="W_perp ~ exp(2 gamma t), itme+5 .. W_perp tepesinin %10'u")
        out["gamma_perp_early"] = e1 and dict(gamma=e1[0], window=[e1[1], e1[2]],
                                              note="Wperp ~ exp(2 gamma t); doygunluk oncesi (egik iki-demet modlari)")
        out["gamma_perp_post"] = e2 and dict(gamma=e2[0], window=[e2[1], e2[2]],
                                             note="doygunluktan sonra Wperp tepesine kadar en dik 10-birimlik pencere")
    # 2D: kink genligi (duzlemsel delik merkezleri t_sat+5'te y-ortalamasindan)
    if d == 2:
        i0 = int(np.argmin(np.abs(t - (out["t_sat"] + 5))))
        prof = phi[i0].mean(axis=1)
        X = S.hole_centers_1d(prof, h[0], 0.3 * pref)
        lam = L[0] / max(1, len(X))
        kink = np.array([S.kink_amplitude(p.astype(float), h[0], X, 0.3 * lam) for p in phi])
        ser["kink"] = kink
        out["holes_planar_X"] = X.tolist()
        out["kink_at"] = {str(int(tt)): float(kink[np.argmin(np.abs(t - tt))]) for tt in (20, 30, 50, 75, 100, 150, 200)
                          if tt <= t[-1] + 0.5}
    # dagilim: f(vx) son anda v=0 civarinda cukur derinligi (tuzaklanmis bolgenin izi) ve genislik
    fb = z["fbins"]
    vc = 0.5 * (fb[1:] + fb[:-1])
    fvx = z["fvx"].astype(float)
    if fvx.size:
        f_end = fvx[-1] / fvx[-1].sum()
        out["vx_rms_end"] = float(np.sqrt(np.sum(f_end * vc ** 2)))
        out["ps_times"] = z["ps_t"].tolist()
    return out, ser


def main():
    names = sorted(os.path.basename(p)[:-4] for p in glob.glob(os.path.join(RUNS, "*.npz")))
    refs = {f: phi_ref(n) for f, n in FAMILY_REF.items()}
    print("phi_ref:", refs, flush=True)
    res, series = {}, {}
    for nm in names:
        o, s = analyze(nm, refs)
        if o is None:
            continue
        res[nm] = o
        for k, v in s.items():
            series[f"{nm}__{k}"] = v
        print(nm, {k: o[k] for k in ("t_sat", "phi_max_sat")}, "th30:", o["th30"]["N_formed"], o["th30"]["N_max_after"],
              o["th30"]["N_end"], "spl", o["th30"]["splits"], "spl_pers", o["th30"]["splits_persistent"],
              "mrg", o["th30"]["merges"],
              "dth", o["th30"]["deaths"], "brn", o["th30"]["births"], "| th50:", o["th50"]["N_formed"],
              o["th50"]["N_max_after"], o["th50"]["N_end"], o["th50"]["splits"], o["th50"]["merges"],
              "| C", o.get("C_at"), "| WE/WEsat", o["WE_over_WEsat_at"],
              "| g_early", o.get("gamma_perp_early"), "g_post", o.get("gamma_perp_post"),
              "g_kick", o.get("gamma_perp_kick_fit"), flush=True)
    # sicrama (bounce) frekansi tahmini: 1D L40 tohumlu kosu, t = 20-40 zaman-ortalamali phi;
    # derin tuzaklanmis elektron icin omega_b ~ k sqrt(dphi/2), dphi = max - min (sinuzoidal kuyu yaklasimi)
    bounce = None
    r = load("1d_L40_seed")
    if r is not None:
        z, _ = r
        t, ph = z["t"], tavg(z["phi"])
        sel = (t >= 18) & (t <= 26)       # 4 delik henuz birlesmeden (k = 2 pi 4 / 40)
        dphi = float(np.mean([p.max() - p.min() for p in ph[sel]]))
        k4 = 2 * np.pi * 4 / 40.0
        bounce = dict(dphi=dphi, k=k4, omega_b=float(k4 * np.sqrt(dphi / 2)), window=[18, 26],
                      note="sinuzoidal kuyu yaklasimi, kaba mertebe")
        print("bounce:", bounce, flush=True)
    with open(os.path.join(RES, "analysis.json"), "w", encoding="utf-8") as fh:
        json.dump(dict(phi_ref=refs, tavg_snapshots=TAVG, bounce=bounce, runs=res), fh, indent=1)
    np.savez_compressed(os.path.join(RES, "analysis_series.npz"), **series)


if __name__ == "__main__":
    main()
