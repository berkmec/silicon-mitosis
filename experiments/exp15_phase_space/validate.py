"""Exp15 dogrulama: (a) soguk Langmuir, (b) soguk iki-demet buyume orani, (c) Landau sonumu, (d) enerji.

Cikti: results/exp15_phase_space/validation.json, validation_series.npz
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
import time
import numpy as np
from scipy.optimize import curve_fit
from scipy.special import ndtri
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pic

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "results", "exp15_phase_space")
OUT = os.path.abspath(OUT)
os.makedirs(OUT, exist_ok=True)
series = {}


def energy_drift(WE, KE):
    Et = np.asarray(WE) + np.asarray(KE)
    return float(np.max(np.abs(Et - Et[0])) / Et[0]), float((Et[-1] - Et[0]) / Et[0])


# ---------------------------------------------------------------- (a) Langmuir
def langmuir(m, Ng=64, L=4 * np.pi, ppc=64, dt=0.1, T=80.0, amp=1e-3, tag=""):
    N = Ng * ppc
    x0 = (np.arange(N) + 0.5) * (L / N)
    k = 2 * np.pi * m / L
    x = x0 + (amp / k) * np.sin(k * x0)
    sim = pic.PIC([L], [Ng], x[None, :], np.zeros((1, N)), dt)
    sim.init_leapfrog()
    nst = int(round(T / dt))
    ts, s, WE, KE = [], [], [], []
    ref = None
    for _ in range(nst):
        o = sim.step(diag=True)
        Ek = np.fft.rfft(o["E"][0])[m]
        if ref is None:
            ref = Ek / abs(Ek)
        ts.append(o["t"]); s.append((Ek * np.conj(ref)).real)
        WE.append(o["WE"]); KE.append(o["KEc"].sum())
    ts, s = np.array(ts), np.array(s)
    # sifir gecislerinden kaba frekans, sonra en kucuk kareler
    zc = np.where(np.diff(np.sign(s)) != 0)[0]
    tz = ts[zc] - s[zc] * (ts[zc + 1] - ts[zc]) / (s[zc + 1] - s[zc])
    w_zc = np.pi / np.mean(np.diff(tz))
    f = lambda t, A, w, p: A * np.cos(w * t + p)
    popt, pcov = curve_fit(f, ts, s, p0=[np.max(np.abs(s)), w_zc, 0.0])
    w_fit, w_err = abs(popt[1]), float(np.sqrt(pcov[1, 1]))
    wg, w_pred = pic.langmuir_freq_pred(k, L / Ng, dt)
    dmax, dend = energy_drift(WE, KE)
    series[f"langmuir_m{m}"] = np.vstack([ts, s])
    return dict(m=m, k=k, k_dx=k * L / Ng, Ng=Ng, ppc=ppc, dt=dt, T=T, amp=amp,
                omega_zero_cross=float(w_zc), omega_fit=float(w_fit), omega_fit_err=w_err,
                omega_pred_grid=wg, omega_pred_grid_leapfrog=w_pred,
                rel_err_vs_wp=float(w_fit - 1.0), rel_err_vs_pred=float(w_fit / w_pred - 1.0),
                energy_drift_max=dmax, energy_drift_end=dend)


# ---------------------------------------------------------------- (b) soguk iki demet
def two_stream_cold(nwave=4, Ng=256, ppcb=64, dt=0.1, T=40.0, eps=1e-6, mmax=8, seed=3):
    kmax = np.sqrt(3.0 / 8.0)          # v0 = 1
    L = nwave * 2 * np.pi / kmax
    Nb = Ng * ppcb
    xb = (np.arange(Nb) + 0.5) * (L / Nb)
    x = np.concatenate([xb, xb + 0.5 * L / Nb])
    v = np.concatenate([np.ones(Nb), -np.ones(Nb)])
    rng = np.random.default_rng(seed)
    ph = rng.uniform(0, 2 * np.pi, mmax)
    for m in range(1, mmax + 1):
        k = 2 * np.pi * m / L
        x = x + (eps / k) * np.sin(k * x + ph[m - 1])
    sim = pic.PIC([L], [Ng], x[None, :], v[None, :], dt)
    sim.init_leapfrog()
    nst = int(round(T / dt))
    ts, amps, WE, KE = [], [], [], []
    for _ in range(nst):
        o = sim.step(diag=True)
        Ek = np.fft.rfft(o["E"][0])[1:mmax + 1]
        ts.append(o["t"]); amps.append(np.abs(Ek)); WE.append(o["WE"]); KE.append(o["KEc"].sum())
    ts, amps = np.array(ts), np.array(amps)
    # dogrusal pencere: t >= 10 ve baskin mod genligi doygunluktan ~2 basamak asagida kalana dek
    WEa = np.array(WE)
    t_end_lin = ts[np.argmax(WEa > 1e-3 * WEa.max())]
    res = []
    for m in range(1, mmax + 1):
        k = 2 * np.pi * m / L
        g_th = pic.cold_two_stream_growth(k)
        sel = (ts >= 10.0) & (ts <= t_end_lin)
        if g_th > 0 and sel.sum() > 10:
            p = np.polyfit(ts[sel], np.log(amps[sel, m - 1]), 1)
            g_meas = float(p[0])
        else:
            g_meas = float("nan")
        res.append(dict(m=m, k=float(k), kv0=float(k), gamma_theory=g_th, gamma_meas=g_meas,
                        rel_err=(g_meas / g_th - 1.0) if g_th > 0 else None))
    dmax, dend = energy_drift(WE, KE)
    series["twostream_t"] = ts
    series["twostream_amps"] = amps
    best = max([r for r in res if r["gamma_theory"] > 0], key=lambda r: r["gamma_theory"])
    # tek-mod tohumlu kosular: her mod ayri, daha uzun dogrusal pencere
    single = []
    for m in range(1, mmax + 1):
        k = 2 * np.pi * m / L
        g_th = pic.cold_two_stream_growth(k)
        if g_th <= 0:
            continue
        x = np.concatenate([xb, xb + 0.5 * L / Nb])
        x = x + (eps / k) * np.sin(k * x)
        sim = pic.PIC([L], [Ng], x[None, :], v[None, :].copy(), dt)
        sim.init_leapfrog()
        tt, aa, W1, K1 = [], [], [], []
        Tm = min(120.0, 10.0 + np.log(3e-3 / (eps / k)) / g_th + 10.0)
        for _ in range(int(round(Tm / dt))):
            o = sim.step(diag=True)
            Ek = np.fft.rfft(o["E"][0])
            tt.append(o["t"]); aa.append(2 * abs(Ek[m]) / Ng); W1.append(o["WE"]); K1.append(o["KEc"].sum())
            if 2 * np.max(np.abs(Ek[1:])) / Ng > 3e-3:
                break
        tt, aa = np.array(tt), np.array(aa)
        # baslangic: sonen ve salinan koklerin katkisi e^-4'un altina insin diye t1 = max(10, 4/gamma)
        t1 = max(10.0, 4.0 / g_th)
        sel = tt >= t1
        g_meas = float(np.polyfit(tt[sel], np.log(aa[sel]), 1)[0])
        g_meas10 = float(np.polyfit(tt[tt >= 10.0], np.log(aa[tt >= 10.0]), 1)[0])
        d1, _ = energy_drift(W1, K1)
        single.append(dict(m=m, k=float(k), gamma_theory=g_th, gamma_meas=g_meas, rel_err=g_meas / g_th - 1.0,
                           gamma_meas_from_t10=g_meas10,
                           fit_window=[t1, float(tt[-1])], energy_drift_max=d1))
    return dict(L=L, Ng=Ng, ppc_per_beam=ppcb, dt=dt, T=T, eps=eps, fit_window=[10.0, float(t_end_lin)],
                modes=res, most_unstable=best, single_mode=single,
                energy_drift_max=dmax, energy_drift_end=dend)


# ---------------------------------------------------------------- (c) Landau
def landau(nbits=20, Ng=64, k=0.5, alpha=0.01, dt=0.1, T=50.0):
    N = 1 << nbits
    L = 2 * np.pi / k
    x0 = (np.arange(N) + 0.5) * (L / N)
    x = x0 + (alpha / k) * np.sin(k * x0)
    v = ndtri(pic.bitrev(nbits, N))           # sessiz baslangic: bit-ters Maxwell (v_th = 1)
    sim = pic.PIC([L], [Ng], x[None, :], v[None, :], dt)
    sim.init_leapfrog()
    nst = int(round(T / dt))
    ts, s, a1, WE, KE = [], [], [], [], []
    ref = None
    for _ in range(nst):
        o = sim.step(diag=True)
        Ek = np.fft.rfft(o["E"][0])[1]
        if ref is None:
            ref = Ek / abs(Ek)
        ts.append(o["t"]); s.append((Ek * np.conj(ref)).real); a1.append(abs(Ek))
        WE.append(o["WE"]); KE.append(o["KEc"].sum())
    ts, s, a1 = np.array(ts), np.array(s), np.array(a1)
    root = pic.landau_root(k)
    # 1) izdusumlu gercek sinyalin |s| tepeleri: gurultu tabani t_floor'u bul (tepeler azalmayi birakana dek)
    sa = np.abs(s)
    pk = np.where((sa[1:-1] > sa[:-2]) & (sa[1:-1] >= sa[2:]))[0] + 1
    pk = pk[(ts[pk] > 6.0)]
    keep = [pk[0]]
    for p in pk[1:]:
        if sa[p] < sa[keep[-1]]:
            keep.append(p)
        else:
            break
    t_floor = float(ts[keep[-1]])
    # tepe yontemi (capraz kontrol): t in [6, 25] (gurultu tabanindan uzak), parabolik alt-adim inceltme
    pk25 = pk[ts[pk] <= 25.0]
    tp, ap = [], []
    for i in pk25:
        a_, b_, c_ = sa[i - 1], sa[i], sa[i + 1]
        den = a_ - 2 * b_ + c_
        off = 0.5 * (a_ - c_) / den if den != 0 else 0.0
        tp.append(ts[i] + off * dt)
        ap.append(b_ - 0.25 * (a_ - c_) * off)
    tp, ap = np.array(tp), np.array(ap)
    gp = np.polyfit(tp, np.log(ap), 1)[0]
    w_peaks = np.pi / np.mean(np.diff(tp))
    keep = pk25
    # 2) A e^{g t} cos(w t + p) en kucuk kareler, [6, t_floor]
    sel = (ts >= 6.0) & (ts <= t_floor)
    f = lambda t, A, g, w, p: A * np.exp(g * t) * np.cos(w * t + p)
    p0 = [np.max(np.abs(s[sel])) * np.exp(0.15 * 6), root.imag, root.real, 0.0]
    best = None
    for p_init in np.linspace(0, 2 * np.pi, 8, endpoint=False):
        try:
            popt, pcov = curve_fit(f, ts[sel], s[sel], p0=[p0[0], p0[1], p0[2], p_init], maxfev=20000)
            r = np.sum((f(ts[sel], *popt) - s[sel]) ** 2)
            if best is None or r < best[0]:
                best = (r, popt, pcov)
        except RuntimeError:
            pass
    _, popt, pcov = best
    dmax, dend = energy_drift(WE, KE)
    series["landau"] = np.vstack([ts, s, a1])
    return dict(N=N, Ng=Ng, k=k, alpha=alpha, dt=dt, T=T, quiet_start="bit-reversed Maxwell + duzgun konum",
                theory_omega=float(root.real), theory_gamma=float(root.imag),
                peaks_gamma=float(gp), peaks_omega=float(w_peaks), n_peaks=int(len(keep)), t_noise_floor=t_floor,
                fit_gamma=float(popt[1]), fit_gamma_err=float(np.sqrt(pcov[1, 1])),
                fit_omega=float(abs(popt[2])), fit_omega_err=float(np.sqrt(pcov[2, 2])),
                fit_window=[6.0, t_floor], energy_drift_max=dmax, energy_drift_end=dend)


if __name__ == "__main__":
    t0 = time.time()
    res = {}
    res["langmuir"] = [langmuir(1), langmuir(4), langmuir(1, dt=0.05)]
    print("Langmuir:", [(r["m"], r["dt"], round(r["omega_fit"], 6), round(r["omega_pred_grid_leapfrog"], 6)) for r in res["langmuir"]], flush=True)
    res["two_stream_cold"] = two_stream_cold()
    print("Two-stream:", [(r["m"], round(r["gamma_theory"], 4), round(r["gamma_meas"], 4)) for r in res["two_stream_cold"]["modes"]], flush=True)
    print("Two-stream single:", [(r["m"], round(r["gamma_theory"], 4), round(r["gamma_meas"], 4), r["fit_window"]) for r in res["two_stream_cold"]["single_mode"]], flush=True)
    res["landau"] = landau()
    lr = res["landau"]
    print("Landau: th", lr["theory_omega"], lr["theory_gamma"], "fit", lr["fit_omega"], lr["fit_gamma"],
          "peaks", lr["peaks_omega"], lr["peaks_gamma"], "t_floor", lr["t_noise_floor"], flush=True)
    # warm two-stream dispersion (deney parametreleri, v0=1, vth=0.2) referans tablosu
    ks = np.linspace(0.05, 1.1, 211)
    gw = np.array([pic.warm_two_stream_growth(k, 1.0, 0.2) for k in ks])
    gc = np.array([pic.cold_two_stream_growth(k) for k in ks])
    series["disp_k"] = ks; series["disp_warm"] = gw; series["disp_cold"] = gc
    i = int(np.nanargmax(gw))
    res["warm_two_stream_ref"] = dict(v0=1.0, vth=0.2, gamma_max=float(gw[i]), k_at_max=float(ks[i]),
                                      lambda_at_max=float(2 * np.pi / ks[i]))
    res["runtime_s"] = time.time() - t0
    with open(os.path.join(OUT, "validation.json"), "w", encoding="utf-8") as fh:
        json.dump(res, fh, indent=1, ensure_ascii=False)
    np.savez_compressed(os.path.join(OUT, "validation_series.npz"), **series)
    print("done", res["runtime_s"])
