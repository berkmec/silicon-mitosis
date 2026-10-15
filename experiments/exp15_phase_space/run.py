"""Exp15 deney kosulari: iki-demet kararsizligi -> elektron delikleri; 1D / 2D / 3D karsilastirma.

Kullanim:  python run.py ad1 ad2 ...   (sirayla, tek surec; dusuk oncelik)
Cikti:     results/exp15_phase_space/runs/<ad>.npz ve <ad>.json
Fizik:     simetrik Maxwell demetleri +-v0 (x yonunde), her biri n = 1/2, izotrop v_th; v0 = 1, v_th = 0.2.
           Uzunluk birimi v0/omega_p, zaman 1/omega_p. B (varsa) +x yonunde, Omega_ce = B.
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
from scipy.special import ndtri
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pic

OUT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..",
                                   "results", "exp15_phase_space", "runs"))
os.makedirs(OUT, exist_ok=True)

V0, VTH, DT, SNAP = 1.0, 0.2, 0.1, 1.0
PS_TIMES = [0, 10, 15, 20, 30, 50, 75, 100, 150, 200, 300, 500, 750, 1000]
CONFIGS = {}


def add(name, **kw):
    base = dict(B=0.0, mseed=None, amp=0.0, smooth=0)
    base.update(kw)
    CONFIGS[name] = base


# NOT: CPU baska oturumlarin isleriyle paylasildigi icin (4 cekirdek dolu) parcacik sayilari/sureler
# ilk plana gore dusuruldu (2D genis ppc 64 -> 32, 3D ppc 32 -> 16, 1D L100 ppc 400 -> 200, t 1000 -> 600).
# 1D: gurultuden baslayan, ~10 dalga boylu kutu (birlesme sayimi), 3 tohum
for s in (1, 2, 3):
    add(f"1d_L100_s{s}", L=(100.0,), N=(512,), ppc=200, nv=1, seed=s, t_end=600.0)
# L = 40 ailesi (4 dalga boyu, m = 4 duzlemsel tohum): 1D, 2D dar (quasi-1D), 2D genis, 2D genis + B, tohumsuz
add("1d_L40_seed", L=(40.0,), N=(128,), ppc=400, nv=1, seed=11, mseed=4, amp=0.05, t_end=300.0)
add("2d_L40_narrow", L=(40.0, 1.25), N=(128, 4), ppc=32, nv=2, seed=11, mseed=4, amp=0.05, t_end=300.0)
add("2d_L40_narrow_ppc64", L=(40.0, 1.25), N=(128, 4), ppc=64, nv=2, seed=11, mseed=4, amp=0.05, t_end=300.0)
add("2d_L40_wide", L=(40.0, 40.0), N=(128, 128), ppc=32, nv=2, seed=11, mseed=4, amp=0.05, t_end=150.0)
add("2d_L40_wide_B2", L=(40.0, 40.0), N=(128, 128), ppc=32, nv=3, seed=11, mseed=4, amp=0.05, t_end=150.0, B=2.0)
add("2d_L40_wide_noseed", L=(40.0, 40.0), N=(128, 128), ppc=32, nv=2, seed=12, t_end=150.0)
# temiz enine test: y'de tam simetrik sessiz baslangic (delikler kesin duzlemsel olusur), doygunluktan sonra
# t = 25'te kucuk rastgele y-itmesi (tum k_y modlarini tohumlar) -> olusmus deliklerin enine kararsizligi
add("2d_L40_wide_quietkick", L=(40.0, 40.0), N=(128, 128), ppc=32, nv=2, seed=13, mseed=4, amp=0.05, t_end=150.0,
    quiet_y=True, kick_t=25.0, kick_amp=0.05)
add("2d_L40_wide_quietkick_B2", L=(40.0, 40.0), N=(128, 128), ppc=32, nv=3, seed=13, mseed=4, amp=0.05,
    t_end=150.0, quiet_y=True, kick_t=25.0, kick_amp=0.05, B=2.0)
# L = 20 ailesi (2 dalga boyu, m = 2 tohum), 3D ile ayni cozunurluk (dx = 0.625)
add("1d_L20_seed", L=(20.0,), N=(32,), ppc=400, nv=1, seed=21, mseed=2, amp=0.05, t_end=150.0)
add("2d_L20_seed", L=(20.0, 20.0), N=(32, 32), ppc=16, nv=2, seed=21, mseed=2, amp=0.05, t_end=150.0)
add("2d_L20_seed_ppc32", L=(20.0, 20.0), N=(32, 32), ppc=32, nv=2, seed=21, mseed=2, amp=0.05, t_end=150.0)
add("3d_L20_seed", L=(20.0, 20.0, 20.0), N=(32, 32, 32), ppc=16, nv=3, seed=21, mseed=2, amp=0.05, t_end=120.0)
add("3d_L20_seed_B2", L=(20.0, 20.0, 20.0), N=(32, 32, 32), ppc=16, nv=3, seed=21, mseed=2, amp=0.05,
    t_end=120.0, B=2.0)
# temiz 2D-3D karsilastirmasi: dik yonlerde tam simetrik baslangic + t = 25'te rastgele dik itme
add("2d_L20_quietkick", L=(20.0, 20.0), N=(32, 32), ppc=16, nv=2, seed=23, mseed=2, amp=0.05, t_end=120.0,
    quiet_y=True, kick_t=25.0, kick_amp=0.05)
add("3d_L20_quietkick", L=(20.0, 20.0, 20.0), N=(32, 32, 32), ppc=16, nv=3, seed=23, mseed=2, amp=0.05,
    t_end=120.0, quiet_y=True, kick_t=25.0, kick_amp=0.05)


def run(name):
    c = CONFIGS[name]
    L, N = np.array(c["L"], float), tuple(c["N"])
    d = len(N)
    npart = int(c["ppc"] * np.prod(N))
    rng = np.random.default_rng(c["seed"])
    if c.get("quiet_y"):
        # taban: x kafesi + bit-ters Maxwell hizlari (sessiz 1D), dik yonlerde hucre basina 1 esit aralikli kopya
        # (M = N[1]*...*N[d-1]). CIC birim bolusumu nedeniyle yogunluk dik koordinatlardan tam bagimsiz;
        # dinamik dik oteleme simetrisini korur.
        M = int(np.prod(N[1:]))
        nbb = npart // M // 2
        nbits = int(round(np.log2(nbb)))
        assert 2 ** nbits == nbb
        q = pic.bitrev(nbits, nbb)
        vxt = VTH * ndtri(q)
        vyt = VTH * ndtri(q[rng.permutation(nbb)])
        xb0 = (np.arange(nbb) + 0.5) * L[0] / nbb
        xb = np.concatenate([xb0, xb0 + 0.5 * L[0] / nbb])[None, :]
        vb = np.zeros((c["nv"], 2 * nbb))
        vb[0] = np.concatenate([V0 + vxt, -V0 + vxt[::-1]])
        vb[1] = np.concatenate([vyt, vyt[::-1]])
        if c["nv"] > 2:
            vb[2] = np.concatenate([vyt[rng.permutation(nbb)], vyt[rng.permutation(nbb)]])
        if c["mseed"]:
            xb = pic.displace_x(xb, L[0], [c["mseed"]], [c["amp"]])
        x = np.empty((d, 2 * nbb * M))
        x[0] = np.tile(xb[0], M)
        grids = np.meshgrid(*[(np.arange(N[a]) + 0.5) * L[a] / N[a] for a in range(1, d)], indexing="ij")
        for a in range(1, d):
            x[a] = np.repeat(grids[a - 1].ravel(), 2 * nbb)
        v = np.tile(vb, (1, M))
        npart = x.shape[1]
    else:
        x, v = pic.beams(L, npart, V0, VTH, c["nv"], rng)
        if c["mseed"]:
            x = pic.displace_x(x, L[0], [c["mseed"]], [c["amp"]])
    kicked = False
    sim = pic.PIC(L, N, x, v, DT, B=c["B"], smooth=c["smooth"])
    del x, v
    sim.init_leapfrog()
    nsteps = int(round(c["t_end"] / DT))
    every = int(round(SNAP / DT))
    nxb = min(N[0], 256)
    vbins = np.linspace(-2.5, 2.5, 101)
    fbins = np.linspace(-3.0, 3.0, 121)
    ts, phis, WE, KE = [], [], [], []
    ps_t, ps_h, fvx, fvp = [], [], [], []
    t0 = time.time()
    for n in range(nsteps + 1):
        if c.get("kick_t") is not None and not kicked and sim.t >= c["kick_t"] - 0.5 * DT:
            for a in range(1, d):
                sim.x[a] += c["kick_amp"] * rng.standard_normal(sim.np)
            sim._wrap(full=True)
            kicked = True
            print(f"[{name}] t={sim.t:.2f}: y-itmesi uygulandi (std {c['kick_amp']})", flush=True)
        diag = (n % every == 0)
        o = sim.step(diag=diag)
        if not diag:
            continue
        t = o["t"]
        ts.append(t)
        phis.append(o["phi"].astype(np.float32))
        WE.append(o["WE"])
        KE.append(o["KEc"])
        if any(abs(t - tp) < 0.5 * DT for tp in PS_TIMES):
            # faz uzayi x-vx (dik koordinatlar uzerinden toplanmis); x leapfrog sonrasi t+dt'de (ihmal edilebilir)
            H, _, _ = np.histogram2d(sim.x[0], sim.v[0], bins=[nxb, vbins], range=[[0, L[0]], [-2.5, 2.5]])
            ps_t.append(t); ps_h.append(H.astype(np.float32))
            fvx.append(np.histogram(sim.v[0], bins=fbins)[0])
            if c["nv"] > 1:
                fvp.append(np.histogram(sim.v[1], bins=fbins)[0])
        if len(ts) % 10 == 1:
            Et = WE[-1] + KE[-1].sum()
            E0 = WE[0] + KE[0].sum()
            print(f"[{name}] t={t:7.1f} WE={WE[-1]:.4e} dE/E={(Et - E0) / E0:+.2e} "
                  f"Kperp={KE[-1][1:].sum() if c['nv'] > 1 else 0:.4e} wall={time.time() - t0:.0f}s", flush=True)
    runtime = time.time() - t0
    ts = np.array(ts); WE = np.array(WE); KE = np.array(KE)
    Et = WE + KE.sum(axis=1)
    drift_max = float(np.max(np.abs(Et - Et[0])) / Et[0])
    drift_end = float((Et[-1] - Et[0]) / Et[0])
    np.savez_compressed(os.path.join(OUT, name + ".npz"), t=ts, phi=np.array(phis), WE=WE, KE=KE,
                        ps_t=np.array(ps_t), ps_h=np.array(ps_h), fvx=np.array(fvx),
                        fvp=np.array(fvp) if fvp else np.zeros(0), vbins=vbins, fbins=fbins,
                        L=L, N=np.array(N))
    meta = dict(name=name, config=c, v0=V0, vth=VTH, dt=DT, snap=SNAP, npart=npart, d=d,
                energy_drift_max=drift_max, energy_drift_end=drift_end, runtime_s=runtime,
                WE_max=float(WE.max()), t_WE_max=float(ts[np.argmax(WE)]))
    with open(os.path.join(OUT, name + ".json"), "w", encoding="utf-8") as fh:
        json.dump(meta, fh, indent=1)
    print(f"[{name}] bitti: {runtime:.0f}s, enerji kaymasi max {drift_max:.2e}", flush=True)
    del sim


if __name__ == "__main__":
    for nm in sys.argv[1:]:
        if os.path.exists(os.path.join(OUT, nm + ".json")) and "--force" not in sys.argv:
            print(f"[{nm}] zaten var, atlandi", flush=True)
            continue
        if nm.startswith("--"):
            continue
        run(nm)
    print("HEPSI BITTI", flush=True)
