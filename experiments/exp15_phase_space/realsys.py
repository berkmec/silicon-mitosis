"""Exp15: normalize PIC birimlerinin gercek sistemlere eslenmesi (mertebe hesabi).

Cikti: results/exp15_phase_space/realsys.json
Ornek parametreler ACIKCA varsayimdir (ozel bir cihaz/olcum iddiasi yok).
"""
import os
import json
import numpy as np

e, me, eps0, kB = 1.602176634e-19, 9.1093837015e-31, 8.8541878128e-12, 1.380649e-23
hbar, hP = 1.054571817e-34, 6.62607015e-34
OUT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..",
                                   "results", "exp15_phase_space", "realsys.json"))


def plasma(n_m3, T_eV, B_T, lnL=10.0, m=me, eps_r=1.0):
    wp = np.sqrt(n_m3 * e ** 2 / (eps0 * eps_r * m))
    lD = np.sqrt(eps0 * eps_r * T_eV * e / (n_m3 * e ** 2))
    vth = np.sqrt(T_eV * e / m)
    wc = e * B_T / m
    # NRL: elektron carpisma frekansi nu_e = 2.91e-6 n[cm^-3] lnL T_eV^-1.5 [1/s] (yalniz eps_r = 1 icin)
    nu = 2.91e-6 * (n_m3 * 1e-6) * lnL * T_eV ** -1.5
    return dict(n_m3=n_m3, T_eV=T_eV, B_T=B_T, omega_p=wp, f_p_Hz=wp / (2 * np.pi), lambda_D_m=lD, v_th_m_s=vth,
                Omega_ce_over_omega_p=wc / wp, nu_e_over_omega_p=nu / wp, mfp_over_lambda_D=(vth / nu) / lD,
                N_D=n_m3 * lD ** 3, hole_size_m_for_10_lambdaD=10 * lD)


def silicon(n_cm3=1e17, mstar=0.26, eps_r=11.7, T=4.0, mu_cm2Vs=(100.0, 1000.0)):
    n = n_cm3 * 1e6
    m = mstar * me
    eps = eps0 * eps_r
    wp = np.sqrt(n * e ** 2 / (eps * m))
    lD = np.sqrt(eps * kB * T / (n * e ** 2))
    out = dict(n_cm3=n_cm3, mstar=mstar, eps_r=eps_r, T_K=T, omega_p=wp, f_p_THz=wp / (2 * np.pi) / 1e12,
               hbar_omega_p_meV=hbar * wp / e * 1e3, kT_meV=kB * T / e * 1e3, lambda_D_classical_nm=lD * 1e9)
    for gv, md, tag in ((1, mstar, "tek_vadi_m*=0.26"), (6, 0.328, "6_vadi_md=0.328")):
        EF = hbar ** 2 * (3 * np.pi ** 2 * n / gv) ** (2.0 / 3.0) / (2 * md * me)
        lTF = np.sqrt(2 * eps * EF / (3 * n * e ** 2))
        out[f"E_F_meV_{tag}"] = EF / e * 1e3
        out[f"T_F_K_{tag}"] = EF / kB
        out[f"lambda_TF_nm_{tag}"] = lTF * 1e9
    ldb = hP / np.sqrt(2 * np.pi * m * kB * T)
    out["lambda_dB_thermal_nm"] = ldb * 1e9
    out["n_lambda_dB3"] = n * ldb ** 3
    out["N_D_classical"] = n * lD ** 3
    a = (3 / (4 * np.pi * n)) ** (1 / 3)
    aB = eps_r / mstar * 0.0529177e-9
    out["r_s"] = a / aB
    out["Gamma_classical"] = e ** 2 / (4 * np.pi * eps * a) / (kB * T)
    out["Mott_n_c_SiP_cm3"] = 3.7e18
    vth = np.sqrt(kB * T / m)
    EF1 = out["E_F_meV_tek_vadi_m*=0.26"] * 1e-3 * e
    vF = np.sqrt(2 * EF1 / m)
    out["v_thermal_m_s"] = vth
    out["v_F_m_s_tek_vadi"] = vF
    for mu in mu_cm2Vs:
        tau = mu * 1e-4 * m / e
        out[f"mu{int(mu)}"] = dict(tau_s=tau, omega_p_tau=wp * tau, mfp_nm_vF=vF * tau * 1e9,
                                   mfp_over_lambda_D=vF * tau / lD)
    # 50 lambda_D'lik delik (bu PIC'teki delikler ~10 v0/omega_p ~ 50 lambda_D,demet) icin gereken mobilite
    L_hole = 50 * lD
    tau_req = L_hole / vF
    out["hole_50lambdaD_nm"] = L_hole * 1e9
    out["mu_required_cm2Vs_for_mfp_eq_hole"] = e * tau_req / m * 1e4
    return out


if __name__ == "__main__":
    res = dict(
        notes="Ornek parametreler varsayimdir; mertebe hesabidir.",
        lab_low_T=plasma(1e16, 3.0, 0.05),
        solar_wind_1AU=plasma(5e6, 10.0, 5e-9, lnL=25.0),
        auroral_region=plasma(1e6, 100.0, 5e-6, lnL=25.0),
        silicon_4K=silicon(),
    )
    with open(OUT, "w", encoding="utf-8") as fh:
        json.dump(res, fh, indent=1)
    for k, v in res.items():
        print(k, v)
