"""Exp11 physical mapping (order-of-magnitude; every input below is an ASSUMPTION or a cited/remembered value).

A (supercapacitor idea): summarise model_a_phys.json (threshold gate capacitance C*(n) for the uniform 2D electron
   layer to become unstable within q <= k_F).
B (dissipative cells): cell radius R ~ rho * l, l = sqrt(D tau) (rho = R_div / l from the 2D simulation).
   Minimum power per chip area to keep a fraction f of a layer covered with driven cells of radius R:
       P/A ~ f * (4/3) * n * E_event * D * rho^2 / R            (every carrier in a cell turns over once per tau = l^2/D)
   Continuum (pattern) validity: N = n (4/3) pi R^3 >= N_min carriers per cell.
   Quantum target at 4 K (isolated sphere in Si): E_C = e^2/(8 pi eps0 eps R) >= 10 kT and
   Delta E = hbar^2 pi^2 / (2 m R^2) >= 10 kT."""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "exp11"
E, KB, EPS0, HBAR, M0 = 1.602176634e-19, 1.380649e-23, 8.8541878128e-12, 1.054571817e-34, 9.1093837015e-31

l_sim = 10.0  # l = sqrt(D/k) = 10 w in the simulations (k = 0.01)
rho, rho_src = 1.0, "varsayılan (simülasyondan okunamadı)"
p = OUT / "b3d.json"
if p.exists():  # 3D binary division: radius of the single droplet in the record just before N = 2
    h3 = json.loads(p.read_text(encoding="utf-8"))["hist"]
    i = next((j for j, x in enumerate(h3) if x["n"] >= 2), None)
    if i:
        rho, rho_src = max(h3[i - 1]["radii"]) / l_sim, "3B Zwicker ikili bölünme (b3d.py), bölünmeden önceki kayıt"
rho2d = None
p = OUT / "b_rho.json"
if p.exists():
    rho2d = json.loads(p.read_text(encoding="utf-8"))["rho"]  # 2D worm break (larger -> needs ~(rho2d/rho)^2 more power)

# ---------------------------------------------------------------- quantum target sizes
eps_si, m_si = 11.7, 0.26
res = {"assumptions": {}, "quantum_target": {}, "B": {}, "A": {}}
for T in (4.0, 300.0):
    R_ec = E**2 / (8 * np.pi * EPS0 * eps_si * 10 * KB * T)
    R_de = np.sqrt(HBAR**2 * np.pi**2 / (2 * m_si * M0 * 10 * KB * T))
    res["quantum_target"][f"{T:g}K"] = dict(R_EC_nm=R_ec * 1e9, R_dE_nm=R_de * 1e9)
    print(f"{T:5.0f} K kuantum hedefi: E_C>=10kT için R <= {R_ec*1e9:.3g} nm; ΔE>=10kT için R <= {R_de*1e9:.3g} nm")

# ---------------------------------------------------------------- B scenarios
f_cov, N_min = 0.3, 100
scen = {
    "e-h sıvısı (n=3,3e18, E=E_g=1,12 eV)": (3.3e18, 1.12),
    "donmuş katkı darbe iyonlaşması (n=1,7e17, E=45 meV)": (1.7e17, 0.045),
}
res["assumptions"] = dict(f_coverage=f_cov, N_min=N_min, rho=rho, rho_source=rho_src, rho_2D_worm=rho2d,
                          power_factor_2D_worm=(rho2d / rho) ** 2 if rho2d else None,
                          D_cm2_s=[0.1, 1.0, 10.0], P_max_W_cm2=[0.1, 1.0, 10.0],
                          note="n, E: e-h sıvısı yoğunluğu literatür (3,33e18), E_g ve 45 meV (B alıcı) bellekten; "
                               "D, P_max, f, N_min varsayım. Joule ısınması dahil değil (alt sınır).")
for name, (n_cm3, E_eV) in scen.items():
    n = n_cm3 * 1e6
    Ev = E_eV * E
    R_N = (3 * N_min / (4 * np.pi * n)) ** (1 / 3)
    rows = []
    for D_cm2 in (0.1, 1.0, 10.0):
        D = D_cm2 * 1e-4
        for Pmax in (0.1, 1.0, 10.0):
            R_P = f_cov * (4 / 3) * n * Ev * D * rho**2 / (Pmax * 1e4)
            P20 = f_cov * (4 / 3) * n * Ev * D * rho**2 / 20e-9 / 1e4  # W/cm^2 needed at R = 20 nm
            tau20 = (20e-9 / rho) ** 2 / D
            rows.append(dict(D_cm2_s=D_cm2, Pmax_W_cm2=Pmax, R_N_nm=R_N * 1e9, R_P_nm=R_P * 1e9,
                             R_min_nm=max(R_N, R_P) * 1e9, P_at_20nm_W_cm2=P20, tau_at_20nm_s=tau20))
    res["B"][name] = rows
    best = min(rows, key=lambda r: r["R_min_nm"])
    print(f"{name}: N>={N_min} için R >= {R_N*1e9:.3g} nm; R=20 nm'de gereken güç "
          f"{min(r['P_at_20nm_W_cm2'] for r in rows):.3g}-{max(r['P_at_20nm_W_cm2'] for r in rows):.3g} W/cm² "
          f"(τ = {rows[3]['tau_at_20nm_s']:.2g} s @ D=1); en iyimser R_min = {best['R_min_nm']:.3g} nm "
          f"(D={best['D_cm2_s']}, P_max={best['Pmax_W_cm2']} W/cm²)")

# ---------------------------------------------------------------- review addition: cells in a 2D layer (MOSFET channel)
# disk of radius R in a layer of areal density n_s: N = n_s pi R^2, P/A = f n_s E D rho^2 / R^2.
# The smallest R allowed by both limits is reached at n_s* = sqrt(N_min P / (pi f E D rho^2)):
#   R* = (N_min f E D rho^2 / (pi P))^(1/4)
C_ox = EPS0 * 3.9 / 4.148e-9  # SKY130 gate oxide, F/m^2
R_gate = np.sqrt(E**2 / (2 * C_ox * 10 * KB * 4.0) / np.pi)  # E_C >= 10 kT at 4 K, parallel-plate gate only
res["quantum_target"]["4K_gated_disk_parallel_plate_nm"] = R_gate * 1e9
print(f"  4 K kuantum hedefi, SKY130 kapılı 2B disk (yalnız paralel plaka, kenar alanı yok): R <= {R_gate*1e9:.3g} nm")
res["B_2D_layer"] = {}
for name, E_eV in (("2B tabaka, bant-bant darbe iyonlaşması (E=E_g=1,12 eV)", 1.12),
                   ("2B tabaka, sığ katkı darbe iyonlaşması (E=45 meV)", 0.045)):
    rows = []
    for D_cm2 in (0.1, 1.0, 10.0):
        for Pmax in (0.1, 1.0, 10.0):
            D, P, Ev = D_cm2 * 1e-4, Pmax * 1e4, E_eV * E
            ns_opt = np.sqrt(N_min * P / (np.pi * f_cov * Ev * D * rho**2))
            R_opt = (N_min * f_cov * Ev * D * rho**2 / (np.pi * P)) ** 0.25
            rows.append(dict(D_cm2_s=D_cm2, Pmax_W_cm2=Pmax, ns_opt_cm2=ns_opt * 1e-4, R_min_nm=R_opt * 1e9))
    res["B_2D_layer"][name] = rows
    best = min(rows, key=lambda r: r["R_min_nm"])
    typ = next(r for r in rows if r["D_cm2_s"] == 1.0 and r["Pmax_W_cm2"] == 1.0)
    print(f"{name}: en iyimser R* = {best['R_min_nm']:.3g} nm (n_s* = {best['ns_opt_cm2']:.2g} cm⁻², D={best['D_cm2_s']}, "
          f"P={best['Pmax_W_cm2']} W/cm²); tipik (D=1, 1 W/cm²) R* = {typ['R_min_nm']:.3g} nm")
    sens = {str(Nm): best["R_min_nm"] * (Nm / N_min) ** 0.25 for Nm in (10, 30, 100, 300)}
    res["B_2D_layer"][name + " | N_min duyarlılığı (en iyimser köşe), nm"] = sens
    print("   N_min duyarlılığı (en iyimser köşe, R* ∝ N_min^1/4):", {k: round(v, 1) for k, v in sens.items()})

# ---------------------------------------------------------------- A summary from model_a_phys.json
a = json.loads((OUT / "model_a_phys.json").read_text(encoding="utf-8"))
for key, cur in a["curves"].items():
    cs = [(n, c) for n, c in zip(cur["n_cm2"], cur["C_star_uF_cm2"]) if c is not None]
    res["A"][key] = dict(
        unstable_any=bool(cs),
        min_C_star=min(c for _, c in cs) if cs else None,
        n_max_at_SKY130=max([n for n, c in cs if c <= a["C_SKY130_uF_cm2"]], default=None),
        n_max_at_10uF=max([n for n, c in cs if c <= 10.0], default=None))
    print(f"A {key:26s}: {res['A'][key]}")
(OUT / "phys.json").write_text(json.dumps(res, indent=1, ensure_ascii=False), encoding="utf-8")
print("rho =", rho, "(", rho_src, ");  2B solucan rho =", rho2d)
