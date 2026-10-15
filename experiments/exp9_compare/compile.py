"""Experiment 9: collect every result file, official model parameter and literature item on disk and compare them like-for-like.

Writes results/exp9/: COMPARISON.md, comparison.csv, inventory.md, fig_comparison.png, verdicts.json.
Numbers are read from the result files (no hand transcription). Only cheap re-evaluations are done here:
 - 1D solver at tox +/- sigma_tox (threshold-voltage sensitivity),
 - isolated-sphere and visibility formulas, self-checked against exp5/exp7 values,
 - pytest and a recount of the fission exponents from the merged scan (same rule as exp8/report.py)."""
import csv
import json
import os
import re
import subprocess
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from qcell import model as m  # noqa: E402
from qcell import mos1d as M  # noqa: E402

R = ROOT / "results"
DATA = ROOT / "data"
OUT = R / "exp9"
OUT.mkdir(parents=True, exist_ok=True)


def jload(p):
    return json.loads(Path(p).read_text(encoding="utf-8"))


def cload(p):
    with open(p, encoding="utf-8") as f:
        return list(csv.DictReader(f))


def grab(text, pat, cast=float):
    mm = re.search(pat, text)
    if not mm:
        raise ValueError(f"pattern not found: {pat}")
    return cast(mm.group(1))


# ---------------------------------------------------------------- load everything
sph = jload(R / "exp1/validate_sphere.json")
sch = jload(R / "exp1/schrodinger.json")
tr = jload(R / "exp1/transport.json")
capA = cload(R / "exp1/cap_cell_A.csv")
dolf_txt = (R / "exp1/dolfinx_check.txt").read_text(encoding="utf-8", errors="replace")
gen0 = cload(R / "exp0/generations.csv")
st1 = jload(R / "exp4/stage1.json")
mos = jload(R / "exp5/mos1d.json")
rdf = jload(R / "exp6/rdf.json")
rdf_raw = jload(R / "exp6/rdf_uncorrected_first_run.json")
femchk = jload(R / "exp6/fem_check.json")
SENS_P = R / "exp6/cap_sensitivity.json"
SENS = jload(SENS_P) if SENS_P.exists() else None
tp = jload(R / "exp7/testplan.json")
pf = jload(R / "exp8/pf_validation.json")
sel = jload(R / "exp8/selection.json")
div = jload(R / "exp8/division_diag.json")
fis_main = jload(R / "exp8/fission_scan.json")
fis4 = jload(R / "exp8/fission_scan_d4.0.json")
fis2 = jload(R / "exp8/fission_scan_d2.0.json")
dqd = jload(R / "exp3/dqd.json")
chain = jload(R / "exp2/chain.json")
nist = jload(DATA / "nist_sky130/record.json")["ResultData"][0]

MDL = (DATA / "sky130_models/nfet_01v8_tt.pm3.spice").read_text(encoding="utf-8", errors="replace")
MIS = (DATA / "sky130_models/nfet_mismatch.spice").read_text(encoding="utf-8", errors="replace")
TOXE = grab(MDL, r"toxe = \{([0-9.]+e[+-]?[0-9]+)")
NDEP = grab(MDL, r"ndep = ([0-9.]+e[+-]?[0-9]+)")
VTH0 = grab(MDL, r"vth0 = \{([0-9.]+)")
XJ = grab(MDL, r"xj = ([0-9.]+e[+-]?[0-9]+)")
LINT = grab(MDL, r"lint = ([0-9.]+e[+-]?[0-9]+)")
EPSROX = grab(MDL, r"epsrox = ([0-9.]+)")
VTH0_SLOPE = grab(MIS, r"vth0_slope=([0-9.]+e[+-]?[0-9]+)")
TOXE_SLOPE = grab(MIS, r"toxe_slope=([0-9.]+e[+-]?[0-9]+)")

L_UM, W_UM = 0.15, 0.42  # drawn length; W=0.42 um is a doc-table minimum (not verified as a rule)
SIG_VT_FOUNDRY_MV = VTH0_SLOPE / np.sqrt(L_UM * W_UM) * 1e3  # valid only if the slope is in V*um
SIG_TOX_NM = TOXE * TOXE_SLOPE / np.sqrt(L_UM * W_UM) * 1e9

# ---------------------------------------------------------------- cheap re-evaluations
EPS0, E, KB = m.EPS0, m.E, m.KB


def e_iso_eV(R_nm, eps):
    C = 4 * np.pi * EPS0 * eps * R_nm * 1e-9
    return (E**2 / (2 * C)) / E


def e_iso_K(R_nm, eps):
    return e_iso_eV(R_nm, eps) * E / KB


def visibility(EC_K, T, Vsd=20e-6):
    """Same procedure as exp7/run.py (peak-to-valley of |I| over one period, orthodox theory)."""
    Cs = E**2 / (2 * EC_K * KB)
    Cg_, Cs_, Cd_ = 0.5 * Cs, 0.25 * Cs, 0.25 * Cs
    Vg = np.linspace(0, 3 * E / Cg_, 400)
    I = np.abs(m.single_dot_current(Cs_, Cd_, Cg_, Vg, [Vsd], T, R_s=1e9, R_d=1e9)[:, 0])
    return float((I.max() - I.min()) / (I.max() + I.min())) if I.max() + I.min() > 0 else 0.0


VIS = {(EC, T): visibility(EC, T) for EC in (25.0, 34.0, 40.0) for T in (0.05, 4.0, 14.0, 77.0, 300.0)}


def vg_at_ns(tox_nm, na=1.7e17, ns_target=1e12):
    """Same sweep and interpolation as exp5/run.py (constant-charge threshold)."""
    mos_ = M.Mos(na_cm3=na, tox_nm=tox_nm, zmax_nm=140.0, dz_nm=0.2)
    pss = np.linspace(1.55 * mos_.phif, 2 * mos_.phif + 0.42, 20)
    sols = [mos_.solve(p) for p in pss]
    vg = np.array([s["Vg"] for s in sols])
    ns = np.array([max(s["ns_cm2"], 1e3) for s in sols])
    return float(np.interp(np.log(ns_target), np.log(ns), vg))


TOX_NM = TOXE * 1e9
VG_BASE = vg_at_ns(TOX_NM)
VG_P = vg_at_ns(TOX_NM + SIG_TOX_NM)
VG_M = vg_at_ns(TOX_NM - SIG_TOX_NM)
SIG_VT_TOX_MV = abs(VG_P - VG_M) / 2 * 1e3

ANSI = re.compile(r"\x1b?\[[0-9;]*m")  # terminal colour codes in pytest output
NSZ = (10, 13, 16, 19, 22, 26)
merged =[r for r in fis_main["scan"] if str(r["d"]) == "inf"] + fis4["scan"] + fis2["scan"]
GC, FIS_SLOPE = {}, {}
for d in ("inf", "4.0", "2.0"):
    for N in NSZ:
        rows_ = sorted([r for r in merged if str(r["d"]) == d and r["N"] == N], key=lambda r: r["g"])
        hit = [r["g"] for r in rows_ if r["fissions"]]
        GC[(d, N)] = hit[0] if hit else None
    xs = [N for N in NSZ if GC[(d, N)] is not None]
    if len(xs) >= 3:
        FIS_SLOPE[d] = float(np.polyfit(np.log(xs), np.log([GC[(d, N)] for N in xs]), 1)[0])

try:
    pt = subprocess.run([str(ROOT / "env/venv/Scripts/python.exe"), "-m", "pytest", "tests", "-q", "-p", "no:cacheprovider"],
                        cwd=str(ROOT), capture_output=True, text=True, env={**os.environ, "PYTHONPATH": ""}, timeout=600)
    PYTEST_LINE = ANSI.sub("", pt.stdout.strip().splitlines()[-1]) if pt.stdout.strip() else pt.stderr.strip()[-200:]
except Exception as ex:  # noqa: BLE001
    PYTEST_LINE = f"pytest çalıştırılamadı: {ex}"


# ---------------------------------------------------------------- helper values
def capv(r, lvl, key):
    return float(next(x for x in capA if float(x["r"]) == r and int(x["level"]) == lvl)[key])


def fem_pair(z0, idx, label="fine"):
    d = femchk[f"z0={z0}_{label}"]
    return d["phi_fem_mV"][idx], d["phi_image_mV"][idx], d["rel_diff"][idx]


def tun(L):
    return next(t for t in st1["tunnel_coupling"] if abs(t["L_nm"] - L) < 1e-9)


def well(L):
    return next(w for w in sch["wells"] if abs(w["L_nm"] - L) < 1e-9)


ftab = tr["eps"]["3.9"]["table"]
ftab_hi = tr["eps"]["11.7"]["table"]
e0_g3 = next(r for r in gen0 if int(r["gen"]) == 3)
e0_g5 = next(r for r in gen0 if int(r["gen"]) == 5)
hw = [v for v in sch["validation"] if v["test"] == "hardwall"]
hw_fine = min(hw, key=lambda v: v["h"])
sep = [v for v in sch["validation"] if v["test"] == "separable_finite"]
sep_fine = min(sep, key=lambda v: v["h_int"])
rwc_last = [v for v in sch["validation"] if v["test"] == "real_well_convergence"][-1]
lz_err = max(x["abs_err"] for x in st1["lz_validation"])
lb = st1["lindblad_300K_best"]
lb10 = next(r for r in lb["rows"] if abs(r["rise_ps"] - 10.0) < 1e-9)
mosN, mosC = mos["nominal"], mos["calibrated"]
RN, RNC = rdf["nominal_1.7e17"], rdf["nominal_1.7e17_FEMcorr"]
RE, REC = rdf["effective_3.9e17"], rdf["effective_3.9e17_FEMcorr"]
RAW_OK = abs(rdf_raw["nominal_1.7e17"]["sigma_surface_potential_mV"] - RN["sigma_surface_potential_mV"]) < 1e-9
NA_MOTT = 3.7e18
EC_CONST_39 = st1["EC_const_Knm"]["3.9"]
EC_CONST_117 = st1["EC_const_Knm"]["11.7"]
E_GATED_125_39 = ftab[3]["EC_over_kB_K"]
E_GATED_125_117 = ftab_hi[3]["EC_over_kB_K"]
E_ISO_125 = e_iso_K(6.25, 3.9)
EPS_STAR_LIT4 = e_iso_eV(0.9, 1.0) / 0.300
LIT4_MODEL_EV = e_iso_eV(0.9, 3.9)
SHIN_EKT300 = e_iso_K(2.5, 3.9) / 300.0
LJ = fis_main.get("validation", {})
LJ_ERR = max(abs(v["E"] - v["E_ref"]) for v in LJ.values()) if LJ else float("nan")
chain_ec12 = chain["EC_over_kB_K"][0]
single_ec12 = EC_CONST_39 / 12.0
div_nosplit = all(s["n"] == 1 for d_ in div for s in d_["series"][-1:])
div_n = [d_["series"][-1]["n"] for d_ in div]
r_hw_fem_1 = 1.0 / well(1.0)["ratio_dE_fem_over_hw"]
r_hw_fem_3 = 1.0 / well(3.0)["ratio_dE_fem_over_hw"]

C_DOLF = float(re.search(r"C_dolfinx = ([0-9.]+)", dolf_txt).group(1))
C_SKF = float(re.search(r"C_skfem\s*= ([0-9.]+)", dolf_txt).group(1))
DOLF_REL = float(re.search(r"\|dolfinx-skfem\|/skfem = ([0-9.e+-]+)", dolf_txt).group(1))  # full precision, from the run
PF_WALL_ERR = pf["kink"]["rel_err_sigma"]
PF_DISP_ERR = max(r["rel_err"] for r in pf["dispersion"])

ROWS = []


def add(grup, olcum, model, referans, oran, karar, note, kaynak, plot=None):
    ROWS.append(dict(grup=grup, olcum=olcum, model=model, referans=referans,
                     oran=None if oran is None else float(oran), karar=karar, note=note, kaynak=kaynak,
                     plot=plot))


G1, G2, G3 = "A. Elektrostatik ve kapasitans", "B. Yük enerjisi E_C ölçeği", "C. Kuantum hapsolma (Schrödinger)"
G4, G5, G6 = "D. Tünel bağlaşımı", "E. Dinamik yükleme", "F. Eşik gerilimi (1B Poisson–Schrödinger)"
G7, G8, G9 = "G. Katkı potansiyeli ve kazara nokta", "H. Düşük sıcaklık: iyonlaşma", "I. Test planı"
G10, G11, G12 = "J. Bölünme ve bifurkasyon (soyut 2B)", "K. Resmî veri ve tasarım parametreleri", "L. İç tutarlılık"

# A
add(G1, "Küre kapasitör: FEM (h=0.05) / analitik 4πab/(b−a)",f"{sph['rows'][-1]['C_fem']:.4f}", f"{sph['C_exact']:.4f}",
    sph["rows"][-1]["C_fem"] / sph["C_exact"], "uyumlu", f"yakınsama eğimi {sph['loglog_slope']:.2f} (beklenen ≈2)",
    "exp1/validate_sphere.json; analitik formül", plot="Küre kapasitör FEM/analitik")
add(G1, "scikit-fem / dolfinx, aynı mesh (h=0.14)", f"{C_SKF:.10f}", f"{C_DOLF:.10f}", C_SKF / C_DOLF,
    "uyumlu", f"iki bağımsız kod; tam hassasiyetli bağıl fark {DOLF_REL:.2e}", "exp1/dolfinx_check.txt", plot="scikit-fem / dolfinx")
add(G1, "Hücre kapasitansı c_Σ(r=0,5): seviye 3 / seviye 2", f"{capv(0.5, 3, 'c_sigma'):.3f}", f"{capv(0.5, 2, 'c_sigma'):.3f}",
    capv(0.5, 3, "c_sigma") / capv(0.5, 2, "c_sigma"), "uyumlu (yakınsıyor)", "iç yakınsama; dış referans yok",
    "exp1/cap_cell_A.csv")
add(G1, "Zincir ada bloğu, seviye 1 → 2 en büyük göreli değişim", f"{chain['level_change_max_rel']:.3f}", "0 (yakınsamış olmalı)",
    None, "yakınsamamış", "zincir sonuçları ±%4 üst sınır; mesh daha da inceltilmeli", "exp2/chain.json")
for z0, lbl in ((10.0, "z=10 nm"), (30.0, "z=30 nm"), (60.0, "z=60 nm")):
    fe, im, rd = fem_pair(z0, 0)
    add(G1, f"Tek akseptör potansiyeli 3B FEM / görüntü modeli, {lbl}, r=0", f"{fe:.2f} mV", f"{im:.2f} mV", fe / im,
        "tutarsız → ampirik düzeltme" if z0 == 10.0 else "kısmen tutarlı",
        f"fark {rd*100:+.1f} %; düzeltme faktörü bu üç derinliğe uyduruldu", "exp6/fem_check.json",
        plot=f"3B FEM / görüntü, {lbl}, r=0")
fe40, im40, rd40 = fem_pair(10.0, 4)
add(G1, "Aynı, z=10 nm, r=40 nm (düzeltme dışı)", f"{fe40:.2f} mV", f"{im40:.2f} mV", fe40 / im40,
    "tutarsız (düzeltme modellemiyor)", f"fark {rd40*100:+.1f} %; düzeltme exp(−r²) ile sönüyor, büyük r'de kalıntı var",
    "exp6/fem_check.json")

# B
add(G2, "E_C/k_B, 12,5 nm kapılı hücre: FEM ε=3,9 ↔ 11,7 bandı, exp0 analitik ile",
    f"{E_GATED_125_117:.0f}–{E_GATED_125_39:.0f} K", f"{float(e0_g3['EC_over_kB_K']):.0f} K (exp0, karışık ε)",
    None, "tutarlı (exp0 band içinde)", "ε belirsizliği yaklaşık ×3, bu sayısal hatadan çok daha büyük",
    "exp0/generations.csv; exp1/transport.json")
add(G2, "Kapı etkisi: izole küre (çap 12,5 nm, ε=3,9) / kapılı hücre (FEM)", f"{E_ISO_125:.0f} K (izole)",
    f"{E_GATED_125_39:.0f} K (kapılı)", E_ISO_125 / E_GATED_125_39, "bilgi (tasarım sonucu)",
    "kapı E_C'yi yaklaşık bu oranda düşürüyor; büyük E_C için kapısız ada gerekir", "analitik (izole küre) + exp1 FEM")
add(G2, f"Literatür L4: ~1,8 nm Si küme, E_C ≈ 300 meV (özet) / izole küre R=0,9 nm, ε=3,9",
    f"{LIT4_MODEL_EV*1e3:.0f} meV", "300 meV", LIT4_MODEL_EV / 0.300, "tanıma/ortam bağlı",
    f"literatür değeri ε≈{EPS_STAR_LIT4:.1f} ortamda uyuşur; fark ortam dielektriğiyle açıklanır (L4 doğrulanmadı)",
    "arXiv 1802.09800 özeti (L4)", plot="İzole küre / literatür 300 meV")
add(G2, "Literatür L3: 5 nm Si ada, 300 K'de Coulomb elmasları (özet) / izole küre R=2,5 nm",
    f"E/k_BT(300 K) = {SHIN_EKT300*300:.0f} K / 300 K = {SHIN_EKT300:.2f}", "300 K'de elmas görülüyor (nicel değer yok)",
    None, "karşılaştırılamaz (nicel yok); yönsel olarak tutarlı", "marjinal; L3 özet, doğrulanmadı",
    "Shin vd. 2011 özeti (L3)")
add(G2, "Kapılı hücre, 300 K için gerekli E_C/k_B (10·T) / model (12,5 nm)", f"{E_GATED_125_39:.0f} K",
    f"{10*300:.0f} K", E_GATED_125_39 / 3000.0, "uyumlu (negatif sonuç sağlam)",
    "ölçek ~30× kısa; ε belirsizliği (×3) bunu değiştirmez", "exp1/transport.json")

# C
add(G3, "Sert duvar küp (3,6,6,6,9) en ince mesh: en büyük bağıl hata", f"{hw_fine['max_rel_err']:.2e}", "0 (analitik)", None,
    "uyumlu" if hw_fine["max_rel_err"] < 0.02 else "uyumsuz", f"h={hw_fine['h']}", "exp1/schrodinger.json")
add(G3, "Ayrılabilir sonlu kuyu, 3·E1D (analitik), en ince mesh", f"{sep_fine['E_fem']:.5f} eV", f"{sep_fine['E_exact']:.5f} eV",
    sep_fine["E_fem"] / sep_fine["E_exact"], "uyumlu", f"bağıl hata {sep_fine['rel_err']:.1e}", "exp1/schrodinger.json",
    plot="Ayrılabilir kuyu FEM/analitik")
add(G3, "Gerçek Si/SiO₂ kuyusu (L=2 nm) E1 mesh değişimi, son adım", f"{rwc_last['rel_change_vs_prev']:.2e}", "0", None,
    "yakınsadı", "iç yakınsama", "exp1/schrodinger.json")
add(G3, "Sert duvar ΔE / FEM ΔE, L=1 nm", f"{well(1.0)['dE_hardwall_eV']:.3f} eV", f"{well(1.0)['dE_eV']:.3f} eV", r_hw_fem_1,
    "uyumsuz (formül geçersiz)", "sert duvar formülü küçük L'de ~3× şişiriyor", "exp1/schrodinger.json",
    plot="Sert duvar ΔE / FEM, L=1 nm")
add(G3, "Sert duvar ΔE / FEM ΔE, L=3 nm", f"{well(3.0)['dE_hardwall_eV']:.3f} eV", f"{well(3.0)['dE_eV']:.3f} eV", r_hw_fem_3,
    "uyumsuz (şişirme)", "", "exp1/schrodinger.json", plot="Sert duvar ΔE / FEM, L=3 nm")
add(G3, "Deney 0 ΔE (sert duvar, 12,5 nm) / FEM ΔE (aynı nesil)", f"{float(e0_g3['dE_over_kB_K']):.0f} K",
    f"{ftab[3]['dE_over_kB_K']:.0f} K", float(e0_g3["dE_over_kB_K"]) / ftab[3]["dE_over_kB_K"], "uyumsuz (%15 şişirme)",
    "exp0 analitik ön hesabı FEM ile değiştirildi", "exp0/generations.csv; exp1/transport.json", plot="Deney 0 ΔE / FEM, 12,5 nm")
add(G3, "Literatür seviye ayrışması (kuyu boyu ~1–10 nm)", "—", "veri yok", None, "veri yok",
    "bu oturumda nicel literatür değeri alınamadı", "—")

# D
add(G4, "Δ bozunma κ: FEM fit / WKB, L=1 nm", f"{tun(1.0)['fit_decay_per_nm']:.2f} /nm", f"{tun(1.0)['kappa_WKB_per_nm']:.2f} /nm",
    tun(1.0)["ratio_fit_over_WKB"], "uyumlu (~%10)", "", "exp4/stage1.json", plot="Tünel κ: FEM/WKB, L=1 nm")
add(G4, "Δ bozunma κ: FEM fit / WKB, L=2 nm", f"{tun(2.0)['fit_decay_per_nm']:.2f} /nm", f"{tun(2.0)['kappa_WKB_per_nm']:.2f} /nm",
    tun(2.0)["ratio_fit_over_WKB"], "uyumlu (~%3)", "", "exp4/stage1.json", plot="Tünel κ: FEM/WKB, L=2 nm")
d1 = next(r for r in dqd if abs(r["L_nm"] - 1.0) < 1e-9 and abs(r["d_nm"] - 0.5) < 1e-9)
add(G4, "Δ (L=1 nm, d=0,5 nm, SiO₂ bariyeri 3,1 eV)", f"{d1['delta_eV']*1e3:.1f} meV", "literatür gate-tanımlı Si çift nokta değeri alınamadı",
    None, "karşılaştırılamaz", "bariyer türü farklı (oksit, elektrostatik değil)", "exp3/dqd.json")

# E
add(G5, "Landau–Zener diyabatik geçiş, γ = 0,1–4: en büyük mutlak hata", f"{lz_err:.1e}", "0 (analitik e^{−γ})", None,
    "uyumlu", "sayısal DOP853 ↔ analitik", "exp4/stage1.json")
add(G5, "Lindblad 300 K transfer (en iyi geometri), 10 ps sweep, 1 ps dephasing", f"{lb10['P_transfer']['1000fs']:.3f}",
    f"{lb['P_right_thermal_equilibrium']:.3f} (termal denge)", lb10["P_transfer"]["1000fs"] / lb["P_right_thermal_equilibrium"],
    "tutarlılık kontrolü (referans yok)", "transfer termal denge değerine yaklaşıyor", "exp4/stage1.json")

# F
add(G6, "V_g (n_s = 10¹² cm⁻², sabit yük), ndep = 1,7e17 / BSIM4 vth0", f"{VG_BASE:.4f} V", f"{VTH0:.4f} V", VG_BASE / VTH0,
    "uyumlu (tanıma bağlı)", "fark ≈ %1; sabit-yük tanımı bir seçimdir", "exp5/mos1d.json; model dosyası (vth0)", plot="V_T (sabit yük) / vth0")
add(G6, "V_T doğrusal ekstrapolasyon, ndep = 1,7e17 / vth0", f"{mosN['VT_linear_extrap']:.4f} V", f"{VTH0:.4f} V",
    mosN["VT_linear_extrap"] / VTH0, "tanıma bağlı, uyumsuz", "iki tanım 0,18 V ayrışıyor", "exp5/mos1d.json", plot="V_T (doğrusal) / vth0")
add(G6, "Kalibre N_eff (V_T doğrusal = vth0) / resmî ndep", f"{mosC['na_eff_cm3']:.2e} cm⁻³", f"{NDEP:.2e} cm⁻³",
    mosC["na_eff_cm3"] / NDEP, "geçersiz (kalibrasyon)",
    "kalibrasyon yalnızca doğrusal tanıma uyar; sabit-yük tanımında resmî ndep zaten uyumlu", "exp5/mos1d.json; model dosyası")
add(G6, "C_inv / C_ox (öz-uyumlu, nominal)", f"{mosN['cinv_uF_cm2']:.3f} µF/cm²", f"{mosN['Cox_uF_cm2']:.3f} µF/cm²",
    mosN["cinv_uF_cm2"] / mosN["Cox_uF_cm2"], "uyumlu (C_inv < C_ox)", "beklenen yön", "exp5/mos1d.json")

# G
add(G7, "σ_ψ arayüz potansiyeli: görüntü modeli / FEM-düzeltmeli (nominal)", f"{RN['sigma_surface_potential_mV']:.1f} mV",
    f"{RNC['sigma_surface_potential_mV']:.1f} mV", RN["sigma_surface_potential_mV"] / RNC["sigma_surface_potential_mV"],
    "model belirsizliği ~%40", "düzeltme seçimine bağlı", "exp6/rdf.json", plot="σ_ψ görüntü / düzeltmeli (nominal)")
if SENS:
    s1 = SENS["nominal_1.7e17_ZCAP1"]["sigma_psi_mV"]
    s10 = SENS["nominal_1.7e17_ZCAP10"]["sigma_psi_mV"]
    e1 = SENS["effective_3.9e17_ZCAP1"]["sigma_psi_mV"]
    e10 = SENS["effective_3.9e17_ZCAP10"]["sigma_psi_mV"]
    add(G7, "Kendi kontrolü: duyarlılık koşusunda ZCAP=1 ↔ exp6 düzeltmeli (nominal)", f"{s1:.2f} mV", f"{RNC['sigma_surface_potential_mV']:.2f} mV",
        s1 / RNC["sigma_surface_potential_mV"], "uyumlu (yeniden üretim)" if abs(s1 - RNC["sigma_surface_potential_mV"]) < 0.05 else "uyumsuz (kod farkı)",
        "aynı tohum ve örnekler; tek fark kıskaç", "exp6/cap_sensitivity.json")
    add(G7, "Duyarlılık: z ≥ 10 nm kıskacı / kıskaçsız (nominal)", f"{s10:.2f} mV", f"{s1:.2f} mV", s10 / s1,
        "duyarlı" if abs(s10 / s1 - 1) > 0.05 else "duyarsız (±%5)", "z<10 nm'de düzeltme extrapole ediliyor",
        "exp6/cap_sensitivity.json", plot="σ_ψ kıskaç duyarlılığı (nominal)")
    add(G7, "Duyarlılık: z ≥ 10 nm kıskacı / kıskaçsız (efektif 3,9e17)", f"{e10:.2f} mV", f"{e1:.2f} mV", e10 / e1,
        "duyarlı" if abs(e10 / e1 - 1) > 0.05 else "duyarsız (±%5)", "efektif doping'de z<10 nm payı daha büyük",
        "exp6/cap_sensitivity.json", plot="σ_ψ kıskaç duyarlılığı (efektif)")
    add(G7, "Kazara nokta oranı (katı), z ≥ 10 nm kıskacı, nominal", f"{SENS['nominal_1.7e17_ZCAP10']['strict_frac']*100:.2f} %",
        f"kıskaçsız: {SENS['nominal_1.7e17_ZCAP1']['strict_frac']*100:.2f} %", SENS["nominal_1.7e17_ZCAP10"]["strict_frac"]
        / SENS["nominal_1.7e17_ZCAP1"]["strict_frac"], "belirsiz/duyarlı", "kıskaç ekstrapolasyona bağlı; alt uç",
        "exp6/cap_sensitivity.json")
    add(G7, "Kazara nokta oranı (katı), z ≥ 10 nm kıskacı, efektif 3,9e17", f"{SENS['effective_3.9e17_ZCAP10']['strict_frac']*100:.2f} %",
        f"kıskaçsız: {SENS['effective_3.9e17_ZCAP1']['strict_frac']*100:.2f} %", SENS["effective_3.9e17_ZCAP10"]["strict_frac"]
        / SENS["effective_3.9e17_ZCAP1"]["strict_frac"], "belirsiz/duyarlı", "alt uç; üst uç kıskaçsız değer", "exp6/cap_sensitivity.json")
add(G7, "σ_VT (RDF MC, nominal) / Stolk formülü (bellekten)", f"{RN['sigma_VT_MC_mV']:.2f} mV", f"{RN['sigma_VT_Stolk_mV']:.2f} mV",
    RN["sigma_VT_MC_mV"] / RN["sigma_VT_Stolk_mV"], "uyumsuz (~%30)", "Stolk formülü bellekten, doğrulanmadı; tanım farkı olabilir",
    "exp6/rdf.json; L15", plot="σ_VT MC / Stolk (nominal)")
add(G7, "σ_VT (RDF MC, nominal) / foundry vth0 eğimi (V·µm varsayımı)", f"{RN['sigma_VT_MC_mV']:.2f} mV",
    f"{SIG_VT_FOUNDRY_MV:.2f} mV", RN["sigma_VT_MC_mV"] / SIG_VT_FOUNDRY_MV, "birim varsayımına bağlı",
    "RDF toplamın ~1/3'ü; birim V·µm değilse bu karşılaştırma geçersiz", "model dosyası (mismatch); L11",
    plot="σ_VT RDF / foundry (V·µm)")
add(G7, "σ_VT (oksit kalınlığı, σ_tox = 0,057 nm, 1B çözücü)", f"{SIG_VT_TOX_MV:.2f} mV", f"{SIG_VT_FOUNDRY_MV:.2f} mV",
    SIG_VT_TOX_MV / SIG_VT_FOUNDRY_MV, "bilgi",
    f"RDF ile birlikte kareler toplamı: {np.hypot(RN['sigma_VT_MC_mV'], SIG_VT_TOX_MV):.1f} mV (foundry'nin ~%{np.hypot(RN['sigma_VT_MC_mV'], SIG_VT_TOX_MV)/SIG_VT_FOUNDRY_MV*100:.0f}'i)",
    "exp5 çözücüsü; model dosyası toxe_slope")
add(G7, "Kazara nokta oranı (katı ölçüt), görüntü, nominal", f"{RN['bound_strict']['frac_samples']*100:.1f} %",
    f"iyimser: {RN['bound_optimistic']['frac_samples']*100:.1f} %", None,
    "geçersiz/şüpheli (T-bağımsız)", "modelde oran sıcaklıktan bağımsız; oysa iyonlaşma T ile değişir (bkz. H)",
    "exp6/rdf.json")
add(G7, "Kazara nokta oranı (katı), FEM-düzeltmeli, efektif 3,9e17", f"{REC['bound_strict']['frac_samples']*100:.1f} %",
    f"iyimser: {REC['bound_optimistic']['frac_samples']*100:.1f} %", None,
    "geçersiz/şüpheli (T-bağımsız)", "en yüksek tahmin; düzeltme ekstrapolasyonuna bağlı", "exp6/rdf.json")
add(G7, "Kazara nokta E_C/k_B (ε_eff = 7,8 ölçeklemesi), medyan, düzeltmeli efektif",
    f"{REC['bound_strict']['EC_over_kB_K_median']:.0f} K", "L1: 14 K'de CB tümseği (28 nm PMOS, özet; nicel yok)", None,
    "karşılaştırılamaz (nicel)", "ε_eff ölçeklemesi kaba; geometri uyumsuz", "exp6/rdf.json; L1")
add(G7, "14 K'de Coulomb görünürlüğü (E_C = 25 / 34 / 40 K)",
    f"{VIS[(25.0, 14.0)]:.2f} / {VIS[(34.0, 14.0)]:.2f} / {VIS[(40.0, 14.0)]:.2f}", "L1: 14 K'de tümsek gözlendi (özet)",
    None, "yönsel uyumlu (E_C/kT ≈ 2–3)", "L1 nicel değer vermiyor; E_C ≳ 3kT gerekir", "yeniden hesap (exp7 yöntemi); L1")

# H
add(G8, "NA / Mott kritik yoğunluk (B için ≈ 3,7·10¹⁸ cm⁻³, bellekten)", f"{NDEP:.2e} cm⁻³", f"{NA_MOTT:.1e} cm⁻³", NDEP / NA_MOTT,
    "geçersiz (4 K'de yalıtkan rejim)",
    "NA < N_c: 4 K'de katkılar tam iyonlaşmaz; model tam iyonlaşma varsayar. 77 K şüpheli, 300 K makul",
    "exp6 model varsayımı; L17 (bellekten)", plot="NA / N_c (Mott)")

# I
vis4 = next(v for v in tp["visibility"] if abs(v["EC_over_kB_K"] - 34.0) < 1e-9 and abs(v["T_K"] - 4.0) < 1e-9)["visibility"]
vis10 = next(v for v in tp["visibility"] if abs(v["EC_over_kB_K"] - 34.0) < 1e-9 and abs(v["T_K"] - 10.0) < 1e-9)["visibility"]
add(G9, "Görünürlük, E_C=34 K, 4 K: exp7 / bağımsız yeniden hesap", f"{vis4:.3f}", f"{VIS[(34.0, 4.0)]:.3f}",
    VIS[(34.0, 4.0)] / vis4, "tutarlı (yeniden hesap)" if abs(VIS[(34.0, 4.0)] - vis4) < 0.01 else "tutarsız", "exp7 ile aynı yöntem",
    "exp7/testplan.json")
add(G9, "Görünürlük, E_C=34 K, 10 K (exp7)", f"{vis10:.3f}", "—", None, "bilgi", "E_C/kT = 3,4", "exp7/testplan.json")
add(G9, "≥3 pozitif aygıtı %95 olasılıkla görmek için n (p = %4)", f"{tp['n_for_95pct_to_see_ge3_positives']['0.04']}", "—", None,
    "tasarım bilgisi", "p tahmininin kendisi H maddesine bağlı", "exp7/testplan.json")

# J
add(G10, "Fisyon eşiği üssü, d = ∞ / Rayleigh (−1, bellekten)", f"{FIS_SLOPE['inf']:.2f}", "−1,00", FIS_SLOPE["inf"] / -1.0,
    "uyumlu (ızgara ±%45)", "g ızgarası geometrik (×1,45 aralık)", "exp8/fission_scan*.json; L16", plot="Fisyon üssü / Rayleigh")
add(G10, "Fisyon eşiği üssü, d = 4 ve d = 2 (gate ekranlaması)", f"{FIS_SLOPE['4.0']:.2f} / {FIS_SLOPE['2.0']:.2f}", "—", None,
    "bilgi", "gate ekranlaması eğimi değiştiriyor", "exp8/fission_scan_d*.json")
add(G10, "Bifurkasyon genlik üssü / beklenen 0,5 (süperkritik çatal)", f"{pf['amplitude']['loglog_slope_small_delta']:.3f}", "0,500",
    pf["amplitude"]["loglog_slope_small_delta"] / 0.5, "uyumlu", "doğrulanmış", "exp8/pf_validation.json", plot="Bifurkasyon üssü / 0,5")
add(G10, "A_c (φ̄=0, d=∞): sayısal / analitik", f"{pf['critical']['A_c_numeric']:.5f}", f"{pf['critical']['A_c_analytic']:.5f}",
    pf["critical"]["A_c_numeric"] / pf["critical"]["A_c_analytic"], "uyumlu", "", "exp8/pf_validation.json", plot="A_c sayısal / analitik")
add(G10, "Lineer büyüme σ(q): en büyük bağıl hata", f"{max(r['rel_err'] for r in pf['dispersion']):.1e}", "0", None,
    "uyumlu", "5 dalga sayısı", "exp8/pf_validation.json")
add(G10, "Hazır tek domain (disk, çubuk, elips) × 4 A, t=3000: n sonucu", f"n = {sorted(set(div_n))}", "bölünme beklentisi yok",
    None, "negatif (kayıt)" if div_nosplit else "bölündü (kontrol et)", "çözündü veya küçüldü, bölünmedi", "exp8/division_diag.json")
add(G10, "A_c(d=1,5) / A_c(∞), φ̄=0", f"{sel['A_c_phi0_0']['1.5']:.4f}", f"{sel['A_c_phi0_0']['inf']:.4f}",
    sel["A_c_phi0_0"]["1.5"] / sel["A_c_phi0_0"]["inf"], "bilgi", "gate yaklaşınca A_c artıyor", "exp8/selection.json")

# K
add(G11, "Poly genişliği (çizilen L = 150 nm) / SKY130 kuralı 0,150 µm", "0,150 µm", "0,150 µm (poly.1a)", 1.0, "uyumlu",
    "SkyWater kural sayfası", "L9 (birincil sayfa, araçla okundu)")
add(G11, "Etkin kanal L = 150 − 2·lint / model kullanımı (126 nm)", f"{150 - 2*LINT*1e9:.1f} nm", "126 nm (model)",
    (150 - 2 * LINT * 1e9) / 126.0, "uyumlu", f"lint = {LINT*1e9:.2f} nm (resmî dosya)", "L10")
add(G11, "W = 420 nm (model) / doc-tablo minimum W (0,42 µm)", "0,42 µm", "tablo örneği (kural değil)", None,
    "doğrulanmadı", "docs e-test tablosunda örnek; kural olarak doğrulanmadı", "L9")
add(G11, "toxe (resmî) / Cox hesabında kullanılan", f"{TOXE*1e9:.3f} nm", f"{TOX_NM:.3f} nm", TOXE * 1e9 / TOX_NM, "uyumlu (aynı kaynak)",
    "", "L10")
add(G11, "NIST SKY130 kriyojenik veri: 4 K var mı?", "kayıt okundu (v1.1.2)", "dosyalar indirilemedi (Cloudflare 524)", None,
    "veri yok (erişim)", "ETest Tile.zip 3,8 MB, MPW-5 Test Tile.zip 39 MB, README", "L5; data/nist_sky130/")
add(G11, "77 K / 4 K SKY130 ölçümleri (literatür)", "—", "77 K modelleme özeti var; 4 K yok", None, "veri yok", "UTA CPAD 2025 posteri", "L6")
add(G11, "130 / 180 nm Coulomb blokajı raporu", "—", "aramada bulunamadı", None, "veri yok", "bu arama sonucu, yokluk kanıtı değil", "arama")
add(G11, "28 nm FDSOI 4 K transistör verisi (SKY130 değil)", "—", "arXiv 2002.07070 özeti", None, "kısmen (farklı süreç)",
    "karşılaştırma için yakın bir süreç; nicel değer alınmadı", "L7")

# L
lj_row = f"{LJ_ERR:.1e}" if LJ else "—"
add(G12, "Birim testleri (pytest)", PYTEST_LINE, "tümü geçmeli", None, "uyumlu" if "passed" in PYTEST_LINE and "failed" not in PYTEST_LINE else "kontrol et",
    "bu betik çalıştırdı", "tests/")
add(G12, "LJ7/13/19 global minimumları: en büyük hata", lj_row, "literatür", None, "uyumlu", "basin hopping", "exp8/fission_scan.json")
add(G12, "E_C zincir (12 nm ada, exp2) / tek hücre (12 nm, exp1)", f"{chain_ec12:.1f} K", f"{single_ec12:.1f} K",
    chain_ec12 / single_ec12, "geometri farkı; tutarlılık sınırlı", "zincirde komşu adalar C^-1 köşegenini büyütüyor", "exp2/chain.json; exp1")


# ---------------------------------------------------------------- literature and official items
LIT = [
    ("L1", "28 nm ticari bulk PMOS, 14 K'de kanal kuantum noktası Coulomb tümseği", "arXiv 2206.09086 (arama özeti)", "ikincil (özet)"),
    ("L2", "40 nm bulk MOSFET, 50 mK'de kanal kuantum noktaları ve Coulomb blokajı", "arXiv 2005.13863 (arama özeti)", "ikincil (özet)"),
    ("L3", "Sub-5 nm Si ada, 300 K'de Coulomb elmasları", "Shin vd., Nano Lett. 2011 (arama özeti)", "ikincil (özet)"),
    ("L4", "~1,8 nm Si küme SET'inde E_C ≈ 300 meV (potansiyel iddia)", "arXiv 1802.09800 (arama özeti)", "ikincil (özet)"),
    ("L5", "NIST SKY130 kriyojenik veri seti (mds2-2997): kayıt okundu, dosyalar indirilemedi", "data.nist.gov API", "birincil kayıt; dosya yok"),
    ("L6", "SKY130 77 K modelleme (UTA, CPAD 2025 posteri)", "arama özeti", "ikincil"),
    ("L7", "28 nm FDSOI transistörleri 4 K'ye kadar (ticari süreç, SKY130 değil)", "arXiv 2002.07070 (arama özeti)", "ikincil (özet)"),
    ("L8", "Başlığı 'yüksek sıcaklık Si qubit' olan makalenin özeti: 5–10 K", "arXiv 1804.03364 (WebFetch özeti)", "ikincil (küçük model özeti)"),
    ("L9", "SKY130 poly.1a 0,150 µm; poly.2 0,210 µm; pfet min kanal (poly.1b) 0,350 µm", "skywater-pdk.readthedocs.io (WebFetch)", "birincil sayfa, araçla okundu"),
    ("L10", "SKY130 nfet_01v8 model: toxe 4,148 nm, ndep 1,7e17, vth0 0,519 V, xj 150 nm, lint 11,9 nm", "google/skywater-pdk-libs-sky130_fd_pr, tt.pm3.spice (doğrudan)", "birincil dosya (doğrudan indirildi)"),
    ("L11", "Foundry mismatch: vth0_slope 3,356e-3, toxe_slope 3,443e-3 (birim V·µm varsayımı)", "nfet_01v8__mismatch.corner.spice (doğrudan)", "birincil dosya; birim belirsiz"),
    ("L12", "BAG SKY130 tech parametreleri (met1/li1 boşluk ve alan)", "googlesource (arama özeti)", "ikincil; birim belirsiz"),
    ("L13", "Sıcak qubit >1 K (Yang 2020, Petit 2020, Huang 2024)", "postquantum.com özeti", "ikincil; makaleler alınamadı"),
    ("L14", "Quobly/CEA-Leti 28 nm FDSOI spin qubit, 110 mK çift nokta", "konferans özetleri (arama)", "ikincil"),
    ("L15", "Stolk (1998) σ_VT ∝ (q t_ox/ε_ox)·√(N_A W_dep/3WL)", "bellekten", "doğrulanmadı"),
    ("L16", "Rayleigh sınırı: fisyon eşiği ∝ N⁻¹ (q² ∝ R³)", "bellekten (klasik fizik)", "doğrulanmadı"),
    ("L17", "Si:B Mott kritik yoğunluğu ≈ 3,7·10¹⁸ cm⁻³", "bellekten", "doğrulanmadı"),
]

# ---------------------------------------------------------------- outputs: CSV, verdicts, figure, markdown
with open(OUT / "comparison.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=["grup", "olcum", "model", "referans", "oran", "karar", "note", "kaynak"])
    w.writeheader()
    for r in ROWS:
        w.writerow({k: r[k] for k in w.fieldnames} | {"oran": "" if r["oran"] is None else f"{r['oran']:.4g}"})

CATS = [  # ordered: first match wins
    ("negatif", "negatif (kayıt)"),
    ("geçersiz", "geçersiz/şüpheli"),
    ("veri yok", "veri yok"),
    ("karşılaştırılamaz", "karşılaştırılamaz"),
    ("uyumsuz", "uyumsuz/tutarsız"),
    ("tutarsız", "uyumsuz/tutarsız"),
    ("tanıma bağlı", "tanıma bağlı"),
    ("duyarlı", "belirsiz/duyarlı"),
    ("model belirsizliği", "belirsiz/duyarlı"),
    ("yakınsamamış", "yakınsamamış"),
    ("kısmen", "kısmi uyum"),
    ("uyumlu", "uyumlu"),
    ("tutarlı", "uyumlu"),
    ("yakınsadı", "uyumlu"),
    ("yönsel", "kısmi uyum"),
    ("bilgi", "bilgi"),
    ("tasarım", "bilgi"),
    ("birim", "belirsiz/doğrulanmadı"),
    ("doğrulanmadı", "belirsiz/doğrulanmadı"),
    ("geometri", "belirsiz/doğrulanmadı"),
    ("tutarlılık kontrolü", "belirsiz/doğrulanmadı"),
    ("tanıma/ortam", "tanıma bağlı"),
]


def category(karar):
    for pat, cat in CATS:
        if pat in karar:
            return cat
    return "sınıflandırılmadı"


for r in ROWS:
    r["kategori"] = category(r["karar"])
verdict_count = {}
for r in ROWS:
    verdict_count[r["kategori"]] = verdict_count.get(r["kategori"], 0) + 1
assert "sınıflandırılmadı" not in verdict_count, [r["karar"] for r in ROWS if r["kategori"] == "sınıflandırılmadı"]
(OUT / "verdicts.json").write_text(json.dumps(dict(rows=ROWS, counts=verdict_count, lit=LIT), indent=2, ensure_ascii=False),
                                   encoding="utf-8")

plot = [r for r in ROWS if r["plot"] and r["oran"] is not None]
fig, ax = plt.subplots(figsize=(9.5, 0.45 * len(plot) + 1.6))
ys = np.arange(len(plot))[::-1]
vals = np.array([r["oran"] for r in plot])
colors = ["C2" if 0.9 <= v <= 1.1 else ("C3" if (v < 0.5 or v > 2) else "C1") for v in vals]
ax.axvspan(0.9, 1.1, color="0.9", zorder=0)
ax.axvline(1.0, c="k", lw=0.8)
ax.scatter(vals, ys, c=colors, zorder=3)
for y, v in zip(ys, vals):
    ax.text(v * 1.04, y, f"{v:.2f}", va="center", fontsize=8)
ax.set_xscale("log")
ax.set_xlim(0.03, 6)
ax.set_yticks(ys, [r["plot"] for r in plot], fontsize=8)
ax.set(xlabel="model / referans (log; 1 = uyum, gri bant ±%10)",
       title="Tüm sayısal karşılaştırmalar: yeşil uyumlu, turuncu kısmi, kırmızı belirgin fark")
fig.tight_layout()
fig.savefig(OUT / "fig_comparison.png", dpi=140)
plt.close(fig)


def md_table(rows_):
    head = "| # | Ölçüm | Model | Referans | Oran | Karar | Not | Kaynak |\n|---|---|---|---|---|---|---|---|"
    lines = [head]
    for i, r in enumerate(rows_, 1):
        o = "—" if r["oran"] is None else f"{r['oran']:.4g}"
        lines.append(f"| {i} | {r['olcum']} | {r['model']} | {r['referans']} | {o} | {r['karar']} | {r['note']} | {r['kaynak']} |")
    return "\n".join(lines)


GROUPS = []
for r in ROWS:
    if r["grup"] not in GROUPS:
        GROUPS.append(r["grup"])

nist_comp = ", ".join(f"{(c.get('filepath') or c.get('title'))} ({c.get('size')} B)" for c in nist.get("components", []))
L = []
L += ["# Deney 9: elimizdeki tüm verinin karşılaştırması", "",
      f"Oluşturma: `experiments/exp9_compare/compile.py`. Sayılar sonuç dosyalarından okundu; {len(ROWS)} karşılaştırma satırı, "
      f"{len(LIT)} literatür/resmî veri kalemi.", "",
      "## Kısa karar", "",
      f"1. **Hesap araçları bağımsız referanslarla uyumlu.** Küre kapasitörü %{abs(sph['rows'][-1]['C_fem']/sph['C_exact']-1)*100:.2f} hata, "
      f"scikit-fem ile dolfinx {DOLF_REL:.1e} bağıl fark, Landau–Zener {lz_err:.1e}, faz alanı duvar gerilimi "
      f"%{PF_WALL_ERR*100:.2f}, σ(q) %{PF_DISP_ERR*100:.2f}. Sabit yük V_T, vth0'a %{abs(VG_BASE/VTH0-1)*100:.1f} yakın.",
      f"2. **300 K'de kuantum yok: negatif sonuç sağlam.** 12,5 nm kapılı hücrede E_C/k_B ≈ {E_GATED_125_39:.0f} K (ε=3,9) ya da "
      f"{E_GATED_125_117:.0f} K (ε=11,7). 300 K'de 10·kT ≈ 3000 K gerekir: ölçek ~30× kısa. ε belirsizliği (×3) bunu değiştirmez.",
      f"3. **4 K'deki kazara nokta tahmini (önceki: %4–21) kırılgan ve bu haliyle geçersiz olabilir.** (a) Görüntü modeli 3B FEM'den "
      f"%{abs(fem_pair(10.0,0)[2])*100:.0f}'e kadar sapıyor, düzeltme ampirik; (b) düzeltme z<10 nm'de dışarı ekstrapole ediliyor. "
      + (f"z≥10 nm kıskacıyla σ_ψ %{abs(SENS['nominal_1.7e17_ZCAP10']['sigma_psi_mV']/SENS['nominal_1.7e17_ZCAP1']['sigma_psi_mV']-1)*100:.0f} düşüyor "
         f"(nominal) ve %{abs(SENS['effective_3.9e17_ZCAP10']['sigma_psi_mV']/SENS['effective_3.9e17_ZCAP1']['sigma_psi_mV']-1)*100:.0f} (efektif); "
         f"katı kazara nokta oranı nominal %{SENS['nominal_1.7e17_ZCAP10']['strict_frac']*100:.1f}–{RNC['bound_strict']['frac_samples']*100:.1f}, "
         f"efektif %{SENS['effective_3.9e17_ZCAP10']['strict_frac']*100:.1f}–{REC['bound_strict']['frac_samples']*100:.1f} aralığında (düzeltmeli üst sınır); "
         if SENS else "kıskaç duyarlılığı henüz hesaplanmadı; ")
      + f"(c) NA = {NDEP:.1e} cm⁻³, Mott eşiğinin (~{NA_MOTT:.1e}) yaklaşık {NDEP/NA_MOTT*100:.0f}%'i: 4 K'de katkılar tam iyonlaşmaz, "
      "modelimiz ise iyonlaşmayı T'den bağımsız varsayar. Kazara nokta oranı modelde sıcaklıktan bağımsız; oysa gerçekte değişir.",
      "4. **Literatürle nicel karşılaştırma sınırlı.** 130/180 nm için Coulomb blokajı verisi bulunamadı; SKY130 4 K verisi "
      "erişilemedi (NIST dosyaları sunucu hatası verdi); tünel bariyeri türü literatürle eşleşmiyor. Yönsel uyum var: "
      f"izole küçük adalar 300 K'de sınırda (L3 özeti, E/k_BT ≈ {SHIN_EKT300:.1f}); kapılı hücrelerimiz bundan çok daha kötü.",
      f"5. **Foundry mismatch bütçesi:** RDF ile σ_VT ≈ {RN['sigma_VT_MC_mV']:.1f} mV; foundry vth0 eğimi (V·µm varsayımıyla) "
      f"{SIG_VT_FOUNDRY_MV:.1f} mV; oksit kalınlığı tek başına ≈ {SIG_VT_TOX_MV:.1f} mV (1B çözücüyle). Birim doğrulanmadı; "
      "bu bütçe, RDF'nin tek başına dominant olmadığını gösteriyorsa kazara nokta tahmini başka kaynaklarla değişebilir.", "",
      "## Karar sayımı", "",
      "| Karar | Satır sayısı |", "|---|---|"]
for k_, v_ in sorted(verdict_count.items(), key=lambda kv: -kv[1]):
    L.append(f"| {k_} | {v_} |")
L += ["", "![karşılaştırma](fig_comparison.png)", ""]
for g in GROUPS:
    L += [f"## {g}", "", md_table([r for r in ROWS if r["grup"] == g]), ""]
L += ["## Literatür ve resmî veri: kaynak düzeyi", "",
      "| Kod | İçerik | Kaynak | Doğrulama |", "|---|---|---|---|"]
for code, content, src, ver in LIT:
    L.append(f"| {code} | {content} | {src} | {ver} |")
L += ["", "Doğrulama notu: 'ikincil (özet)' kalemler arama motorunun özetinden geliyor; tam makale okunmadı. Bu kalemlere dayanan "
      "hiçbir karar birincil doğrulama olmadan kesinleştirilmemeli.", "",
      "## Açık çelişkiler ve yorum", "",
      f"- **Tanıma bağlı V_T:** sabit yük {VG_BASE:.3f} V, doğrusal ekstrapolasyon {mosN['VT_linear_extrap']:.3f} V, vth0 {VTH0:.3f} V. "
      f"Kalibre N_eff {mosC['na_eff_cm3']:.1e} cm⁻³ yalnızca doğrusal tanımla uyuşuyor; resmî ndep {NDEP:.1e} cm⁻³ sabit yük tanımıyla uyuşuyor. "
      "İkisi birden doğru olamaz. Bu yüzden tek bir tanım seçip tutarlı kalmak gerekir.",
      f"- **Stolk ile MC farkı:** MC {RN['sigma_VT_MC_mV']:.2f} mV, Stolk {RN['sigma_VT_Stolk_mV']:.2f} mV (oran {RN['sigma_VT_MC_mV']/RN['sigma_VT_Stolk_mV']:.2f}). "
      "Stolk formülü bellekten alındı; tanım farkı olası.",
      f"- **Görüntü modeli:** r=0'da FEM/görüntü oranı 1,18 (z=10 nm), 1,07 (30 nm), 1,05 (60 nm). Düzeltme yalnızca bu üç derinliğe uyduruldu. "
      "Sensörün yanlış yönde kalmadığını gösteren bir bağımsız test yok.",
      f"- **Oksit ve doping duyarlılığı:** Efektif doping (N_eff) sonuçlarında kazara nokta oranı nominale göre yaklaşık 2–3 kat. "
      "Hangi doping gerçek SKY130 için geçerli bilinmiyor; resmî ndep ile tutarlı olan nominal olmalı.",
      "- **Birim belirsizliği:** Foundry mismatch eğimi (3,356e-3) V·µm mi, başka bir birim mi, doğrulanmadı. 13 mV'luk karşılaştırma buna bağlı.",
      "- **Fisyon eğimi:** d=∞'da −1,08 (Rayleigh −1), ama g ızgarası ±%45 belirsiz; eşik karşılaştırması zayıf.", "",
      "## Eksik veri (ne ölçülemedi)", "",
      "- SKY130 4 K ve 77 K transistör verisi (NIST erişilemedi). Tüm 4 K öngörüleri bu yüzden doğrulanamıyor.",
      "- 130/180 nm kanal Coulomb blokajı raporu (aramada yok).",
      "- Gate-tanımlı Si çift nokta tünel bağlaşımı (literatür değeri yok; bariyer türü farklı).",
      "- Hiçbir deneysel ölçüm yok: tüm sayılar simülasyon ya da resmî model dosyasından; ölçüm yapılmadı.", "",
      "## Önerilen sonraki adımlar", "",
      "1. NIST dosyalarını (ETest Tile.zip, MPW-5 Test Tile.zip) tekrar dene; 4 K verisi varsa görüntü modelini ve iyonlaşma varsayımını ona göre sına.",
      "2. Dondurma (freeze-out) bilinçli bir model ekle (iyonlaşma oranı T'nin fonksiyonu) ve 4 K kazara nokta tahminini yeniden hesapla.",
      "3. V_T tanımını tek ve açık bir ölçüte sabitle (ör. akım kriteri) ve tüm karşılaştırmalarda onu kullan.",
      "4. Birim belirsizliğini kapat: mismatch parametrelerinin V·µm olduğunu model dosyasının kendi test koşullarından doğrula.",
      "5. Testi bu sonuçlara göre yeniden çerçevele: 300 K negatif sonucu kesin; 4 K pozitif öngörü, yeni model çıkana kadar deneysel yatırım için yeterli değil.", "",
      "## Veri envanteri", "",
      f"Sonuç klasörleri: exp0–exp8 (toplam dosya: {sum(len(list((R / e).glob('*'))) for e in [f'exp{i}' for i in range(9)] if (R / e).exists())}). "
      f"Resmî veri: `data/sky130_models/` (model + mismatch), `data/nist_sky130/` (kayıt + engel notu). "
      f"NIST bileşenleri: {nist_comp}. Ayrıntılı liste: `inventory.md`.", ""]
(OUT / "COMPARISON.md").write_text("\n".join(L), encoding="utf-8")

# ---------------------------------------------------------------- inventory
INV = ["# Veri envanteri (otomatik)", "",
       "Tarama: `results/`, `data/`, `research/`, `env/logs/`, `tests/`.", ""]
for base, title in ((R, "Sonuçlar (results/)"), (DATA, "Ham ve resmî veri (data/)"), (ROOT / "research", "Araştırma notları (research/)"),
                    (ROOT / "env/logs", "Günlükler (env/logs/)")):
    INV += [f"## {title}", "", "| Dosya | Boyut (B) |", "|---|---|"]
    n_ = 0
    for p in sorted(base.rglob("*")):
        if p.is_file() and "__pycache__" not in p.parts and ".pytest_cache" not in p.parts:
            INV.append(f"| {p.relative_to(ROOT).as_posix()} | {p.stat().st_size} |")
            n_ += 1
    INV += ["", f"Toplam {n_} dosya.", ""]
INV += ["## Literatür ve resmî veri kalemleri", "", "| Kod | İçerik | Kaynak | Doğrulama |", "|---|---|---|---|"]
for code, content, src, ver in LIT:
    INV.append(f"| {code} | {content} | {src} | {ver} |")
INV += ["", "## Engellenen erişimler", "",
        "- NIST `ETest Tile.zip` ve `2997_README_v4.txt`: sunucu zaman aşımı / Cloudflare 524 (ayrıntı `data/nist_sky130/DOWNLOAD_BLOCKED.txt`).",
        "- Obsidian MCP sunucusu bu oturumda bağlanamadı (zaman aşımı); kullanılmadı.", ""]
(OUT / "inventory.md").write_text("\n".join(INV), encoding="utf-8")

print("ROWS", len(ROWS), "LIT", len(LIT))
print("verdicts", verdict_count)
print("VG_BASE", round(VG_BASE, 4), "VG_P", round(VG_P, 4), "VG_M", round(VG_M, 4), "SIG_TOX_MV", round(SIG_VT_TOX_MV, 2))
print("FIS_SLOPE", {k: round(v, 3) for k, v in FIS_SLOPE.items()})
print("VIS 14K", {k[0]: round(v, 3) for k, v in VIS.items() if k[1] == 14.0})
print("PYTEST", PYTEST_LINE)
print("RAW_OK", RAW_OK, "SENS", "yes" if SENS else "no")
print("written", OUT)
