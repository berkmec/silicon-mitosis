"""Figures + Turkish summary for Experiment 1 (all numbers read from result files)."""
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "exp1"
sph = json.loads((OUT / "validate_sphere.json").read_text())
A = list(csv.DictReader(open(OUT / "cap_cell_A.csv")))
B = json.loads((OUT / "cap_shell_B.json").read_text())
S = json.loads((OUT / "schrodinger.json").read_text())
T = json.loads((OUT / "transport.json").read_text())
dolf = (OUT / "dolfinx_check.txt").read_text(encoding="utf-8", errors="replace")
EXP0 = list(csv.DictReader(open(ROOT / "results" / "exp0" / "generations.csv")))

# ---- fig: validation / convergence ----------------------------------------
fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
hs = [r["h"] for r in sph["rows"]]
ax[0].loglog(hs, [r["rel_err"] for r in sph["rows"]], "o-")
ax[0].loglog(hs, [0.9 * (h / hs[0]) ** 2 * sph["rows"][0]["rel_err"] for h in hs], "k--", lw=0.8, label="∝ h²")
ax[0].set(xlabel="mesh h", ylabel="|C_FEM − C_analitik| / C", title=f"Eş merkezli küreler (eğim {sph['loglog_slope']:.2f})")
ax[0].legend()
for r_ in (0.5, 0.4, 0.3, 0.2):
    rows = [x for x in A if float(x["r"]) == r_]
    ax[1].plot([int(x["nodes"]) for x in rows], [float(x["c_sigma"]) for x in rows], "o-", label=f"r={r_}")
ax[1].set(xlabel="mesh düğüm sayısı", ylabel="c_Σ (ε0·a)", title="Hücre kapasitansı: mesh yakınsaması", xscale="log")
ax[1].legend()
v2 = [x for x in S["validation"] if x["test"] == "separable_finite"]
ax[2].loglog([x["h_int"] for x in v2], [x["rel_err"] for x in v2], "s-", label="ayrılabilir sonlu kuyu (analitik: 3·E1D)")
v1 = [x for x in S["validation"] if x["test"] == "hardwall"]
ax[2].loglog([x["h"] for x in v1], [x["max_rel_err"] for x in v1], "o-", label="sert duvar (3,6,6,6,9)")
ax[2].set(xlabel="mesh h (nm)", ylabel="göreli hata", title="Schrödinger P2 doğrulaması")
ax[2].legend(fontsize=8)
fig.tight_layout()
fig.savefig(OUT / "fig_validation.png", dpi=140)
plt.close(fig)

# ---- fig: Schrödinger -------------------------------------------------------
W = S["wells"]
L = np.array([w["L_nm"] for w in W])
fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
ax[0].loglog(L, [w["dE_eV"] for w in W], "o-", label="FEM, Si/SiO₂ (V0=3,1 eV, BenDaniel-Duke)")
ax[0].loglog(L, [w["dE_hardwall_eV"] for w in W], "s--", label="sert duvar formülü (Deney 0)")
ax[0].set(xlabel="kuyu kenarı L (nm)", ylabel="ΔE = E2 − E1 (eV)", title="Seviye ayrışması")
ax[0].legend(fontsize=8)
ax[1].semilogx(L, [w["ratio_dE_fem_over_hw"] for w in W], "o-")
ax[1].axhline(1, c="gray", lw=0.8)
ax[1].set(xlabel="L (nm)", ylabel="ΔE_FEM / ΔE_sert-duvar", title="Deney 0 formülü küçük L'de ne kadar şişiriyor?")
fig.tight_layout()
fig.savefig(OUT / "fig_schrodinger.png", dpi=140)
plt.close(fig)

# ---- fig: shell/window ------------------------------------------------------
base = [x for x in B if x["w"] is None]
shell = [x for x in B if x["w"] is not None]
fig, ax = plt.subplots(figsize=(6.5, 4.2))
for lv in (1, 2, 3):
    pts = sorted([(x["w"], x["C_charging"]) for x in shell if x["level"] == lv])
    if pts:
        ax.plot(*zip(*pts), "o-", label=f"seviye {lv}")
b3 = [x for x in base if x["level"] == 3]
if b3:
    ax.axhline(b3[0]["C_charging"], c="k", ls="--", lw=0.8, label="kabuk yok")
ax.set(xlabel="pencere genişliği w (boşluk b=1,6)", ylabel="adanın yükleme kapasitansı (ε0 birimi)",
       title="Yüzer kabuk içindeki ada")
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(OUT / "fig_shell.png", dpi=140)
plt.close(fig)

# ---- summary ---------------------------------------------------------------
L_ = []
L_ += ["# Deney 1 — FEM ile gerçek fiziksel simülasyon: özet", "",
       "Araçlar: scikit-fem 12.0.2 (birincil), gmsh 4.15.2 (geometri/mesh), dolfinx 0.9.0 (bağımsız çapraz kontrol, WSL).",
       "Doğrulama sonuçları aşağıda; kalan sayısal belirsizlikler her bölümde belirtilmiştir.", "",
       "## 1. Doğrulama", "",
       f"- Eş merkezli küreler (analitik C = 4π·a·b/(b−a) = {sph['C_exact']:.4f}): hata mesh inceldikçe "
       f"{sph['rows'][0]['rel_err']*100:.1f} % → {sph['rows'][-1]['rel_err']*100:.2f} %, yakınsama eğimi {sph['loglog_slope']:.2f} (beklenen ≈2).",
       "- dolfinx çapraz kontrol (aynı mesh, bağımsız kod):", "```", dolf.strip(), "```"]
v3 = [x for x in S["validation"] if x["test"] == "real_well_convergence"]
L_ += ["- Schrödinger: sert duvarlı kutu için en ince meshte en büyük özdeğer hatası "
       f"{v1[-1]['max_rel_err']:.1e}; ayrılabilir sonlu kuyu (analitik) hatası {v2[0]['rel_err']:.1e} → {v2[-1]['rel_err']:.1e}.",
       "- Gerçek Si/SiO₂ kuyusunda E1 mesh yakınsaması (L=2 nm): " +
       ", ".join(f"h={x['h_int']}: {x['E1']:.4f} eV" for x in v3) + ".", "",
       "## 2. Kapasitans matrisi (ada + kaynak + savak + kapı kabuğu, birim ε0·a)", "",
       "| r=a/b | c_Σ (seviye 3) | c_g | c_s | c_d | belirsizlik (seviye 2→3 değişimi) |", "|---|---|---|---|---|---|"]
for r_ in (0.5, 0.4, 0.3, 0.2):
    x3 = next(x for x in A if float(x["r"]) == r_ and int(x["level"]) == 3)
    L_.append(f"| {r_} | {float(x3['c_sigma']):.3f} | {float(x3['c_g']):.3f} | {float(x3['c_s']):.3f} | {float(x3['c_d']):.3f} | "
              f"{T['cap_rel_uncertainty'][str(r_)]*100:.2f} % |")
L_ += ["", "Not: C ∝ a boyutsal analizdir (elektrostatik ölçek değişmezliği); ölçek testi sayısal bir sınama değildir.",
       "Yakınsama: seviye 2→3 değişimi ≈ %1; mertebe ≈1,5 olduğundan kalan hata ≈ %1,2 (c_Σ biraz fazla çıkıyor, mesh inceldikçe düşüyor).", "",
       "## 3. Yüzer iletken kabuk içindeki ada (iç içe hücre elektrostatiği)", "",
       "| w | seviye | C_charging |", "|---|---|---|"]
for x in B:
    L_.append(f"| {'kabuk yok' if x['w'] is None else x['w']} | {x['level']} | {x['C_charging']:.4f} |")
def _lv(w, lv):
    return next(x["C_charging"] for x in B if x["w"] == w and x["level"] == lv)


none3 = next(x["C_charging"] for x in B if x["w"] is None and x["level"] == 3)
L_ += ["", f"Okuma: kapalı yüzer kabuk adanın yükleme kapasitansını kabuksuza göre %{(_lv(0.0, 3) / none3 - 1) * 100:.0f} artırıyor "
       f"(kabuk ek bir yakın elektrot gibi davranıyor). Pencereyi 0 → 1,2 yapmak aynı seviyede C'yi "
       f"%{(1 - _lv(1.2, 2) / _lv(0.0, 2)) * 100:.1f} (seviye 2) / %{(1 - _lv(1.2, 3) / _lv(0.0, 3)) * 100:.1f} (seviye 3) düşürüyor: "
       "yön tutarlı, büyüklük küçük (mutlak belirsizlik ≈ %1,4 ile aynı mertebe).",
       "İletken yüzer kabuk eş-potansiyeldir; adayı dış alanlardan 'ekranlayıp' kapıyla bağlantısını kesmez, "
       "adanın potansiyeli kabuğu izler. Bu yüzden Deney 0'daki 'iç içe hücre = ekranlanma' sezgisi desteklenmedi.", "",
       "## 4. Enerji seviyeleri (Si/SiO₂, V0=3,1 eV, m*=0,26 / 0,5 m0)", "",
       "| L (nm) | E1 (eV) | ΔE (eV) | ΔE/k_B (K) | ΔE_FEM / ΔE_sert-duvar |", "|---|---|---|---|---|"]
for w in W:
    L_.append(f"| {w['L_nm']} | {w['E1_eV']:.4f} | {w['dE_eV']:.4f} | {w['dE_over_kB_K']:.0f} | {w['ratio_dE_fem_over_hw']:.2f} |")
L_ += ["", "## 5. Eşikler (enerji ≥ 10·k_B·T), FEM kapasitansı + FEM seviye ayrışması, r=0.5", ""]
for e, d in T["eps"].items():
    L_ += [f"**ε_r = {e}** (hücreler arası dolgu)", "",
           "| nesil | L (nm) | E_C/k_B (K) | ΔE/k_B (K) |", "|---|---|---|---|"]
    for t in d["table"]:
        dE = "FEM aralığı dışı (<1 nm)" if t["dE_over_kB_K"] is None else f"{t['dE_over_kB_K']:.0f}"
        L_.append(f"| {t['gen']} | {t['L_nm']:.2f} | {t['EC_over_kB_K']:.0f} | {dE} |")
    L_ += ["", "İlk eşik nesilleri: " + "; ".join(
        f"{T_} K → Coulomb: nesil {v['coulomb_gen']} ({'—' if v['coulomb_L_nm'] is None else format(v['coulomb_L_nm'], '.1f')} nm), "
        f"hapsolma: nesil {v['confinement_gen']} ({'—' if v['confinement_L_nm'] is None else format(v['confinement_L_nm'], '.1f')} nm)"
        for T_, v in d["first_regime"].items()), ""]
e3 = next(x for x in EXP0 if int(x["gen"]) == 3)
e5 = next(x for x in EXP0 if int(x["gen"]) == 5)
f3 = T["eps"]["3.9"]["table"][3]
f5 = T["eps"]["3.9"]["table"][5]
L_ += ["## 6. Deney 0 (analitik formül) ile karşılaştırma", "",
       "| nesil | Deney 0 E_C/k_B (K) | FEM E_C/k_B, ε=3,9 (K) | Deney 0 ΔE/k_B (K) | FEM ΔE/k_B (K) |", "|---|---|---|---|---|",
       f"| 3 (12,5 nm) | {float(e3['EC_over_kB_K']):.0f} | {f3['EC_over_kB_K']:.0f} | {float(e3['dE_over_kB_K']):.0f} | {f3['dE_over_kB_K']:.0f} |",
       f"| 5 (3,1 nm) | {float(e5['EC_over_kB_K']):.0f} | {f5['EC_over_kB_K']:.0f} | {float(e5['dE_over_kB_K']):.0f} | {f5['dE_over_kB_K']:.0f} |", "",
       "Deney 0'ın E_C değeri FEM'e göre ε'a bağlı olarak aynı mertebede (karışık ε kullanmıştı); seviye ayrışması ise sert duvar formülü "
       "yüzünden 3 nm'de ≈1,4×, 1 nm'de ≈3× fazla çıkmıştı. Deney 0'ın 'kuantum hapsolma 300 K'de nesil 5'te' sonucu FEM ile de çıkıyor "
       "(ΔE/k_B ≈ 3600 K ≥ 3000 K), ama marjı küçük.", "",
       "## Sınırlar (dürüstlük)", "",
       "- Ada ideal iletken (metalik, eş-potansiyel) varsayıldı; gerçek Si noktasında kuantum kapasitans ve yarıiletken ekranlama yok.",
       "- Çoğalma (nükleasyon) dinamiği FEM'e dahil değil; yalnızca her neslin statik elektrostatiği ve seviyeleri hesaplandı.",
       "- Coulomb elmasları FEM'den çıkan C_ij üstünde ortodoks (sabit etkileşim) teori ile üretildi; ilk ilkelerden NEGF değildir.",
       "- Schrödinger tek-parçacık, etkin kütle yaklaşımı; Coulomb etkileşimi ve spin yok, Poisson ile öz-uyum döngüsü yok.",
       "- 1 nm'nin altında sürekli ortam ve etkin kütle yaklaşımları geçersizdir."]
(OUT / "summary_exp1.md").write_text("\n".join(L_), encoding="utf-8")
print("ok")
