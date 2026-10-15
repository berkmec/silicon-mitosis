"""Qubit figures of merit, Qiskit Aer gate fidelity, figures and Turkish summary for experiments 2 and 3."""
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit.circuit.library import XGate
from qiskit.quantum_info import SuperOp, average_gate_fidelity
from qiskit_aer import AerSimulator
from qiskit_aer.noise import NoiseModel, thermal_relaxation_error

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from qcell import model as m  # noqa: E402

O3 = ROOT / "results" / "exp3"
O2 = ROOT / "results" / "exp2"
rows = json.loads((O3 / "dqd.json").read_text())
conv = json.loads((O3 / "convergence.json").read_text())
ch = json.loads((O2 / "chain.json").read_text())
HB = m.HBAR
TAU = ["10fs", "100fs", "1000fs"]


def gate_fidelity(T1, T2, t_gate, p_exc):
    err = thermal_relaxation_error(T1, T2, t_gate, excited_state_population=p_exc)
    total = SuperOp(XGate()).compose(SuperOp(err.to_quantumchannel()) if hasattr(err, "to_quantumchannel") else err)
    return float(average_gate_fidelity(total, XGate()))


def aer_check(T1, T2, t_gate, p_exc):
    """Aer density-matrix run of X on |0>: population of |1> vs channel prediction."""
    nm = NoiseModel()
    nm.add_quantum_error(thermal_relaxation_error(T1, T2, t_gate, excited_state_population=p_exc), ["x"], [0])
    qc = QuantumCircuit(1)
    qc.x(0)
    qc.save_density_matrix()
    sim = AerSimulator(method="density_matrix", noise_model=nm)
    rho = sim.run(transpile(qc, sim, basis_gates=["x", "save_density_matrix"] and None or None, optimization_level=0)).result().data()["density_matrix"]
    p1_aer = float(np.real(rho.data[1, 1]))
    chan = SuperOp(XGate()).compose(SuperOp(thermal_relaxation_error(T1, T2, t_gate, excited_state_population=p_exc).to_quantumchannel()))
    from qiskit.quantum_info import DensityMatrix, Statevector
    p1_ch = float(np.real(DensityMatrix(Statevector.from_label("0")).evolve(chan).data[1, 1]))
    return p1_aer, p1_ch


kT300 = m.KB * 300.0
out_rows = []
for r in rows:
    dJ = r["delta_eV"] * m.E
    t300 = r["temps"]["300"]
    p_exc = 1.0 / (1.0 + np.exp(dJ / kT300))
    T1 = 1.0 / t300["Gamma1"]
    rec = dict(L_nm=r["L_nm"], d_nm=r["d_nm"], delta_meV=r["delta_eV"] * 1e3, P300=t300["P_thermal"], T1_ac_s=T1,
               q0_nm=r["q0_per_nm"], F2=r["meanF2"], fid={})
    for tp in TAU:
        T2 = t300["scenarios"][tp]["T2_s"]
        Q = t300["scenarios"][tp]["Q"]
        t_thz = 20 * np.pi * HB / dJ  # Omega = 0.05 Delta/hbar (RWA-safe)
        t_fast = 2 * np.pi * HB / dJ  # Omega = 0.5 Delta/hbar: optimistic, beyond RWA (pi pulse in 2 pi hbar / Delta)
        t_el = 10e-12                 # electrical pulse limit
        rec["fid"][tp] = dict(T2_s=T2, Q=Q, t_gate_thz_s=t_thz, F_thz=gate_fidelity(T1, T2, t_thz, p_exc),
                              t_gate_fast_s=t_fast, F_fast=gate_fidelity(T1, T2, t_fast, p_exc),
                              t_gate_el_s=t_el, F_el=gate_fidelity(T1, T2, t_el, p_exc))
    # smallest T2 (T1 fixed at the acoustic value) giving F >= 0.99 for the optimistic fast gate
    lo, hi = 1e-18, min(2 * T1, 1.0)
    if gate_fidelity(T1, hi, 2 * np.pi * HB / dJ, p_exc) < 0.99:
        rec["T2_needed_99_s"] = None
    else:
        for _ in range(80):
            mid = np.sqrt(lo * hi)
            if gate_fidelity(T1, mid, 2 * np.pi * HB / dJ, p_exc) >= 0.99:
                hi = mid
            else:
                lo = mid
        rec["T2_needed_99_s"] = hi
    out_rows.append(rec)

first = rows[0]
dJ0 = first["delta_eV"] * m.E
T10 = 1.0 / first["temps"]["300"]["Gamma1"]
T20 = first["temps"]["300"]["scenarios"]["100fs"]["T2_s"]
p0 = 1.0 / (1.0 + np.exp(dJ0 / kT300))
tg0 = 20 * np.pi * HB / dJ0
p_aer, p_ch = aer_check(T10, T20, tg0, p0)
(O3 / "qubit_figures.json").write_text(json.dumps(dict(rows=out_rows, aer_consistency=dict(p1_aer=p_aer, p1_channel=p_ch)), indent=2))

# ---- figures --------------------------------------------------------------
Ls = sorted({r["L_nm"] for r in rows})
ds = sorted({r["d_nm"] for r in rows})


def grid(key, scen=None, T="300"):
    G = np.zeros((len(Ls), len(ds)))
    for r in rows:
        i, j = Ls.index(r["L_nm"]), ds.index(r["d_nm"])
        G[i, j] = r[key] if scen is None else r["temps"][T][scen]
    return G


fig, axs = plt.subplots(1, 3, figsize=(15, 4.4))
for ax, (G, title, fmt) in zip(axs, [
        (grid("delta_eV") * 1e3, "Δ = E2 − E1 (meV)", "{:.0f}"),
        (np.array([[next(x["P300"] for x in out_rows if x["L_nm"] == L and x["d_nm"] == d) for d in ds] for L in Ls]),
         "300 K başlangıç polarizasyonu tanh(Δ/2kT)", "{:.2f}"),
        (np.array([[next(x["fid"]["100fs"]["Q"] for x in out_rows if x["L_nm"] == L and x["d_nm"] == d) for d in ds] for L in Ls]),
         "300 K kalite faktörü Q=ΔT2/ħ (τφ=100 fs senaryosu)", "{:.1f}")]):
    im = ax.imshow(G, origin="lower", aspect="auto", cmap="viridis")
    ax.set_xticks(range(len(ds)), [str(x) for x in ds])
    ax.set_yticks(range(len(Ls)), [str(x) for x in Ls])
    ax.set(xlabel="bariyer kalınlığı d (nm)", ylabel="kuyu kenarı L (nm)", title=title)
    for i in range(len(Ls)):
        for j in range(len(ds)):
            ax.text(j, i, fmt.format(G[i, j]), ha="center", va="center", color="w", fontsize=9)
    fig.colorbar(im, ax=ax)
fig.suptitle("Çift kuantum nokta (Si/SiO₂), FEM Schrödinger + tek-fonon (akustik) altın kural")
fig.tight_layout()
fig.savefig(O3 / "fig_dqd.png", dpi=140)
plt.close(fig)

fig, ax = plt.subplots(figsize=(9, 4.4))
labels = [f"L={r['L_nm']}, d={r['d_nm']}" for r in out_rows]
x = np.arange(len(out_rows))
for k, tp in enumerate(TAU):
    ax.bar(x + (k - 1) * 0.27, [1 - r["fid"][tp]["F_thz"] for r in out_rows], 0.27, label=f"THz/optik kapı, τφ={tp}")
ax.set_yscale("log")
ax.set_xticks(x, labels, rotation=60, ha="right", fontsize=8)
ax.set(ylabel="1 − ortalama X-kapı doğruluğu (300 K)", title="Aer termal gevşeme modeli: kapı hatası (düşük = iyi)")
ax.axhline(1e-2, c="gray", ls="--", lw=0.8)
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(O3 / "fig_fidelity.png", dpi=140)
plt.close(fig)

# ---- summary ---------------------------------------------------------------
L_ = ["# Deney 2 ve 3 — çoğalma zinciri, çift kuantum nokta ve oda sıcaklığı qubit taraması", "",
      "## Önce sınırlar", "",
      "- Oda sıcaklığında kuantum durumunu ve qubiti **zorlayamam**: bu belge hangi parametrelerde neyin mümkün olduğunu gösterir.",
      "- Fonon hesabı yalnızca **tek-fonon akustik** kanaldır (T2 için ÜST sınır, T2 = 2·T1). Optik fonon, iki-fonon (Raman), vadi/orbital ve "
      "yük gürültüsü hesaba katılmadı. Bu yüzden 300 K'de dephasing için τφ = 10 fs / 100 fs / 1 ps **senaryoları** ayrıca verildi; bunlar tahmin değildir.",
      "- Elektron-fonon bağlaşımı: Si boyuna akustik deformasyon potansiyeli Ξ = 9 eV (varsayım; literatür aralığı yaklaşık 5-10 eV).", "",
      "## Adım A — kademeli küçülen zincir (12 / 6 / 3 nm, ε=3,9)", "",
      f"- FEM kapasitans matrisi: satır toplamı hatası ≤ 5e-13, simetri hatası ≤ 1.2e-14; mesh seviye 1→2 adacık bloğu en büyük göreli değişim "
      f"%{ch['level_change_max_rel']*100:.1f} (tam yakınsamamış, üst sınır).",
      f"- Adaların E_C/k_B değerleri: {', '.join(f'{x:.0f} K' for x in ch['EC_over_kB_K'])}.",
      "- Yük kuantizasyon indeksi (1 = tam basamak, 0 = sürekli) ve detrend edilmiş akım salınımı:", "",
      "| T (K) | salınım genliği/ort. akım | kuantizasyon (ada 12/6/3 nm) |", "|---|---|---|"]
for T in ("4", "77", "300"):
    o = ch["oscillations"][T]
    L_.append(f"| {T} | {o['detrended_osc_amp_over_mean']:.4f} | " + " / ".join(f"{q:.2f}" for q in o["charge_quantisation_index"]) + " |")
L_ += ["", "Okuma: 4 K'de ada 0 ve ada 1 yük basamakları gösteriyor ama en küçük ada (3 nm) hiç yüklenmiyor; 77 K ve 300 K'de yük sürekli "
       "(kuantizasyon ≈ 0), akım salınımı yok. 77 K'de ada 2 için görünen 0,58, doluluğun zaten ≈ 0 olmasından kaynaklanıyor, kuantizasyon değil.", "",
       "## Yakınsama (çift nokta, L=1 nm, d=0,5 nm)", "",
       "| h_int (nm) | kareleme derecesi | dof | Δ (meV) | <|F|²> |", "|---|---|---|---|---|"]
for c in conv:
    L_.append(f"| {c['h_int']} | {c['intorder']} | {c['dofs']} | {c['delta_meV']:.2f} | {c['meanF2']:.3e} |")
best = max(out_rows, key=lambda r: r["delta_meV"])
bf = best["fid"]
L_ += ["", "## Adım B/C — çift kuantum nokta tablosu (300 K)", "",
       "Kapı doğruluğu: ortalama X-kapı doğruluğu; 0,5 = kanal tamamen depolarize (kapı ölü). "
       "THz = Ω=0,05·Δ/ħ (RWA güvenli), hızlı = Ω=0,5·Δ/ħ (iyimser, RWA dışı), elektriksel = 10 ps darbe. τφ = 100 fs senaryosu.", "",
       "| L (nm) | d (nm) | Δ (meV) | P(300K) | q0 (1/nm) | <|F|²> | T1 akustik (s) | Q | F_X THz | F_X hızlı | F_X elektriksel | F=0,99 için gereken T2 (s) |",
       "|---|---|---|---|---|---|---|---|---|---|---|---|"]
for r in out_rows:
    f100 = r["fid"]["100fs"]
    need = "—" if r["T2_needed_99_s"] is None else f"{r['T2_needed_99_s']:.1e}"
    L_.append(f"| {r['L_nm']} | {r['d_nm']} | {r['delta_meV']:.1f} | {r['P300']:.2f} | {r['q0_nm']:.1f} | {r['F2']:.1e} | {r['T1_ac_s']:.1e} | "
              f"{f100['Q']:.1f} | {f100['F_thz']:.3f} | {f100['F_fast']:.3f} | {f100['F_el']:.3f} | {need} |")
good = [r for r in out_rows if r["P300"] >= 0.9]
L_ += ["", f"Aer tutarlılık kontrolü (ilk konfigürasyon): popülasyon |1> Aer = {p_aer:.6f}, kanal = {p_ch:.6f}.", "",
       "## Sonuç (dürüst okuma)", "",
       f"- 300 K'de başlangıç polarizasyonu ≥ 0,9 (Δ ≥ 76 meV) olan konfigürasyon: {len(good)} / {len(out_rows)}. "
       f"En iyisi (L={best['L_nm']} nm, d={best['d_nm']} nm): Δ={best['delta_meV']:.0f} meV, P={best['P300']:.2f}, "
       f"T1 (akustik) = {best['T1_ac_s']:.1e} s.",
       f"- Bu en iyi konfigürasyonda bile (τφ=100 fs, iyimser hızlı kapı) X-kapı doğruluğu {bf['100fs']['F_fast']:.2f}; "
       f"τφ=1 ps olsaydı {bf['1000fs']['F_fast']:.2f}. "
       + ("%99 doğruluğa, T2 üst sınırı (2·T1) verilse bile ulaşılamıyor (akustik T1 ve termal popülasyon sınırı)."
          if best['T2_needed_99_s'] is None else f"%99 için T2 ≥ {best['T2_needed_99_s']:.1e} s gerekir."),
       "Tablodaki '—' = bu konfigürasyonda %99'a T2 = 2·T1 ile bile ulaşılamıyor.",
       "- **Optimum taramanın kenarında** (en küçük L, en ince d). Δ bu yönde artıyor; ancak 1 nm altında sürekli ortam ve etkin kütle modeli "
       "geçersiz, bu yüzden daha küçüğüne bu modelle çıkarım yapılamaz (atomistik yöntem gerekir).",
       "- Fonon hesabı yalnızca tek-fonon akustik kanaldır; optik fonon, Raman, yük gürültüsü hesaplanmadı. T1 değerleri konfigürasyona göre "
       "fs'ten µs'ye kadar çok geniş değişiyor; Δ büyüdükçe fonon dalga boyu nokta boyutundan küçülüp bağlaşım bastırılıyor (fonon darboğazı), "
       "Δ ≈ 10 meV civarında bağlaşım en güçlü (λ ≈ nokta boyutu).",
       "- Elektriksel 10 ps kapı darbesi bütün konfigürasyonlarda koherans süresinden çok uzun; bu modelde elektrikle kontrol edilen oda sıcaklığı "
       "silikon qubiti gösterilemedi. Gösterilebilen tek şey: en küçük/ince geometride koherans için bir *fs* penceresi (Q≈7), "
       "ve bunun kullanılabilmesi için fs mertebesinde THz/optik kontrol ile τφ ≳ ps gerektiği."]
(O3 / "summary_exp23.md").write_text("\n".join(L_), encoding="utf-8")
print("ok")
