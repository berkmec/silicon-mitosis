# The quantum attempt — currently unsuccessful

Research by **counfield** — <https://counfield.com/>

This folder documents the original aim of the project and why, at present, it has **not succeeded**. The code and data are in
`experiments/exp0 … exp10` and `results/exp0 … exp10`; the detailed reports and the roadmap are in Turkish
(`research/roadmap/ROADMAP.md`, `results/exp9/COMPARISON.md`).

## 1. Objective

The starting idea was a *self-replicating transistor–supercapacitor cell* that, under electrical drive, divides into ever smaller
copies until it reaches a quantum regime. The concrete goal was:
- **Coulomb blockade or qubit behaviour** on the **standard 130 nm SkyWater SKY130 process**;
- **at room temperature**, without cryogenics.

## 2. Approach

We first built validated numerical tools and then confronted each requirement with the process limits:

| stage | experiment | what was computed | validation |
|---|---|---|---|
| 0 | Exp0 | analytic toy model of a recursive cell (replaced later by real tools) | — |
| 1 | Exp1 | 3D FEM capacitance of island + leads + gate shell; Schrödinger levels in finite wells | sphere capacitor 0.33 %; scikit-fem vs dolfinx 8×10⁻¹⁶ |
| 1 | Exp2 | chains of islands (12 / 6 / 3 nm): charge steps vs temperature | single-island limit |
| 1 | Exp3 | double quantum dot at 300 K: polarisation, acoustic-phonon T1, gate fidelity | Qiskit Aer simulation agrees with the channel calculation |
| 1 | Exp4 | tunnel coupling (FEM vs. WKB), Landau–Zener with Lindblad damping, SKY130 scale gap | tunnel splitting FEM/WKB 1.03–1.11; Landau–Zener vs. analytic 1.2×10⁻³ |
| 2a | Exp5 | 1D Poisson–Schrödinger of the SKY130 MOS stack (official model parameters) | Airy triangular well < 0.5 % |
| 2b | Exp6 | random-dopant potential in a 126 × 420 nm channel: accidental quantum dots | single charge + images (analytic) |
| 2c | Exp7 | statistical test plan for an accidental-dot measurement (≥ 3 positives need 51–313 devices) | — |
| 2c | Exp8 | screened-Coulomb Cahn–Hilliard (bifurcation) and Coulomb-fission models; a prepared domain did not divide | phase-field interface tension 0.41 % |
| — | Exp9 | comparison of all data (64 rows, 17 literature items) | — |
| — | Exp10 | 2D electron-gas compressibility of the SKY130 inversion layer | Hartree–Fock analytic |

## 3. Principal findings

**3.1 Scale gap.** Coulomb blockade needs $E_C \ge 10\,k_BT$. From the FEM capacitance of the cell geometry, the largest admissible island edge is:

| T | ε = 3.9 | ε = 11.7 |
|---|---|---|
| 300 K | 0.39 nm | 0.13 nm |
| 77 K | 1.5 nm | 0.5 nm |
| 4 K | 29 nm | 9.8 nm |
| 1 K | 117 nm | 39 nm |

The smallest SKY130 gate is 150 nm, so room-temperature operation is short by a factor of **about 380–1150**. A 150 nm island satisfies the criterion only below ≈ 0.8 K, which requires a dilution refrigerator.

**3.2 No coupling between gate-defined cells.**
- For a 20 meV barrier in silicon, the WKB tunnel amplitude falls to about $10^{-29}$ across the 210 nm poly-to-poly spacing.
- Cells defined by SKY130 gates are therefore uncoupled.

**3.3 No coherence at room temperature.**
- Even in the most favourable geometry (L = 1 nm, d = 0.5 nm) at 300 K:
  - initial polarisation ≤ 0.74;
  - gate fidelity ≤ 0.78, even optimistically;
  - fidelity 0.5 for a 10 ps electrical gate.
- Only the single-phonon acoustic channel was included, so the true T2 is lower.

**3.4 The accidental-dot prediction was retracted.** A random-dopant Monte Carlo first suggested that 4–21 % of minimum devices would show Coulomb blockade at 4 K. A later review found three problems:
- **Edge artefact:** excluding the four box edges gives 0/80.
- **Over-correction:** the FEM correction exceeded the physical limit of 1.5.
- **Missing physics:** screening and the floating body at 4 K were absent.

The prediction is withdrawn (`results/exp11/review_exp6_exp10_exp11.md`).

**3.5 Inversion-layer compressibility (Exp10).**
- The inversion layer has negative compressibility at low density ($n_s \lesssim 1.5$–$3\times10^{12}$ cm⁻²).
- The gate-screened Hartree term is 11–21 times larger, so the layer does not spontaneously split into dots.
- An earlier claim that a nearby gate would enable splitting was an artefact of an unscreened correlation proxy and is retracted.

**3.6 No cryogenic process data.**
- The NIST SKY130 cryogenic dataset could not be downloaded: the server returned errors (`data/nist_sky130/`).
- Coulomb blockade has been reported in 28 nm (14 K) and 40 nm (50 mK) commercial devices; we read these reports from their
  abstracts. Our search found no report for the 130/180 nm nodes.

## 4. Verdict

> **Currently unsuccessful.** A room-temperature quantum device, or a qubit, on standard 130 nm SKY130 CMOS is incompatible with the physics computed here: the scale, coupling and coherence arguments are independent of one another.

The verdict could change only through a different, much harder experiment:
- fabricate minimum-size SKY130 devices (e.g. a multi-project wafer);
- measure them at ≤ 1 K, looking for single-electron (Coulomb-blockade) behaviour;
- treat this as a search for single-electron transistor behaviour, not for qubits.

## 5. What came out of the attempt

The tools and the negative results led to the classical line of the project:
- dissipative, self-dividing **current filaments** in silicon (Exp12–13, Exp17);
- dimension-dependent division in active matter (Exp11, Exp14).

See the main [README](../README.md).

## 6. Figures from this stage

- FEM charge-stability diagrams: [`results/exp1/fig_diamonds_fem.png`](../results/exp1/fig_diamonds_fem.png)
- Schrödinger levels: [`results/exp1/fig_schrodinger.png`](../results/exp1/fig_schrodinger.png)
- Chain of islands vs temperature: [`results/exp2/fig_chain.png`](../results/exp2/fig_chain.png)
- Double dot / fidelity at 300 K: [`results/exp3/fig_dqd.png`](../results/exp3/fig_dqd.png), [`results/exp3/fig_fidelity.png`](../results/exp3/fig_fidelity.png)
- MOS 1D Poisson–Schrödinger: [`results/exp5/fig_mos1d.png`](../results/exp5/fig_mos1d.png)
- Random-dopant potential (prediction retracted): [`results/exp6/fig_rdf.png`](../results/exp6/fig_rdf.png)
- Comparison of all data: [`results/exp9/fig_comparison.png`](../results/exp9/fig_comparison.png)

---

## Türkçe özet

**Hedef:** standart 130 nm SKY130 sürecinde, oda sıcaklığında kuantum davranışı (Coulomb blokajı / qubit).

**Sonuç: şimdilik başarısız.** Üç bağımsız sebep var:
1. **Ölçek:** oda sıcaklığı için 0,13–0,39 nm'lik ada gerekiyor, en küçük kapı ise 150 nm; fark yaklaşık 380–1150 kat.
2. **Bağlaşım:** 210 nm'lik boşlukta tünelleşme ~10⁻²⁹ çarpanına düşüyor; hücreler birbirine bağlanamıyor.
3. **Koherans:** 300 K'de polarizasyon en fazla 0,74, kapı doğruluğu en fazla 0,78.

Diğer bulgular:
- Kazara kuantum noktası tahmini (%4–21) bir kenar artefaktı çıktığı için **geri çekildi**.
- Kriyojenik SKY130 verisi (NIST) indirilemedi.
- Tek gerçekçi yol: en küçük boyutlu SKY130 cihazlarını ürettirip ≤ 1 K'de ölçmek ve tek elektron davranışı aramak.

Bu çalışmanın araçları ve olumsuz sonuçları, projenin klasik koluna yol açtı: silikonda kendiliğinden bölünen akım filamanları (ana README).
