# Deney 1 — FEM ile gerçek fiziksel simülasyon: özet

Araçlar: scikit-fem 12.0.2 (birincil), gmsh 4.15.2 (geometri/mesh), dolfinx 0.9.0 (bağımsız çapraz kontrol, WSL).
Doğrulama sonuçları aşağıda; kalan sayısal belirsizlikler her bölümde belirtilmiştir.

## 1. Doğrulama

- Eş merkezli küreler (analitik C = 4π·a·b/(b−a) = 12.5664): hata mesh inceldikçe 12.0 % → 0.33 %, yakınsama eğimi 1.91 (beklenen ≈2).
- dolfinx çapraz kontrol (aynı mesh, bağımsız kod):
```
dolfinx 0.9.0: nodes=1597 dofs_in=206 dofs_out=790
C_dolfinx = 12.8519479305
C_skfem   = 12.8519479305
|dolfinx-skfem|/skfem = 8.293e-16
C_analytic = 12.5663706144  (dolfinx hatası 2.273e-02)
```
- Schrödinger: sert duvarlı kutu için en ince meshte en büyük özdeğer hatası 1.7e-03; ayrılabilir sonlu kuyu (analitik) hatası 6.0e-03 → 7.4e-04.
- Gerçek Si/SiO₂ kuyusunda E1 mesh yakınsaması (L=2 nm): h=0.35: 0.6395 eV, h=0.25: 0.6371 eV, h=0.18: 0.6354 eV, h=0.13: 0.6350 eV.

## 2. Kapasitans matrisi (ada + kaynak + savak + kapı kabuğu, birim ε0·a)

| r=a/b | c_Σ (seviye 3) | c_g | c_s | c_d | belirsizlik (seviye 2→3 değişimi) |
|---|---|---|---|---|---|
| 0.5 | 22.989 | 17.549 | 2.719 | 2.721 | 0.97 % |
| 0.4 | 19.847 | 13.834 | 3.006 | 3.006 | 1.05 % |
| 0.3 | 17.743 | 11.324 | 3.211 | 3.209 | 1.12 % |
| 0.2 | 16.267 | 9.544 | 3.361 | 3.362 | 1.22 % |

Not: C ∝ a boyutsal analizdir (elektrostatik ölçek değişmezliği); ölçek testi sayısal bir sınama değildir.
Yakınsama: seviye 2→3 değişimi ≈ %1; mertebe ≈1,5 olduğundan kalan hata ≈ %1,2 (c_Σ biraz fazla çıkıyor, mesh inceldikçe düşüyor).

## 3. Yüzer iletken kabuk içindeki ada (iç içe hücre elektrostatiği)

| w | seviye | C_charging |
|---|---|---|
| kabuk yok | 1 | 5.0943 |
| kabuk yok | 2 | 4.8747 |
| kabuk yok | 3 | 4.7673 |
| 0.0 | 1 | 5.5783 |
| 0.0 | 2 | 5.4239 |
| 0.0 | 3 | 5.3493 |
| 0.2 | 2 | 5.4235 |
| 0.4 | 2 | 5.4233 |
| 0.8 | 2 | 5.4019 |
| 1.2 | 1 | 5.5347 |
| 1.2 | 2 | 5.3697 |
| 1.2 | 3 | 5.2928 |

Okuma: kapalı yüzer kabuk adanın yükleme kapasitansını kabuksuza göre %12 artırıyor (kabuk ek bir yakın elektrot gibi davranıyor). Pencereyi 0 → 1,2 yapmak aynı seviyede C'yi %1.0 (seviye 2) / %1.1 (seviye 3) düşürüyor: yön tutarlı, büyüklük küçük (mutlak belirsizlik ≈ %1,4 ile aynı mertebe).
İletken yüzer kabuk eş-potansiyeldir; adayı dış alanlardan 'ekranlayıp' kapıyla bağlantısını kesmez, adanın potansiyeli kabuğu izler. Bu yüzden Deney 0'daki 'iç içe hücre = ekranlanma' sezgisi desteklenmedi.

## 4. Enerji seviyeleri (Si/SiO₂, V0=3,1 eV, m*=0,26 / 0,5 m0)

| L (nm) | E1 (eV) | ΔE (eV) | ΔE/k_B (K) | ΔE_FEM / ΔE_sert-duvar |
|---|---|---|---|---|
| 1.0 | 1.6356 | 1.4543 | 16876 | 0.34 |
| 1.5 | 0.9700 | 0.9524 | 11052 | 0.49 |
| 2.0 | 0.6371 | 0.6346 | 7364 | 0.58 |
| 3.0 | 0.3340 | 0.3340 | 3876 | 0.69 |
| 4.0 | 0.2049 | 0.2050 | 2379 | 0.76 |
| 6.0 | 0.0997 | 0.0998 | 1158 | 0.83 |
| 8.0 | 0.0587 | 0.0588 | 682 | 0.87 |

## 5. Eşikler (enerji ≥ 10·k_B·T), FEM kapasitansı + FEM seviye ayrışması, r=0.5

**ε_r = 3.9** (hücreler arası dolgu)

| nesil | L (nm) | E_C/k_B (K) | ΔE/k_B (K) |
|---|---|---|---|
| 0 | 100.00 | 12 | 4 |
| 1 | 50.00 | 23 | 17 |
| 2 | 25.00 | 47 | 70 |
| 3 | 12.50 | 94 | 279 |
| 4 | 6.25 | 187 | 1074 |
| 5 | 3.12 | 375 | 3617 |
| 6 | 1.56 | 749 | 10433 |
| 7 | 0.78 | 1499 | FEM aralığı dışı (<1 nm) |
| 8 | 0.39 | 2998 | FEM aralığı dışı (<1 nm) |

İlk eşik nesilleri: 300 K → Coulomb: nesil None (— nm), hapsolma: nesil 5 (3.1 nm); 77 K → Coulomb: nesil 7 (0.8 nm), hapsolma: nesil 4 (6.2 nm); 4 K → Coulomb: nesil 2 (25.0 nm), hapsolma: nesil 2 (25.0 nm)

**ε_r = 11.7** (hücreler arası dolgu)

| nesil | L (nm) | E_C/k_B (K) | ΔE/k_B (K) |
|---|---|---|---|
| 0 | 100.00 | 4 | 4 |
| 1 | 50.00 | 8 | 17 |
| 2 | 25.00 | 16 | 70 |
| 3 | 12.50 | 31 | 279 |
| 4 | 6.25 | 62 | 1074 |
| 5 | 3.12 | 125 | 3617 |
| 6 | 1.56 | 250 | 10433 |
| 7 | 0.78 | 500 | FEM aralığı dışı (<1 nm) |
| 8 | 0.39 | 999 | FEM aralığı dışı (<1 nm) |

İlk eşik nesilleri: 300 K → Coulomb: nesil None (— nm), hapsolma: nesil 5 (3.1 nm); 77 K → Coulomb: nesil 8 (0.4 nm), hapsolma: nesil 4 (6.2 nm); 4 K → Coulomb: nesil 4 (6.2 nm), hapsolma: nesil 2 (25.0 nm)

## 6. Deney 0 (analitik formül) ile karşılaştırma

| nesil | Deney 0 E_C/k_B (K) | FEM E_C/k_B, ε=3,9 (K) | Deney 0 ΔE/k_B (K) | FEM ΔE/k_B (K) |
|---|---|---|---|---|
| 3 (12,5 nm) | 75 | 94 | 322 | 279 |
| 5 (3,1 nm) | 299 | 375 | 5156 | 3617 |

Deney 0'ın E_C değeri FEM'e göre ε'a bağlı olarak aynı mertebede (karışık ε kullanmıştı); seviye ayrışması ise sert duvar formülü yüzünden 3 nm'de ≈1,4×, 1 nm'de ≈3× fazla çıkmıştı. Deney 0'ın 'kuantum hapsolma 300 K'de nesil 5'te' sonucu FEM ile de çıkıyor (ΔE/k_B ≈ 3600 K ≥ 3000 K), ama marjı küçük.

## Sınırlar (dürüstlük)

- Ada ideal iletken (metalik, eş-potansiyel) varsayıldı; gerçek Si noktasında kuantum kapasitans ve yarıiletken ekranlama yok.
- Çoğalma (nükleasyon) dinamiği FEM'e dahil değil; yalnızca her neslin statik elektrostatiği ve seviyeleri hesaplandı.
- Coulomb elmasları FEM'den çıkan C_ij üstünde ortodoks (sabit etkileşim) teori ile üretildi; ilk ilkelerden NEGF değildir.
- Schrödinger tek-parçacık, etkin kütle yaklaşımı; Coulomb etkileşimi ve spin yok, Poisson ile öz-uyum döngüsü yok.
- 1 nm'nin altında sürekli ortam ve etkin kütle yaklaşımları geçersizdir.