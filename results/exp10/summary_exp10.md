# Deney 10 — φ'nin fiziksel karşılığı: SKY130 inversiyon tabakası

> **Düzeltme (2026-10-11 incelemesi, `results/exp11/review_exp6_exp10_exp11.md`):**
> 1. Aşağıdaki "Yeniden kontrol 3" kapı haritası (`check.py`, `with_corr=True`) kapı ekranlı **olmayan** korelasyon vekiline dayanıyor.
>    Kapı Hartree'yi ve değişimi ekranlıyor, korelasyonu ekranlamıyor. Bu model içi bir tutarsızlık. Tutarlı (kapı ekranlı) HF'de tabaka
>    C = 0,05–3000 µF/cm² ve n_s = 10⁹–10¹³ cm⁻² aralığının **her yerinde kararlı** (S/K_kin ≥ +1,006). Bu yüzden "yakın kapı bölünmeyi
>    açıyor / süperkapasitör sezgisi doğrulandı" cümlesi **geri alındı**. Ana tablodaki kararsız varyantların (A, B, D, E) hepsinde
>    değişim-korelasyon Hartree'den az ekranlanmış; tek tutarlı varyant C ve o kararlı.
> 2. Haritadaki "yok" hücreleri tarama tabanının (10¹⁰) ürünü. Vekille n_max şöyle: SiO₂ 4,148 nm için 3,4·10⁹; 2 nm için 9,8·10⁹;
>    ε = 20, 4,148 nm için 7,2·10⁹ cm⁻². Kaba ızgara n_max'ı biraz düşük veriyor (ör. 2,0 → 2,4·10¹⁰).
> 3. "Hartree |K|'nın 11–27 katı": E varyantında Fang–Howard çarpanı unutulmuş, doğrusu **11–21**.
> 4. ε_c = (11,7 + ε_kapı)/2 seçimi burada belgelenmemişti ve ε = 20 sonuçlarını belirliyor. Tabakanın kendi kalınlığı seri yaklaşık
>    4,2 µF/cm² ekliyor; etkin kapı bağlaşımı bunu geçemez.
> 5. Kapının iki zıt rolü var. Düzgün tabakanın ayrışmasında dengeleyici Hartree'yi ekranlar; var olan hücrenin bölünmesinde iticiliği
>    ekranlar. Bu rapor yalnız birincisini inceler.

Yorum: φ ≈ (n − n₀)/n₀, inversiyon tabakasındaki elektron yoğunluğunun normalize sapması.
Exp8 modelindeki iki terim bu yorumda şuna karşılık gelir:
- çift kuyunun eğriliği f''(φ̄) ↔ n₀²·K(n₀), K = d²(E/A)/dn² (kinetik + değişim + korelasyon, yerel ters kompresibilite)
- Coulomb terimi g(q) ↔ n₀²·V_H(q), V_H(q) = e² / (ε₀ q (ε_Si + ε_ox coth(q t_ox)))  (Si | SiO₂ | kapı, tabaka arayüzde)

Parametreler: m* = 0,19, ε_Si = 11,7, ε_ox = 3,9, t_ox = 4,148 nm (resmî SKY130 model dosyası), ε̄ = 7,8.

## Doğrulama (tests/test_egas2d.py, 6 test)
- Değişim integrali, ekransız çekirdekte analitik −(4/3π)e²k_F/(4πε₀ε) ile 1e-8 içinde.
- Hartree atom birimlerinde makaledeki −4√2/(3π r_s) formuyla 1e-8 içinde.
- Hartree–Fock kritik r_s = π/√2 (g = 2) sayısal olarak yeniden üretildi.
- Kapı çekirdeği iki limitte doğru (t_ox → ∞: 2D ekransız; q → 0: e²/C_ox).
- Korelasyon: Attaccalite vd., PRL 88, 256601 (2002), Tablo II (ζ = 0), PDF metninden okundu
  (`data/literature/attaccalite2002.txt`). Toplam enerji makalenin DMC Tablo I değerleriyle r_s = 2, 5, 10, 20'de %1,5 içinde.

## Sonuç
| Varyant | K < 0 eşiği n_c (cm⁻²) | r_s,c | Hartree dahil kararsızlık (q ≤ k_F) | n_s = 10¹²'de K/K_kin | n_s = 10¹²'de V_H(k_F)/|K| |
|---|---|---|---|---|---|
| A tek vadi, ε̄, xc (DMC) — referans | 1,63e12 | 2,03 | yalnız n ≤ 1,8e10 | −0,31 | 11 |
| B Si (g=4), HF, ε̄ | 2,74e12 | 1,57 | yalnız n ≤ 3,2e9 | −0,65 | 14 |
| C Si, HF, kapı ekranlı değişim | 2,57e12 | 1,62 | yok | −0,53 | 17 |
| D Si, HF kapı + korelasyon (vekil) | 3,10e12 | 1,47 | yalnız n ≤ 3,2e9 | −0,81 | 11 |
| E D + Fang–Howard kalınlık | 1,52e12 | 2,11 | yalnız n ≤ 3,2e9 | −0,33 | 27 |

Okuma:
1. **φ'nin çift kuyusu fiziksel olarak var.** SKY130 eşik civarında (n_s ≈ 10¹² cm⁻², r_s ≈ 2,6) yerel ters kompresibilite
   bütün varyantlarda negatif. Yani φ = (n − n₀)/n₀ yorumu, exp8'deki f''(φ̄) < 0 koşulunu gerçekten sağlıyor.
2. **Ama bölünme yok.** Kapı ekranlı Hartree terimi k_F'ye kadar her dalga boyunda |K|'nın 11–27 katı. Düzgün tabaka kararlı.
   Exp8 diliyle: Coulomb gücü kritik değerin çok üzerinde, desen (hücre) oluşmuyor.
3. Kararsızlık yalnızca n ≲ 3·10⁹–2·10¹⁰ cm⁻² (r_s ≳ 20–50) için çıkıyor. Bu, SKY130'un çalışma aralığının dışında ve
   sıvı (homojen gaz) tanımının kendisinin çöktüğü Wigner kristali bölgesine yakın; bu model orada güvenilir değil.

## Sınırlar
- Korelasyon parametrizasyonu tek vadi içindir; Si (iki vadi) için yalnızca vekil (D, E). Değişim (B, C) iki vadi için tam.
- Korelasyonda kapı ekranlaması yok (yalnız değişimde var). Çift kapılı ekranlı xc için arXiv 2410.07445 ilgili, okunmadı.
- Yerel (LDA tipi) enerji q ≳ k_F'de geçersiz; kararlılık yalnız q ≤ k_F için değerlendirildi. Tam cevap için yerel alan
  düzeltmeli statik dielektrik fonksiyon gerekir.
- Arayüz pürüzü, katkı potansiyeli (exp6) ve vadi yarılması dahil değil.
- Fang–Howard form faktörü hafızadan, yalnız duyarlılık için.

## Yeniden kontrol (kullanıcı itirazı üzerine, `experiments/exp10_phi_physics/check.py`)
1. **Hesap hatası yok.** V_H(k_F)/K_kin kodda 10,3883, bağımsız analitik Thomas–Fermi oranı g/(a_B*·k_F) = 10,3883. Coulomb'un
   baskın olması bir hata değil, standart fizik: q_TF ≈ 10·k_F.
2. **q ≤ k_F sınırı sonucu üretmiyor.** SKY130 yığınında q ≤ 2k_F ve 4k_F'de de bölünme yok.
3. **Bölünme kapı çok yaklaşınca çıkıyor** (q ≤ k_F, tarama 1e10–1e13 cm⁻²):

| Kapı dielektriği | t_ox = 4,148 nm (SKY130) | 2 nm | 1 nm | 0,5 nm |
|---|---|---|---|---|
| SiO₂ (ε = 3,9) | yok | yok | n ≤ 2·10¹⁰ | n ≤ 5·10¹⁰ |
| yüksek-κ (ε = 20) | yok | n ≤ 1,3·10¹⁰ | n ≤ 2,5·10¹⁰ | n ≤ 3,2·10¹⁰ |

   Bu bölge r_s ≈ 12–26 (güçlü korelasyonlu, Wigner kristaline yakın). Exp8'in "kapı yaklaşınca bölünme başlar" bulgusuyla
   nitel olarak uyumlu: "süperkapasitör" sezgisi (yakın, güçlü kapı) gerçekten bölünmeyi açıyor.
4. **Uyarı:** En kararsız dalga sayısı her durumda tam q = k_F'de, yani benim geçerlilik kesme sınırımda. Bu yüzden çıkan hücre
   boyutu (λ ≈ 160–350 nm) güvenilir değil; doğru boyut için q'ya bağlı yanıt (Lindhard + yerel alan düzeltmesi) gerekir.
