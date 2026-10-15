# Deney 11: Hangisi daha güçlü bölüyor? Süperkapasitör fikri (A), disipatif yapı (B) ve kazara nokta (C)

Sürüm 2 (2026-10-11, incelemeden sonra). İnceleme raporu: `review_exp6_exp10_exp11.md`.
Kod: `experiments/exp11_compare_division/` (scan, compare, phys, model_a_phys, b3d, gs_ramp, gs_track, extend,
checks_review, conv, pilot), `src/qcell/active.py`, `src/qcell/grayscott.py`, testler `tests/test_active.py`.
Şekiller: `fig_exp11_sim.png`, `fig_exp11_cells.png`, `fig_exp11_phys.png`.
Veri: `runs/*.json|npz` (A, A2, A3, B, GS ve uzatmalar), `classification.json`, `stats.json`, `comparison.csv`,
`phys.json`, `model_a_phys.json`, `b3d*.json`, `gs_ramp.json`, `gs_track.json`, `checks_review.json`.

**Sürüm 1'den farklar.**
- "B açık ara güçlü" hükmü geri alındı.
- 3B bölünmenin disipasyona özgü olmadığı gösterildi.
- 2B sınıflandırma şekil ölçütüyle düzeltildi.
- A için şarj şiddeti tarandı ve dambıl koşusu uzatıldı.
- Gray–Scott sayımı 5/42 oldu.
- Fizik tahmini E_g'ye ve 2B tabakaya göre düzeltildi.
- Coulomb teriminin işareti düzeltildi.
- C (Aşama 2b kazara noktası) karşılaştırmaya eklendi.

## 0. Kısa cevap

1. **Tek bir kazanan yok.** Bölünmeyi üç şey belirliyor: boyut, uzun menzilli itme ve otokatalitik aktivatör–inhibitör yapısı.
   "Disipatif mi süperkapasitör mü" sorusu bunların yerine geçmiyor.
   - **Aynı 2B çekirdek:** büyüme→ikili bölünme A'da 0/28, B'de (aktif damla) 0/30. Kalıcı bölünme 0/28'e karşı 2/30 ve bu fark
     istatistiksel olarak anlamsız (Fisher p = 0,49). İkisi de uzama, solucan ve labirent üretiyor.
   - **3B:** aktif damla Δν = −0,002 … +0,004 aralığının hepsinde iki eşit kıza bölünüyor ve bölünme yarıçapı Δν boyunca sürekli
     değişiyor. Δν = 0 bunun içinde. Δν = 0 ise sızıntılı şarjla, yani nötrleştirici arka plandaki ekransız 3B Coulomb'lu yüklü
     damlayla aynı denklem. Yani 3B'de "süperkapasitör" de bölünüyor; nükleer fisyonun sıvı damla resmine benzer bir durum.
   - **Gray–Scott (otokatalitik):** 2B'de de çoğalıyor. Taranan (F, k) kutusunun 5/42'sinde var, %5 gürültüyle 6/6 dayanıklı ve
     yeni lekeler bölünmeyle doğuyor (doğum 0).
2. **Bölünme kaçınılmaz değil.** Her durumda eşik ve dar bir bant var.
3. **MOSFET kanalı 2B ve kapılı.** Orada süperkapasitör ve aktif damla hücre benzeri bölünme vermiyor; yalnız Gray–Scott türü
   otokatalitik sistemler aday. Silikonda düzgün elektron tabakası tutarlı modelde (kapı ekranlı HF) hiçbir kapıda ayrışmıyor.
4. **Kuantum hedefi.** 4 K'de E_C ≥ 10 kT için hücre yarıçapı geometriye bağlı olarak R ≈ 6–20 nm olmalı. Disipatif hücrenin silikon
   karşılığında olay başına enerji ≈ E_g (kink mekanizması) alınınca en iyimser durum ≈ 123 nm (2B tabaka), ≈ 1,6 µm (3B), tipik
   ≈ 0,4 µm. Kuantuma ulaşmıyor.
5. **C (kazara nokta, Aşama 2b) hücre üretir ama bölünmez.** Önceki "%4–21, a ≈ 15–23 nm" tahmini kenar artefaktıydı ve geri çekildi.
   Yalnız rastgele katkı potansiyelinde çukurlar sık (%39–95) ama büyük (a ≈ 35–41 nm, E_C/kT ≈ 3,5–4). Ekranlama yok; sayı verilemez.

## 1. Fikirden denkleme

Ortak çekirdek: Cahn–Hilliard faz ayrışması, Zwicker vd. 2017 eki Denk. (S.1)–(S.5) birimlerinde. Konsantrasyon Δc, uzunluk w,
zaman t₀ = w²/D birimindedir.

∂ₜc = ∇²[ c(c−1)(2c−1) − ¼∇²c + G ∗ (c − ⟨c⟩) ] + s(c)   (Coulomb terimi itici; sürüm 1'de işaret yanlış yazılmıştı, kod doğruydu)

| | A: süperkapasitör | B1: aktif damla | B2: Gray–Scott |
|---|---|---|---|
| Hücre | yüklü taşıyıcı domeni | yakıtla sürülen damla | leke (otokatalitik) |
| Uzun menzilli terim | G(q) = A(1 − e^{−2qd})/q (kapı d uzaklıkta) | etkin k/q² (Δν = 0'da tam) | — |
| Sürüş | şarj s = J veya s = J − λc | üretim ν₊ − kc / yıkım −ν₋ − k(c−1) | besleme F, söndürme k |

Parametreler: ε = ν₊/k, Δν = ν₋ + ν₊ − k, k = 0,01. Δν = 0'da s(c) = ν₊ − kc tek bir doğrudur. Bu, λ = k çekirdekli bir
Ohta–Kawasaki gradyan akışı (Glotzer vd. 1995) ve sızıntılı şarjla aynı denklem (test edildi).

## 2. Doğrulama

| Kontrol | Sonuç |
|---|---|
| Birim testler | 33/33 |
| Izgara | dx 0,5 → 0,33: ΔR ≤ 0,07 |
| Zaman adımı | gecikme hatası dt·S ile kabaca doğrusal (dt·S = 0,1 / 0,2 / 0,4 / 1,0); durağan yarıçap aynı (9,424) |
| float32 / float64 | \|Δc\| = 1,7·10⁻⁵ (t = 500), yarıçap aynı (`checks_review.json`) |
| Kutu | L = 64 → 128: ΔR = 0,30 (t = 4000) ve zamanla büyüyor. Geç zaman yakınsaması gösterilmedi |
| 3B, makale Fig. 3A hızları | t = 1750 tek (R = 11,64, en-boy 2,2) → t = 2000 iki eşit kız (R = 9,92 / 9,93; kızların en-boy oranı doğumda 2,6, sonra 5,1) |
| Gray–Scott, Pearson noktası (256²) | 1 → en çok 187 → 172 leke (doğum-ölüm dengesi) |

## 3. Sonuçlar

### 3.1 B (disipatif)

| Alt model | Taranan | Sonuç |
|---|---|---|
| Aktif damla 2B (ε 0,08–0,20 × Δν −0,006…+0,004, L = 128 w) | 30 | 6 çözündü, 12 yuvarlak durağan, 9 solucan/labirent (R_eq/r_iç 2,9–8), 1 solucan kopması (ε = 0,20, Δν = −0,002, t ≈ 19 000; solucan kutudan uzun), 2 büyük başlangıç damlası ayrıldı (1'i kalıcı). **Büyüme→ikili bölünme 0/30.** Δν = 0 çizgisinde 0/5, ama ε = 0,20 noktası T'de hâlâ değişiyordu |
| Aktif damla 2B, makale hızları (L = 64 / 128 / 192) | 3 | hepsi solucan; solucan boyu kutuyla ölçekleniyor |
| Aktif damla 3B (ε = 0,20; Δν −0,004 … +0,004; 64³, L = 48 w) | 6 | Δν −0,002 … +0,004 (Δν = 0 dahil): t ≈ 2000'de iki eşit kız; R_kız 11,1 → 6,1 (Δν arttıkça küçülüyor). Δν = −0,004: t = 4000'de tek, R = 19,8 (kutu sınırlı). **Düzeltme (Exp14, 2026-10-15):** bu yapı kompakt bir damla değil, içi boş bir sıvı kabuk. Orta düzlemde merkez konsantrasyonu t = 1000'de 0,78 iken t = 2000'de 0,03'e düşüyor ve kapalı bir seyreltik bölge oluşuyor. Eski sınıflandırma yalnız kovaryans en-boy oranına baktığı için bunu görmedi |
| Gray–Scott (F 0,010–0,058 × k 0,045–0,067) | 42 | **5 çoğalan leke**, 6 labirent, 16 düzgün, 14 söndü, 1 az leke |
| Gray–Scott dayanıklılık | 6 | 6/6 çoğaldı, N_max 73–84 |
| Gray–Scott ölçek değişimi (D/2 ≡ kutu √2 büyük) ve soy takibi | 4 evre + takip | N 83 → 142 → 293 → 587. Takipte doğum 0: yeni lekeler bölünmeyle doğuyor. Ölçek değişiminde N'nin artması kaçınılmaz, gösterilen şey artışın bölünmeyle gerçekleşmesi |

2B ile 3B farkı (yorum): 3B'de ince bir silindir Rayleigh–Plateau kararsızlığıyla kürelere kopar, 2B'de bir şerit kopunca çevresi azalmaz.
(Exp14'te boyut sürekli değişken yapılarak sınandı: bölünme yalnız bir boyut penceresinde oluyor, 2 < d ≲ 3,5–4. Ayrıntı:
`results/exp14_dimension/summary_exp14.md`.)
Aynı çözücü 3B'de ikiye bölüyor, 2B'de solucan üretiyor. 2B'deki tek kopma kutu boyuna bağlı olabilir; L ≥ 256 ile sınanmadı.

### 3.2 A (süperkapasitör)

| Alt grup | Taranan | Sonuç |
|---|---|---|
| Güçlü Coulomb (A = 0,02 / 0,08; d = ∞ / 3 / 1) + A = 0 kontrolü | 15 | 0 kalıcı bölünme; hücre 6/15 koşuda yaşadı (Coulomb yüklü hücreyi dağıtıyor) |
| Zayıf Coulomb (A = 0,001 / 0,003 / 0,01) | 13 | 0 kalıcı bölünme; hücre hep yaşadı; şarjla R ≈ 46'ya kadar yuvarlak büyüdü. Dengede A = 0,01 dambıl oldu; t = 60 000'e uzatıldı, bölünmedi, solucana döndü |
| Şarj şiddeti (J = 10⁻⁴ / 5·10⁻⁴ / 2·10⁻³) | 12 | 0 bölünme; 12/12 ortamdan doğum (uzak alan spinodali aştı, 13–41 yeni hücre). Coulomb'suz kontrolde de aynı. Yani çoğalma şarjın aşırı doygunluğundan, Coulomb'dan değil |

Kapı iki zıt rol oynuyor:
- Düzgün tabakayı domenlere ayırmayı kolaylaştırıyor, çünkü dengeleyici Hartree terimini ekranlıyor.
- Var olan hücrenin ikiye bölünmesini engelliyor, çünkü bölünmeyi iten itmeyi ekranlıyor.

Fiziksel karşılık (Exp10 çözücüsü, `model_a_phys.json`): bu hesap (i) sürecine, yani düzgün tabakanın ayrışmasına aittir.

| Model | Sonuç |
|---|---|
| Tutarlı, kapı ekranlı HF | C = 0,05–3000 µF/cm², n_s = 10⁹–10¹³ cm⁻²: **hiçbir yerde kararsızlık yok** (S/K_kin ≥ +1,006) |
| HF + kapı ekranlı olmayan korelasyon vekili (tutarsız) | kararsızlık çıkıyor (n_s ≥ 10¹⁰'da C* ≥ 1,76 µF/cm², SKY130 C_ox'ta n_s ≤ 3,4·10⁹). Bu, tutarsızlığın ürünü |

### 3.3 C (kazara nokta, Aşama 2b; bölünme değil, düzensizlikten doğan hücre)

Düzeltilmiş Monte Carlo (`results/exp6/review_rerun.json`, aynı 80 örnek):

| Kural | Nominal katkı | Etkin katkı |
|---|---|---|
| Orijinal (yalnız x kenarları dışlanıyor, kubbe dahil) | 3/80 | 3/80 |
| Dört kenar dışlanıyor, kubbe dahil | 0/80 | 0/80 |
| Dört kenar dışlanıyor, yalnız rastgele RDF | 31/80 (a medyan 41 nm) | 76/80 (35 nm) |

Ekranlama, tam çekirdek ve gerçek S/D/STI geometrisi olmadan bu sayılar güvenilir değil.

## 4. Karşılaştırma

| Ölçüt | A: süperkapasitör | B1: aktif damla | B2: Gray–Scott | C: kazara nokta |
|---|---|---|---|---|
| 2B'de tek hücre kalıcı bölündü mü? | 0/28 (+12 ortamdan doğum) | 2/30 (ikili 0) | evet, 5/42 | bölünme yok (hücre düzensizlikten) |
| 3B'de bölünme | ekransız Coulomb ≡ Δν = 0: **evet** | evet (Δν'nin geniş aralığında) | — | — |
| Gürültüye dayanıklılık | — | — | 6/6 | — |
| Hücre boyu (silikon karşılığı) | silikonda ayrışma yok | 2B'de filaman | en iyimser 123 nm, tipik ≈ 0,4 µm | a ≈ 35–41 nm (ekransız, güvenilmez) |
| 4 K kuantum ölçütü (R ≈ 6–20 nm) | — | — | ulaşmıyor | sınırda / belirsiz (E_C/kT ≈ 3,5–17) |

## 5. Silikona çeviri (mertebe tahmini, `phys.py`)

Hücre boyu R ≈ ρℓ, burada ℓ = √(Dτ) ve ρ = 1,16 (3B ikili bölünmeden). Hücre başına en az N_min = 100 taşıyıcı ve gereken güç bütçesi
iki ayrı alt sınır koyuyor.

| Geometri ve olay enerjisi | En iyimser R (D = 0,1 cm²/s, 10 W/cm²) | Tipik R (D = 1, 1 W/cm²) | N_min duyarlılığı |
|---|---|---|---|
| 2B tabaka, E = E_g (kink mekanizması) | 123 nm | 390 nm | N_min = 10 → 69 nm |
| 2B tabaka, E = 45 meV (sığ katkı; kink literatürü desteklemiyor) | 55 nm | 175 nm | N_min = 10 → 31 nm |
| 3B, E = E_g | ≈ 1,6 µm | — | — |

4 K kuantum hedefi (E_C ≥ 10 kT):
- Si içindeki ada (FEM sabiti, hacim eşdeğeri) R ≈ 6–9 nm;
- izole Si küresi R ≤ 17,9 nm;
- SKY130 kapılı disk (yalnız plaka) R ≤ 29,8 nm (kenar alanı eklenince daha küçük).

Disipatif hücrenin sınırı tünel bariyeri değil, taşıyıcılar sıcak. Boy tutsa bile Coulomb blokajı kendiliğinden gelmez.

## 6. "Bölünme kaçınılmaz" itirazına cevap

Haklı olduğun kısım: Exp8 ve Exp10'daki "bölünme yok" sonuçlarım dengeye aitti. Sürülen sistemler (Gray–Scott; 3B aktif damla)
bölünüyor. Ama inceleme iki şey daha gösterdi:
- 3B'de denge eşdeğeri bir uzun menzilli itme (ekransız Coulomb) de bölüyor. Bölünme yalnız disipasyona özgü değil.
- 2B kapılı tabakada (MOSFET) ne süperkapasitör ne aktif damla hücre benzeri bölünüyor.

Gray–Scott'ta bile çoğalma taranan uzayın yalnız 5/42'sinde var. Doğru ifade: "bölünme belirli bir sürüş, boyut ve etkileşim
menzili bandında mümkün; kaçınılmaz değil."

## 7. Geçerlilik sınırları

- 2B koşularda periyodik kutu L = 128 w ve T ≤ 2·10⁴ t₀; koşular durağanlıkta erken durduruldu.
  - Solucan kopması kutu boyunun ürünü olabilir.
  - Başlangıç yarıçapı R₀ = 8 olduğundan "çözündü" sonucu kritik yarıçapa bağlı olabilir.
  - Zwicker taraması tek tohumla yapıldı.
- 3B: 64³, dx = 0,75 w (makalede 0,5), periyodik kutu, tek tohum. Kutu ikinci kuşak için ve Δν = −0,004 için küçük.
- A modelinde 2B tabaka çekirdeği, nötrleştirici arka plan ve düzgün yük enjeksiyonu var; kuantum ve Fermi istatistiği yok.
- Fizik tahmini mertebe düzeyinde; D, P_max, f ve N_min varsayım. Joule ısınması dahil değil.
- C: ekranlama yok, çekirdek görüntü modelinde (tam çekirdek yaklaşık 1,3 kat büyük σ verir), gövde topraklı varsayıldı (4 K'de değil).

## 8. Sonraki adım

Exp12: 4 K'de geniş bir SKY130 nFET'te yüzer (donmuş) gövde ve darbe iyonlaşması (E ≈ E_g) ile kink bistabilitesi ve genişlik boyunca
akım filamanı. Kodu `experiments/exp12_filament/`. Bu, silikonda gerçekçi tek disipatif desen adayı ve klasik.

## Kaynaklar

- Zwicker vd., Nat. Phys. 13, 408 (2017), [arXiv 1603.01571](https://arxiv.org/abs/1603.01571)
- Pearson, Science 261, 189 (1993), [arXiv patt-sol/9304003](https://arxiv.org/abs/patt-sol/9304003)
- Glotzer, Di Marzio, Muthukumar, PRL 74, 2034 (1995), [APS](https://link.aps.org/doi/10.1103/PhysRevLett.74.2034); [cond-mat/0606011](https://arxiv.org/pdf/cond-mat/0606011)
- 4 K CMOS kink: [2104.08467](https://arxiv.org/pdf/2104.08467), [1811.09987](https://arxiv.org/pdf/1811.09987), [2003.00207](https://arxiv.org/pdf/2003.00207), [1809.09013](https://arxiv.org/pdf/1809.09013)
- Si e-h damlaları: [Caltech tezi](https://thesis.caltech.edu/17754/), [Smith & Wolfe](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.51.7521), [0806.2050](https://arxiv.org/pdf/0806.2050), [Yoffa & Adler](https://research.ibm.com/publications/electron-hole-drops-in-silicon)
- Akım filamanları: [Hüpper vd. 1993](https://link.aps.org/doi/10.1103/PhysRevB.47.15515), [Schöll & Drasdo 1990](https://link.springer.com/article/10.1007/BF01309347)
- EDL kapı: [Fujimoto & Awaga 2013](https://doi.org/10.1039/c3cp50755f)
