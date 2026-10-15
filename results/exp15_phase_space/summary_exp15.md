# Deney 15: Faz uzayı yapıları (elektron delikleri) boyuta göre bölünüyor mu?

Doğrulanmış elektrostatik PIC (1D1V / 2D2V / 3D3V). Sürüm 1, 2026-10-15.

Not: Çalışmayı bir alt ajan yaptı. Rapor metni ajanın teslim mesajından ana oturum tarafından kaydedildi. Doğrulama sayıları `validation.json` ile karşılaştırıldı.

- **Kod:** `experiments/exp15_phase_space/`
  - `pic.py`: CIC yükleme, FFT-Poisson, leapfrog/Boris, dispersiyon kökleri
  - `validate.py`, `run.py` (tüm koşu yapılandırmaları)
  - `structures.py`: periyodik etiketleme, birleşme ve bölünme takibi
  - `analyze.py`, `figures.py`, `realsys.py`, `bench.py`
- **Testler:** `tests/test_exp15.py` (14 test, yaklaşık 1–2 s).
- **Veri:** `validation.json`, `validation_series.npz`, `runs/<koşu>.json|npz`, `analysis.json`, `analysis_series.npz`, `realsys.json`.
- **Şekiller:** `fig_validation.png`, `fig_1d_merging.png`, `fig_2d_maps.png`, `fig_3d.png`, `fig_dimension.png`, `fig_phase_space.png`.

## 0. Kısa cevap

1. **Yapıların akıbeti boyuta bağlı, ama hiçbir boyutta kalıcı bölünme yok.**
   - **1D ve quasi-1D (dar kutulu 2D):** delikler kalıcı ve yalnız birleşiyor.
     - Gürültüden başlayan üç 1D koşusunda 7–8 delik oluşuyor, t = 600'de 1–2 delik kalıyor.
     - Birleşme sayıları 7, 5 ve 6. Bölünme iki eşikte de 0.
     - Geç zamanda W_E doygunluk değerinin %35–47'si.
   - **2D (geniş kutu) ve 3D, B = 0:** düzlemsel delikler kıvrılıyor (kink), parçalanıyor ve 20–60 ω_p⁻¹ içinde sönüyor.
     - Geç zamanda W_E/W_E,doy = 0,006–0,04.
     - Başlangıç kinetik enerjisinin %30–43'ü dik ısınmaya gidiyor.
     - Ham bölünme olayları var (2D'de 13, 3D'de 4), ama 5 ω_p⁻¹ boyunca ayrı yaşayan kalıcı bölünme 0–1.
     - Yani olan bölünme değil, parçalanıp yok olma.
   - **Temiz enine test:** başlangıç dik yönde tam simetrik, t = 25'te küçük rastgele dik itme veriliyor.
     - Oluşmuş düzlemsel delikler enine kararsız.
     - γ_⊥ ≈ 0,08 (2D, L = 40), 0,12 (2D, L = 20) ve 0,10 (3D, L = 20) ω_p. Bu, ≈ 0,15–0,23 ω_b'ye karşılık geliyor.
   - **B ∥ x, Ω_ce = 2 ω_p ≈ 3,8 ω_b:** kırılma bastırılıyor. Temiz testte t = 150'de C = 0,97 ve W_E/W_E,doy = 0,34. B = 0'da aynı değerler 0,08 ve 0,008.
2. **Eşik 1D ile 2D arasında; 2D ile 3D arasında nitel fark yok.**
   - Bu, [Morse & Nielson 1969](https://doi.org/10.1103/PhysRevLett.23.1087) ile uyumlu.
   - Faz-alanı bulgusunun ("2D'de bölünmez, 3D'de bölünür") plazmada birebir karşılığı yok.
3. **Doğrulama:**
   - Langmuir frekansı, CIC + leapfrog öngörüsüne 1,2·10⁻⁷ içinde uyuyor.
   - İki-demet en kararsız modu için γ = 0,3544; dispersiyon kökü 0,35355 (+%0,2).
   - Landau sönümü γ = −0,1538 (kuram −0,1534); ω = 1,4150 (kuram 1,4157).
4. **Silikon (n = 10¹⁷ cm⁻³, 4 K):** çarpışmasız Vlasov–Poisson geçerli değil.
   - Elektronlar donörlerde donmuş.
   - Gaz kuantum ve dejenere: λ_dB ≈ 73 nm ≫ λ_D ≈ 1,5 nm.
   - Çarpışmalar önemli: ω_pτ ≈ 0,15–1,5.

## 1. Soru

- **Faz-alanı bulgusu (Deney 11):** aktif damla 2D'de solucana dönüşüyor ama bölünmüyor; 3D'de ikiye bölünüyor.
- **Bu deneyin sistemi:** gerçek bir fiziksel sistem, çarpışmasız Vlasov–Poisson plazması. Faz uzayı 6 boyutlu (3 konum + 3 hız).
- **Yapılar:** simetrik iki-demet kararsızlığının ürettiği elektron delikleri (faz uzayı girdapları, pozitif potansiyel tümsekleri).
- **Bilinen fizik:** 1D'de girdaplar kalıcıdır ve birleşir. 2D ve 3D'de enine kararsızlıkla bozulurlar (Morse & Nielson 1969; Hutchinson [2018](https://arxiv.org/pdf/1804.08594), [2019a](https://arxiv.org/pdf/1901.06717), [2019b](https://arxiv.org/pdf/1906.04065)).
- **Katkı:** yeni fizik değil. Doğrulanmış bir yeniden üretim ve aynı ölçütlerle yapılmış bir boyut karşılaştırması: yapı sayısı, birleşme ve bölünme olayları, enine uyum, dik ısınma.

## 2. Model ve yöntem

**Birimler:**
- ω_p = 1 (n₀ = e = m = ε₀ = 1). Elektronlar için q/m = −1. İyonlar hareketsiz ve düzgün dağılmış.
- Hız birimi v₀, uzunluk birimi v₀/ω_p.

**Sayısal şema:**
- Poisson denklemi ∇²φ = n_e − 1, periyodik kutuda spektral FFT ile çözülüyor; türevde Nyquist bileşeni sıfırlanıyor.
- Yükleme ve alan toplama CIC ile (2^d köşe, `np.bincount`) ve aynı ağırlıklarla yapılıyor. Bu yüzden şema momentumu koruyor; öz-kuvvetin sıfır olduğu testte gösterildi.
- Zaman adımı leapfrog ile, dt = 0,1. B ∥ x durumunda Boris itmesi kullanılıyor (2D3V/3D3V).
- 1D, 2D ve 3D aynı kod yolundan geçiyor.

**Plazma:**
- İki Maxwell demeti: ±v₀ (x yönünde), her birinin yoğunluğu ½, her yönde v_th = 0,2 v₀.
- Sıcak dispersiyondan: γ_max = 0,335, k = 0,63, λ ≈ 10 v₀/ω_p (soğuk sınır 0,354).

**Başlangıç koşulları:**
- **(i) Gürültü:** rastgele konumlar.
- **(ii) Tohum:** (i) + düzlemsel yoğunluk bozulması δn/n = 0,05; L = 40 için m = 4, L = 20 için m = 2.
- **(iii) Sessiz-dik + itme (quietkick):**
  - x yönünde kafes konumları ve bit-ters Maxwell hızları.
  - Her dik hücreye bir kopya. CIC birim bölüşümü sayesinde yoğunluk dik koordinatlardan tam bağımsız ve dinamik bu simetriyi koruyor (W_⊥ ~ 10⁻¹⁴, testte gösterildi).
  - t = 25'te parçacıklara σ = 0,05 rastgele dik yer değiştirme veriliyor.

**Yapı tanımı:**
- φ önce 7 ω_p⁻¹ kayan zaman ortalamasından geçiriliyor (≈ 1,1 Langmuir periyodu); bu, salınımın yarattığı titreşimi bastırıyor. Ardından Gauss filtresi uygulanıyor (σ = 1 hücre).
- Yapı = φ > θ bölgesi. Komşuluk yüz üzerinden, periyodik sınırlar birleştiriliyor, en küçük boyut ≈ (1 v₀/ω_p)^d.
- Eşikler θ = 0,3 ve 0,5 φ_ref. φ_ref, aynı ailenin 1D koşusunda t = 10–40 aralığındaki maksimum: L40 için 1,44, L20 için 1,23.
- **Olay tanımları:** ardışık karelerdeki örtüşmeye bakılıyor.
  - Bir yapı ≥ 2 yapıya ayrılırsa bölünme sayılıyor.
  - ≥ 2 yapı tek yapıya dönüşürse birleşme sayılıyor.
  - Kalıcı bölünme: en az iki parçanın torunları 5 ω_p⁻¹ boyunca boş değil ve ayrık.
- **Eşik gerektirmeyen ölçütler:**
  - W_E/W_E,doy
  - Enine uyum C = W(k_⊥ = 0)/W, burada W ∝ Σk²|φ_k|²
  - Dik ısınma ΔK_⊥

## 3. Doğrulama

| Test | Ölçülen | Beklenen | Fark |
|---|---|---|---|
| (a) Langmuir, k dx = 0,098, dt = 0,1; L = 4π, 64 hücre, 64 ppc | 0,999613 | CIC + leapfrog öngörüsü 0,999613 | ω_p'den −3,9·10⁻⁴; öngörüden 1,2·10⁻⁷ |
| (a) Langmuir, k dx = 0,39 | 0,987567 | 0,987567 | ω_p'den −1,2·10⁻²; öngörüden 4·10⁻⁸ |
| (a) Langmuir, dt = 0,05 | 0,999301 | 0,999301 | 1,2·10⁻⁷ |
| (b) Soğuk iki demet, m = 4 (k v₀ = 0,6124); L = 41,04, 256 hücre, demet başına 64 ppc, tek-mod tohum 10⁻⁶ | 0,3544 (çok-mod koşusunda 0,3549) | 1 = ½[(ω−kv₀)⁻² + (ω+kv₀)⁻²] kökü: 0,35355 | +%0,2 (çok-mod +%0,4) |
| (b) Diğer modlar: m = 1, 2, 3, 5, 6 | 0,1514 / 0,2578 / 0,3320 / 0,3252 / 0,2173 | 0,1464 / 0,2602 / 0,3298 / 0,3263 / 0,2195 | +3,5 / −0,9 / +0,7 / −0,3 / −1,0 % |
| (c) Landau, 2²⁰ parçacık, sessiz başlangıç, 64 hücre, α = 0,01 | uyum (t = 6–36): γ = −0,1538 ± 0,0002, ω = 1,4150 ± 0,0002; tepe yöntemi (t = 6–25): γ = −0,1540, ω = 1,409 | Z fonksiyonu kökü: γ = −0,15336, ω = 1,41566 | +%0,3 / −%0,05 |
| (d) Enerji: Langmuir / iki demet / Landau | 5,0·10⁻³ (dt = 0,1), 1,2·10⁻³ (dt = 0,05) / 5,8·10⁻⁵ / 1,0·10⁻⁶ | — | Langmuir'deki sapma sınırlı bir salınım; dt² ile ölçekleniyor (4 kat), sürüklenme değil |
| Birim testler | 14/14 | — | — |

Notlar:
- Çok-mod koşusunda m = 1 modu ölçülemedi (0,058); sebep geçici rejim. Tek-mod koşusunda bu mod +%3,5 farkla ölçülüyor (pencere t ≥ 4/γ).
- Landau tepe yöntemi gürültü tabanına (t ≈ 36) kadar uzatılınca yaklaşık %5 sapıyor. Bu yüzden pencere t ≤ 25 ile sınırlandı.

## 4. Koşular

| Koşu | Faz uzayı | Kutu (v₀/ω_p) | Izgara | Parçacık (ppc) | Başlangıç | Ω_ce/ω_p | t_son | Enerji kayması (maks) | Süre |
|---|---|---|---|---|---|---|---|---|---|
| 1d_L100_s1/s2/s3 | 1D1V | 100 | 512 | 102 400 (200) | gürültü | 0 | 600 | 2,1 / 2,5 / 2,2·10⁻⁴ | 24–48 s |
| 1d_L40_seed | 1D1V | 40 | 128 | 51 200 (400) | tohum m = 4 | 0 | 300 | 2,6·10⁻⁴ | 13 s |
| 2d_L40_narrow (ve _ppc64) | 2D2V | 40×1,25 | 128×4 | 16 384 (32) / 32 768 (64) | tohum m = 4 | 0 | 300 | 5,0 / 5,3·10⁻⁴ | 17 / 31 s |
| 2d_L40_wide | 2D2V | 40×40 | 128² | 524 288 (32) | tohum m = 4 | 0 | 150 | 3,0·10⁻⁴ | 623 s |
| 2d_L40_wide_B2 | 2D3V | 40×40 | 128² | 524 288 (32) | tohum m = 4 | 2 | 150 | 3,5·10⁻⁴ | 258 s |
| 2d_L40_wide_quietkick | 2D2V | 40×40 | 128² | 524 288 (32) | sessiz-dik + itme | 0 | 150 | 2,8·10⁻⁴ | 203 s |
| 2d_L40_wide_quietkick_B2 | 2D3V | 40×40 | 128² | 524 288 (32) | sessiz-dik + itme | 2 | 150 | 2,7·10⁻⁴ | 196 s |
| 1d_L20_seed | 1D1V | 20 | 32 | 12 800 (400) | tohum m = 2 | 0 | 150 | 6,3·10⁻⁴ | 4 s |
| 2d_L20_seed (ve _ppc32) | 2D2V | 20×20 | 32² | 16 384 (16) / 32 768 (32) | tohum m = 2 | 0 | 150 | 4,3 / 2,3·10⁻³ | 8 / 22 s |
| 2d_L20_quietkick | 2D2V | 20×20 | 32² | 16 384 (16) | sessiz-dik + itme | 0 | 120 | 5,5·10⁻³ | 4 s |
| 3d_L20_seed | 3D3V | 20³ | 32³ | 524 288 (16) | tohum m = 2 | 0 | 120 | 2,6·10⁻³ | 649 s |
| 3d_L20_quietkick | 3D3V | 20³ | 32³ | 524 288 (16) | sessiz-dik + itme | 0 | 120 | 5,3·10⁻³ | 467 s |

- Hücre boyu: L40 ailesinde dx = 0,3125 (λ_D/dx = 0,64); L20 ailesinde ve 3D'de dx = 0,625 (λ_D/dx = 0,32).
- Koşulmayanlar: `2d_L40_wide_noseed` ve `3d_L20_seed_B2` (`run.py` içinde tanımlı; zaman yetmedi).

## 5. Sonuçlar

### 5.1 Ana tablo

Eşik θ = 0,3 φ_ref. Olaylar t_doy + 5'ten sonra sayıldı. ΔK_⊥ toplam başlangıç kinetik enerjisine oranlandı.

| Koşu | N: oluşan → maks → son | Bölünme ham / kalıcı | Birleşme | Kaybolma | W_E/W_E,doy (son) | C(100) | ΔK_⊥/K₀ |
|---|---|---|---|---|---|---|---|
| 1d_L100_s1 | 7 → 8 → 1 | 0 / 0 | 7 | 2 | 0,36 | — | — |
| 1d_L100_s2 | 7 → 7 → 2 | 0 / 0 | 5 | 0 | 0,37 | — | — |
| 1d_L100_s3 | 8 → 8 → 2 | 0 / 0 | 6 | 0 | 0,35 | — | — |
| 1d_L40_seed | 4 → 4 → 2 | 0 / 0 | 2 | 0 | 0,47 | — | — |
| 1d_L20_seed | 2 → 2 → 1 | 0 / 0 | 1 | 0 | 0,36 | — | — |
| 2d_L40_narrow | 3 → 3 → 2 | 0 / 0 | 1 | 0 | 0,46 | 0,997 | 0,007 |
| 2d_L40_wide | 5 → 8 → 0 | 13 / 0 | 1 | 19 | 0,006 | 0,16 | 0,41 |
| 2d_L40_wide_quietkick | 4 → 6 → 0 | 9 / 1 | 4 | 10 | 0,008 | 0,32 | 0,34 |
| 2d_L20_seed | 2 → 2 → 0 | 1 / 1 | 0 | 4 | 0,043 | 0,22 | 0,40 |
| 2d_L20_quietkick | 2 → 2 → 0 | 0 / 0 | 1 | 1 | 0,030 | 0,28 | 0,31 |
| 3d_L20_seed | 2 → 4 → 0 | 4 / 0 | 0 | 6 | 0,011 | 0,13 | 0,43 |
| 3d_L20_quietkick | 2 → 2 → 0 | 0 / 0 | 1 | 1 | 0,014 | 0,19 | 0,30 |
| 2d_L40_wide_B2 | 4 → 4 → 1 | 2 / 0 | 5 | 0 | 0,23 | 0,77 | 0,056 |
| 2d_L40_wide_quietkick_B2 | 4 → 4 → 2 | 0 / 0 | 2 | 0 | 0,34 | 0,98 | 0,010 |

- θ = 0,5'te de 1D ve quasi-1D koşularında bölünme yok.
- 1d_L100_s1'de sayının 7'den 8'e çıkması bir bölünme değil; eşiği geçen yeni bir delik (2 doğum).

### 5.2 1D

`fig_1d_merging.png`:
- φ(x,t) haritasında ağaç biçimli bir birleşme görülüyor.
- t ≈ 15'te 7–8 delik var, t ≈ 480'de son büyük birleşme oluyor.
- Faz uzayında t = 30'da 7 girdap, t = 500'de tek büyük girdap var.
- Dar 2D kutu (L_y = 1,25) 1D ile aynı davranıyor: C ≥ 0,997, ΔK_⊥ %0,7.

### 5.3 2D

**Gürültülü + tohum (`fig_2d_maps.png`, üst sıra):**
- t = 15'te dört düzlemsel delik var, ama şimdiden dalgalı.
- Eğik iki-demet modları gürültüden başlayıp doğrusal evrede γ ≈ 0,20 ile büyüyor; düzlemsel tohumun önü yalnız ~3,6 e-katı. Bu yüzden kırılma, deliğin enine kararsızlığıyla eğik modların karışımı.
- t = 30'da ağ biçimli parçalar var (C = 0,12). t ≥ 50'de yapılar sönük.

**Temiz test (`2d_L40_wide_quietkick`, ikinci sıra):**
- t ≤ 50'ye kadar tam düzlemsel kalıyor (C ≥ 0,995). Bu sürede 1D'deki gibi 4 → 2 birleşme oluyor.
- İtmeden sonra W_⊥, t = 30–61 arasında 124 kat büyüyor: γ_⊥ = 0,079. (En dik 10 birimlik pencere 0,12 veriyor; bu değer doğrusal olmayan hızlanmayı da içeriyor.)
- t = 75'te y boyunca 3–4 dalgalık bir kink görülüyor. t = 100'de yapı parçalanmış (W_E/W_E,doy = 0,04). t = 150'de 0,008'e iniyor.
- Sonuç: kırılma için eğik modlar gerekli değil. Bu, Hutchinson 2018 ile nitel olarak uyumlu.

### 5.4 3D

**`3d_L20_seed` (`fig_3d.png`):** t = 15'te iki düzlemsel levha var. t = 30'da yamalara bölünmüş (C = 0,04). t ≥ 50'de sönük.

**`3d_L20_quietkick`:**
- t = 50'ye kadar düzlemsel.
- W_⊥, t = 30–50 arasında 42 kat büyüyor: γ_⊥ = 0,10. Aynı kutudaki 2D değeri 0,12.
- t = 80'de levha kıvrılmış ve y–z kesitinde yamalı. t = 120'de sönmüş (max φ 0,18).

**Aynı kutuda 2D ile 3D:**

| Ölçüt | 2D | 3D |
|---|---|---|
| γ_⊥ (temiz test) | 0,12 | 0,10 |
| C(100) | 0,28 | 0,19 |
| W_E/W_E,doy(100) | 0,042 | 0,018 |
| Gürültülü koşu, W_E/W_E,doy (son) | 0,043 | 0,011 |
| Gürültülü koşu, ΔK_⊥/K₀ | 0,40 | 0,43 |

3D'de sönme biraz daha hızlı, ama nitel fark yok. Morse & Nielson da "1D ile 2D arasında büyük nitel fark, 2D ile 3D arasında küçük nicel fark" bulmuştu.

### 5.5 Manyetik alan (Ω_ce = 2 ω_p)

- **Sıçrama frekansı:** Δφ ≈ 1,44 ve k = 0,63'ten ω_b ≈ k√(Δφ/2) ≈ 0,53 ω_p (sinüzoidal kuyu, kaba tahmin). Buna göre Ω_ce/ω_b ≈ 3,8.
- **Temiz test:**
  - W_⊥, t = 30–45 arasında 13 kat büyüyor (γ ≈ 0,09), sonra toplam alan enerjisinin yaklaşık %3'ünde doyuyor.
  - t = 150'de C = 0,97 ve W_E/W_E,doy = 0,34.
  - Yapı sayısı 4 → 2: yalnız birleşme var, bölünme yok.
- **Gürültülü:**
  - C = 0,68–0,86 arasında kalıyor ve W_E/W_E,doy(150) = 0,23.
  - θ = 0,5'te 7 ham bölünme var, bunlardan 1'i kalıcı.
- **Yorum:**
  - Hutchinson'ın manyetik kararlılaşmasıyla nitel olarak uyumlu. Kuvvetli B plazmayı fiilen 1D'ye döndürüyor.
  - Eşik taranmadı.
  - Küçük ve doyan W_⊥ büyümesinin kaynağı ayrılmadı. Olası kaynaklar: itmeye manyetize tepki, eğik modlar ya da Hutchinson 2019b'deki çok yavaş yüksek-B kararsızlığı.

## 6. Faz-alanı bulgusuyla (Deney 11) karşılaştırma

| | Faz-alanı (Deney 11) | Vlasov–Poisson (bu deney) |
|---|---|---|
| Uzamış yapı | 2D'de solucan, 3D'de silindir | 1D'de nokta, 2D'de çizgi, 3D'de levha delik |
| Düşük boyut | 2D: bölünme yok | 1D: bölünme yok, yalnız birleşme |
| Yüksek boyut | 3D: ikili bölünme, kızlar yaşıyor | 2D/3D: kink → parçalanma → sönme; kalıcı bölünme 0–1 |
| Eşik | 2D ↔ 3D | 1D ↔ 2D |
| Mekanizma | yüzey gerilimi + sürüş | tuzaklı parçacık "jetting" etkisi ile elektrik "gerilimi" dengesi (Hutchinson 2018) |
| Enerji | sürekli yakıt | demet serbest enerjisi bir kez tükeniyor |
| Bastırma | — | kuvvetli B ∥ demet |

Ortak olan yalnız uzamış bir yapının boyuta bağlı enine kararsızlığı. Eşik boyutu ve sonuç farklı: burada sonuç üreme değil, yok olma. Projenin "hücre gibi bölünen yapı" arayışı açısından, çarpışmasız plazma delikleri 2D/3D'de bölünüp yaşayan bir aday değil; kuvvetli B altında ise yalnız birleşiyorlar.

## 7. Gerçek sistemler

Değerler mertebe hesabıdır; parametreler varsayımdır.

| Sistem | n | T_e | B | f_p | λ_D | Ω_ce/ω_p | ν/ω_p | ℓ/λ_D | N_D |
|---|---|---|---|---|---|---|---|---|---|
| Laboratuvar | 10¹⁶ m⁻³ | 3 eV | 0,05 T | 0,90 GHz | 0,13 mm | 1,6 | 1·10⁻⁵ | 1·10⁵ | 2·10⁴ |
| Güneş rüzgârı (1 AU) | 5 cm⁻³ | 10 eV | 5 nT | 20 kHz | 10,5 m | 0,007 | 9·10⁻¹¹ | 1·10¹⁰ | 6·10⁹ |
| Auroral bölge | 1 cm⁻³ | 100 eV | 5 µT | 9,0 kHz | 74 m | 16 | 1·10⁻¹² | 8·10¹¹ | 4·10¹¹ |

Bu üç sistemde Vlasov–Poisson çok iyi geçerli. Delik boyutu (~25–50 λ_D) laboratuvarda mm, uzayda yüzlerce metre mertebesinde. Ω_ce/ω_p = 0,007 olan durumda bu deneye göre delikler 2D/3D'de kırılmalı; 16 olan durumda (ve laboratuvar örneğinde 1,6) kararlı kalmaları beklenir.

**Silikon** (n = 10¹⁷ cm⁻³, m* = 0,26 m₀, ε = 11,7, 4 K):

- **Plazma ölçekleri:**
  - ω_p = 1,02·10¹³ s⁻¹ (f_p = 1,63 THz).
  - ħω_p = 6,7 meV ≫ kT = 0,34 meV.
  - Klasik λ_D = 1,5 nm.
- **Dejenerelik:**
  - Fermi enerjisi: tek vadi için E_F = 3,0 meV (T_F = 35 K); 6 vadi için 0,73 meV (T_F = 8,4 K). Yani 4 K'de gaz dejenere.
  - λ_TF = 3,6 / 1,8 nm.
- **Kuantum ve bağlaşım:**
  - Termal λ_dB = 73 nm, n λ_dB³ = 39.
  - N_D = 3·10⁻⁴, r_s = 5,6, klasik Γ = 27.
- **Donma:** Si:P için Mott sınırı ≈ 3,7·10¹⁸ cm⁻³ (literatür değeri; kaynağı burada doğrulanmadı). 10¹⁷'de 4 K'de taşıyıcılar donörlere donar; serbest elektron ancak enjeksiyon, aydınlatma ya da kapı birikimiyle olur.
- **Çarpışmalar:**
  - μ = 100–1000 cm²/Vs varsayımıyla τ = 1,5·10⁻¹⁴–1,5·10⁻¹³ s, ω_pτ = 0,15–1,5 ve serbest yol ℓ ≈ 0,9–9,5 nm.
  - ℓ'nin delik ölçeğine (~75 nm) ulaşması için μ ≳ 8000 cm²/Vs gerekir.
- **Sonuç:** Vlasov–Poisson burada uygulanamaz. Dört sebep var:
  1. donma;
  2. kuantum ve dejenerelik (Wigner/Hartree tanımı gerekir);
  3. çarpışmalar;
  4. N_D ≪ 1 (güçlü bağlaşım).

## 8. Sınırlar

1. **CPU:** dört çekirdek başka oturumların işleriyle doluydu.
   - ppc düşürüldü: 2D 64 → 32, 3D 32 → 16. Süre t ≤ 150 ile sınırlandı.
   - Izgara yakınsaması yapılmadı.
   - ppc duyarlılığına yalnız quasi-1D (32/64) ve L20-2D (16/32) koşularında bakıldı; ikisinde de nitel sonuç aynı.
2. **3D kutu küçük:** 2 dalga boyu, λ_D/dx = 0,32, 16 ppc. Enerji kayması 2,6–5,3·10⁻³ (2D L40'ta ~3·10⁻⁴).
3. **Gürültülü koşular:** kırılmaya eğik modlar karışıyor (γ ≈ 0,15–0,20). Temiz itme testi bu iki etkiyi ayırıyor.
4. **Yapı sayıları:** eşiğe, zaman ortalamasına ve en küçük boyuta duyarlı. "Kalıcı" ölçütündeki 5 ω_p⁻¹ bir seçim. Eşik gerektirmeyen ölçütler aynı sonucu veriyor.
5. **γ_⊥:** değerler tek koşulardan; kutuya ve pencere seçimine bağlı olarak 0,08–0,12 arasında.
   - Hutchinson'ın kuramıyla nicel karşılaştırma yapılmadı.
   - Buradaki delikler onun Maxwell-denge deliklerinden farklı: iki-demet ürünü ve derin (ψ ≈ 0,7 m v₀²/e ≫ T_demet = 0,04 m v₀²).
6. **B testi:** tek bir B değeri; eşik taraması ve 3D+B koşusu yok.
7. **Gerçek sistem sayıları:** mertebe hesabı; örnek plazma parametreleri ve silikon mobilitesi varsayım.

## 9. Kaynaklar

- I. H. Hutchinson, [Transverse instability of electron phase-space holes in multi-dimensional Maxwellian plasmas](https://arxiv.org/pdf/1804.08594), arXiv:1804.08594 (2018).
- I. H. Hutchinson, [Transverse instability magnetic field thresholds of electron phase-space holes](https://arxiv.org/pdf/1901.06717), Phys. Rev. E 99, 053209 (2019).
- I. H. Hutchinson, [Electron phase-space hole transverse instability at high magnetic field](https://arxiv.org/pdf/1906.04065), J. Plasma Phys. (2019), [doi:10.1017/S0022377819000564](https://doi.org/10.1017/S0022377819000564).
- R. L. Morse, C. W. Nielson, [One-, Two-, and Three-Dimensional Numerical Simulation of Two-Beam Plasmas](https://doi.org/10.1103/PhysRevLett.23.1087), Phys. Rev. Lett. 23, 1087 (1969).
- M. Oppenheim, D. L. Newman, M. V. Goldman, [Evolution of Electron Phase-Space Holes in a 2D Magnetized Plasma](https://doi.org/10.1103/PhysRevLett.83.2344), Phys. Rev. Lett. 83, 2344 (1999).
