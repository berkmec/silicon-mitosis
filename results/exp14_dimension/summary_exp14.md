# Deney 14: Bölünme boyuta bağlı mı? d = 2…6 boyutta aktif Cahn–Hilliard ve 6B hiperküp

Tarih: 2026-10-15.

- **Kod:**
  - `src/qcell/active_sharp.py`: doğrusal keskin arayüz kuramı, gerçel d.
  - `src/qcell/active_sharp_nl.py`: kimyasal potansiyelde doğrusal olmayan kuram.
  - `src/qcell/active_axi.py`: O(d−1) simetrik (z, r) çözücü, gerçel d.
  - `experiments/exp14_dimension/`: `axi_run`, `cart_run` (yalın d boyutlu Kartezyen çözücü), `hyper6d`, `passive_rp`, `theory`, `window_theory`, `hyper_scan`, `figures`.
- **Testler:** `tests/test_exp14.py`.
- **Şekil:** `fig_exp14.png`.
- **Veri:** `axi/*.json|npz`, `outcomes.json`, `theory.json`, `window_theory.json`, `passive_rp.json`, `hyper6d.json`, `convergence_h.json`, `hyper_scan.json`.

Not: `experiments/exp14_6d_cell`, `results/exp14` ve `src/qcell/cell6d.py` başka bir oturuma ait ve bu deneyin parçası değil. Orada "6B" 3B ızgarada 6 alan anlamına geliyor, bölünme elle yazılmış bir kural ve 5 testin 5'i başarısız. Dokunulmadı.

> **Düzeltme (2026-10-15, yayın öncesi denetim).**
> - **"Doğrusal etkin damla kuramı %27 fazla" ifadesi kesin değildir.** 9 damlada faz alanı / doğrusal kuram oranı 0,70–0,77'dir
>   (doğrulama tablosu).
> - **3B ve 6B'nin ayrık sürüş pencerelerinde bölündüğü sonucu `axi/Z_nmk*_d3.json` ve `axi/Z_nmk*_d6.json` koşularına dayanır**
>   (Δν = (nmk − 0,8) · 0,01).
>   - Bu koşulardan beşi önceki `outcomes.json`'da ve şekilde yoktu: Z_nmk0.5_d3, Z_nmk2_d3, Z_nmk3_d3, Z_nmk1.2_d6, Z_nmk1.6_d6.
>   - `figures.py` düzeltildi; Δν artık formülden hesaplanıyor.
>   - `outcomes.json` ve `fig_exp14.png` yeniden üretildi; önceki 24 satır değişmedi.
> - **6B hiperküp doğrulaması yalnız büyüme evresini kapsar** (t ≤ 40).
> - **"Not" satırındaki başka oturuma ait dosyalar** yayımlanan depoya dahil edilmedi.

## 0. Kısa cevap

1. **Evet, bölünme boyuta bağlı. Ama tekdüze değil: bölünme yalnız bir boyut penceresinde oluyor.**
   - Alt sınırın altında (d ≈ 2–2,2) damla şerit olup uzuyor.
   - Üst sınırın üstünde damla bölünemeden içten boşalıp **sıvı kabuğa** dönüşüyor.
   - Pencere iç tüketim hızı Δν arttıkça yukarı kayıyor ve genişliyor:

   | Δν | bölünme penceresi (d) |
   |---|---|
   | −0,004 | 2,25 – 2,5 (3B'de kabuk) |
   | −0,001 (Zwicker 2017) | 2,3 – 3,5 (3,75'ten itibaren kabuk) |
   | +0,004 | 2,5 – 4 (5'te kabuk) |
   | +0,012 ve +0,022 | **6B'de ikiye bölünme** |

2. **6B için:**
   - Zwicker kimyasında 6B damla bölünmüyor, t = 500'de içi boşalıp kabuk oluyor.
   - Güçlü iç tüketimde (Δν ≥ +0,012) 6B damla uzayıp dambıl olup iki kompakt kıza bölünüyor.
   - Yani yüksek boyutta bölünme mümkün, ama daha güçlü sürüş istiyor.
3. **Gerçek 6B hiperküpte Cahn–Hilliard koşuldu:** 16⁶ = 16,8 milyon nokta, aktif damla büyümesi.
   - O(5) simetrik indirgeme bu koşuyu t ≤ 40 boyunca her noktada ≤ %2,5 farkla (h = 0,5 ile ~%1) yeniden üretiyor.
   - 6B bölünme ve kabuk sonuçları bu doğrulanmış indirgemeden geliyor.
   - Tam 6B ızgarada bölünme koşmak bu makinede mümkün değil (§5).
4. **Saf geometri (pasif CH):** Rayleigh–Plateau kopma bandı kR < √(d−2).
   - 2B'de bant yok, 6B'de iki kat geniş.
   - d = 3–6'da ölçülen eşikler kuramın %2,5–4 altında (yayvan arayüz düzeltmesi).
5. **Doğrusal olmayan keskin arayüz kuramı** (kimyasal potansiyel μ ile):
   - 2B faz alanı damla yarıçaplarını %1–6 hatayla veriyor. Doğrusal "etkin damla" kuramı %27 fazla veriyordu.
   - Pencerenin üst sınırını açıklıyor: d büyüdükçe boşalma yarıçapı R_boşalma yavaş büyüyor (8,5 → 15,4), uzama hızı ise düşüyor.
   - Alt sınırı sonsuz tüp boğulma eşiğiyle **açıklayamıyor**. Eşik d_c = 2,56 / 2,76 / 2,96 (Δν = −0,002 / −0,001 / 0), Δν ile doğrusal artıyor. Simülasyonda ise bölünme d = 2,1–2,5'te başlıyor. Alttaki bölünme uzamış damlanın uçlarındaki kopmayla oluyor (uç kopması, 3 parça).
6. **Kuantumla ilgisi yok;** tamamen klasik bir faz ayrışması/reaksiyon modeli.

## 1. Yöntem

- **Model:** projenin aktif damla modeli (`qcell.active`, Zwicker vd. 2017 birimleri).
  - Denklem: ∂_t c = ∇²[c(c−1)(2c−1) − ∇²c/4] + s(c).
  - Kimya: ε = ν₊/k = 0,20, k = 0,01; Δν = ν₋ + ν₊ − k taranan iç tüketim fazlası.
- **O(d−1) simetrik indirgeme:**
  - Eksen etrafında simetrik şekiller için d boyutlu denklem, (z, r) yarı düzleminde ∇² = ∂_zz + ∂_rr + (d−2)/r ∂_r ile tam.
  - d herhangi bir gerçel sayı olabilir; bu boyutsal süreklileme d_c'yi bulmaya yarıyor.
  - Ayrıklaştırma: z'de Fourier; r'de hücre merkezli sonlu hacim (yüz alanı r^(d−2), tam hücre hacmi). Kütle makine hassasiyetinde korunuyor.
  - Zaman adımı ActiveCH'deki yarı örtük adım.
  - **Sınır:** yalnız ilk olay (bölünme ya da boşalma) güvenilir. Sonrasında kızlar simetri yüzünden eksene dik uzayamıyor, basık büyüyüp halkaya dönüşüyor.
- **Kartezyen d boyutlu çözücü (`cart_run.py`):**
  - ActiveCH ile aynı denklem ve aynı adım; float32, meshgrid'siz.
  - Periyodik kümeleri union-find ile sayıyor. 16⁶ kutu < 1 GB bellek istiyor.
- **Keskin arayüz kuramları:**
  - **Doğrusal:** seyreltik ve yoğun fazda doğrusallaştırılmış tepkime-difüzyon, Gibbs–Thomson μ = H/12, ara yüz hızı = akı farkı. d boyutta modifiye Bessel fonksiyonları r^(−ν)I_(ν+l), r^(−ν)K_(ν+l), ν = m/2 − 1. m = d damla, m = d−1 tüp.
  - **Doğrusal olmayan:** her faz kimyasal potansiyelde çözülüyor (0 = ∇²μ + s(c(μ)), c(μ) tam dallardan).
    - Uzak alan c = 0,20 spinodale (0,211) çok yakın. Orada etkin difüzyon f'' = 0,04, yani doğrusal kuramın varsaydığının 25'te biri. Doğrusal olmayan kuram bunu kendiliğinden içeriyor.
    - Yoğun dalın var olmaması "içi boşalma" ölçütünü veriyor.

## 2. Doğrulama

| Kontrol | Sonuç |
|---|---|
| Kütle korunumu (eksenel, d = 2–6) | bağıl kayma ≤ 10⁻¹³ |
| d = 2 eksenel ↔ ayna simetrik Kartezyen 2B | fark 0,0080 (h = 0,5) → 0,0020 (h = 0,25), 2. mertebe yakınsama |
| Pasif RP eşiği (tüp R = 8) | ölçülen kR_c = 0,962 / 1,372 / 1,686 / 1,949; kuram √(d−2) = 1 / 1,414 / 1,732 / 2 (d = 3–6) |
| Kuram özdeşlikleri | σ₁ (öteleme) = 0 (≤ 10⁻¹²); σ₀ = dv/dR (10⁻⁶ doğrusal, 10⁻⁴ doğrusal olmayan); pasif sınırda tüp eşiği tam √(d−2) |
| Doğrusal olmayan kuram ↔ 2B faz alanı (Exp11, 5 kararlı damla) | oran 0,99–1,06. Doğrusal kuram 9 damlada 0,70–0,77 |
| d = 3 eksenel ↔ Kartezyen 3B 64³ (Exp11 b3d), Δν = −0,001 | ilk bölünme t = 1850 / 1750–2000; kızlar 9,7 / 9,9 |
| aynı, Δν = −0,004 | ikisi de içi boş kabuk. Kartezyende merkez c 0,78 → 0,03; eksenelde t = 1200 |
| Izgara yakınsaması h = 0,5 → 0,25 | d = 3: bölünme t = 1825 / 1875, R = 9,65 / 9,60. d = 4: boşalma t = 750 / 750, R = 13,63 / 13,65 |
| Kaba ızgara uyarısı | h = 1'de d = 4 yanlışlıkla bölünüyor (t ≈ 2000, 3 parça). Arayüz (w = 1) çözülmeden sonuç tersine dönebiliyor |
| **6B hiperküp 16⁶ ↔ O(5) indirgeme** (Zwicker kimyası, R₀ = 4,5, t ≤ 40) | R_eş farkı her noktada ≤ %2,5 (h = 1); h = 0,5 ile t = 40'ta 6,741 / 6,771 |
| Testler | `tests/test_exp14.py` 12/12 |

## 3. Sonuçlar

### 3.1 Boyut ve iç tüketime göre ilk olay (`outcomes.json`)

Eksenel, h = 0,5, R₀ = 5, uzama tohumu ecc = 0,05; alan teorik boyutun birkaç katı.

| Δν | d → sonuç (t, parça sayısı) |
|---|---|
| −0,004 | 2: şerit (t = 12 000'e kadar olay yok) · 2,25: bölünme (2450; 2) · 2,5: bölünme (2350; 2) · 3: **kabuk** (1200) |
| −0,001 | 2, 2,1, 2,2: şerit · 2,3: bölünme (5000; 3) · 2,4: (3250; 3) · 2,5: (2900; 3) · 3: (1850; 2) · 3,5: (1750; 2) · 3,75: **kabuk** (1000) · 4: kabuk (1000) · 5: kabuk (≤ 1000) · 6: kabuk (500) |
| +0,004 | 2,25: yuvarlak, kararlı · 2,5: bölünme (9850; 2) · 2,75: (2200; 2) · 3: (6650; 3) · 4: (1650; 2) · 5: **kabuk** (750) |
| +0,012 | 6: **bölünme** (1400; 2) |
| +0,022 | 6: **bölünme** (3400; 2, kızlar R = 6,8) |

### 3.2 Pencerenin üst sınırı: neden kabuk?

Doğrusal olmayan kurama göre (Δν = −0,001):
- 3B'de durağan kompakt damla yok. Damla büyürken iki şey yarışıyor: uzayıp bölünmek ya da çekirdeğin spinodalin altına düşmesi.
- Yoğun çekirdeğin var olabildiği en büyük yarıçap R_boşalma = 8,5 / 9,6 / 10,6 / 11,5 / 12,3 / 13,9 / 15,4 (d = 2 / 2,5 / 3 / 3,5 / 4 / 5 / 6).
- Büyüyen damlanın göreli uzama hızı σ₂ − v/R, R = 10'da: +1,7·10⁻³ (3B), 0 (4B), −1,5·10⁻³ (5B), −3,1·10⁻³ (6B).
- **Okuma:**
  - 3B'de damla R ≈ 10'da yavaşlarken (v ≈ 4·10⁻⁴) uzama güçlü; boşalmadan bölünüyor.
  - d ≥ 4'te uzama zayıf ya da negatif ve damla hızla büyüyor; önce boşalıyor.
  - Faz alanı aynı sınırı veriyor: 3,5'te bölünme, 3,75'te kabuk.
- Δν arttıkça (güçlü iç tüketim) damlalar küçük kalıyor ve boşalma geç geliyor. Pencere bu yüzden yukarı kayıyor; Δν ≥ +0,012'de 6B'yi içine alıyor.

### 3.3 Pencerenin alt sınırı

- 2B'de (ayna simetrik) uzayan damla bir şerit oluyor ve kopmuyor: d = 2 ve 2,1–2,2'de t = 12 000'e kadar olay yok.
- Doğrusal olmayan kurama göre sonsuz tüpün boğulma eşiği d_c = 2,56 (Δν = −0,002), 2,76 (Δν = −0,001) ve 2,96 (Δν = 0). Δν = +0,004'te 3B tüp bu kurama göre kararlı. Bu, yalıtık 3B damlanın uzun tüp olup geç bölünmesiyle (§3.4) tutarlı.
- Simülasyonda bölünme daha önce, d = 2,3'te başlıyor: şişkin uçların arkasından kopma ve ortada bir uydu, yani 3 parça.
- Sonuç: alt sınırı sonlu uzunlukta uç kopması belirliyor, sonsuz tüp kararsızlığı değil. Kuram yalnız üst sınır için tahmin gücüne sahip.

### 3.4 Kutu boyu etkisi (Exp11'e not)

- Δν = +0,004'te 3B Kartezyen 48'lik kutuda bölünme t ≈ 2000'deydi. Eksenel çözücü aynı periyodik boyla (Lz = 48) bunu yeniden üretiyor: t = 2225, kızlar 6,23 / b3d 6,14.
- Uzun alanda (Lz = 88) ise damla uzun bir tüp oluyor (boy ~80) ve ancak t = 6625'te üçe ayrılıyor.
- Yani zayıf sürüşte küçük periyodik kutu bölünmeyi öne çekiyor. Δν = −0,001'de kutu boyu sonucu değiştirmiyor (1800 / 1850).

### 3.5 Gerçek 6B ızgara neden yalnız doğrulama için?

- **Pasif CH (yüzey gerilimi):** RP kopması için R/L < √(d−2)/2π = 0,32 gerekiyor. 6B'de bu tüpün hacim oranını ≤ %1,7'ye sınırlıyor.
  - Bu kadar az yoğun faz, Gibbs–Thomson konsantrasyonunu (~0,08) karşılayamıyor ve tüp tamamen çözülüyor.
  - Eksenel denemede: R = 4 tüp, L = 16 eşdeğeri kutu; 3B/4B'de kararlı, 5B/6B'de t = 100–200'de çözüldü ("boyut laneti").
- **Aktif CH:**
  - Doğrusal olmayan kuramla 64 kimya tarandı; hiçbirinde 16⁶–20⁶ kutuya sığan (R ≤ 7) kararlı 6B damla yok.
  - Sebep iki sınır. Küçük R'de Gibbs–Thomson engeli var: μ_I = (d−1)/(12R), R ≤ 3'te seyreltik spinodali aşıyor. Büyük R'de çekirdek boşalıyor.
  - Hızlı kimya (k = 0,1–0,2) ile küçültmek de olmuyor: dış besleme boyu √(f''(ε)/k) arayüzden ince kalıyor ve damlalar çözülüyor.
  - Zwicker kimyasında 6B yapılar 15–25 birim. Tam ızgara için ≥ (60–100)⁶ ≈ 10¹¹–10¹² nokta gerekir.
- **Bu yüzden:** 6B'de bölünme ve kabuk sonuçları O(5) indirgemeden. İndirgeme gerçek 6B hiperküpte (16⁶) doğrulandı; ama orada simetri bozan modlar (eksen dışı uzama, kabuğun parçalanması) yok.

## 4. Literatür ve yenilik

- **Bilinenler:**
  - d boyutta Rayleigh–Plateau, kR < √(D−2): [Cardoso & Dias 2006](https://arxiv.org/pdf/hep-th/0602017), [hep-th/0610004](https://arxiv.org/pdf/hep-th/0610004).
  - 2B'de aktif damlanın bölünmeyip dambılda takılması: [2505.20028](https://arxiv.org/pdf/2505.20028), [Bauermann tezi 2022](https://www.pks.mpg.de/fileadmin/user_upload/MPIPKS/group_pages/BiologicalPhysics/dissertations2/JonathanBauermann2022.pdf).
  - 3B'de kabuk ile bölünme rekabeti ve kabuğun çekirdek spinodal kararsızlığından doğması: [Bauermann vd. 2023, PRR](https://arxiv.org/pdf/2306.10852); deney [Bergmann vd. 2023, Nat. Commun.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10582082/).
  - Matematiksel analiz: [Garcke vd. 2025](https://arxiv.org/pdf/2503.09581), [Garcke vd. 2026](https://arxiv.org/abs/2601.11155). Radyal çözümler, kararlılık, kabuk tipi çözümler, bölünme. Genel d'yi ele alıp almadıkları özetten anlaşılmıyor; tam metin kontrol edilmedi.
  - Model ve 3B bölünme: [Zwicker vd. 2017](https://arxiv.org/abs/1603.01571).
- **Muhtemelen yeni:**
  - (i) Boyutun sürekli değişken yapılmasıyla bölünme penceresi ve pencerenin Δν ile yukarı kayması.
  - (ii) Yüksek boyutta bölünmenin daha güçlü sürüş istemesi.
  - (iii) d boyutta μ biçiminde doğrusal olmayan keskin arayüz kuramı (2B faz alanında %1–6).
- **Yayınlanabilirlik:** kısa bir not olarak mümkün görünüyor, ama önce şunlar yapılmalı:
  - Garcke 2026 ve Bauermann 2023'ün tam metni okunmalı (genel d ve faz diyagramları);
  - ε taranmalı;
  - tohum ve gürültü sağlamlığı sınanmalı;
  - en azından 4B'de tam Kartezyen bir kontrol yapılmalı (h = 0,5 ile 112⁴ ≈ 1,6·10⁸ nokta; bu makinede yok).

## 5. Sınırlar

- Eksenel indirgeme yalnız O(d−1) simetrik şekilleri içeriyor. İlk olaydan sonrası, eksen dışı uzama, kabuğun simetri bozan parçalanması ve sarmal/yılanlanma modları kapsam dışı.
- 2B Kartezyende (Exp11) solucanlar kıvrılıyor ve Δν = −0,004'te ağ biçiminde kopuyor. Bu, ayna simetrik d = 2 sonucundan farklı bir yol.
- Tek ε (0,20) ve tek k (0,01). Pencere sınırları 0,1–0,25 d aralığına kadar çözüldü. Tek tohum.
- Kaba ızgara (h = 1) d = 4'te yanlış sonuç veriyor. Bütün sonuçlar h = 0,5'te; d = 3 ve 4 h = 0,25'te doğrulandı.
- 2B'de (Δν = +0,004, d = 2,25) "olay yok" sonucu sınırda: damla küçük (R = 4,5) ve yuvarlak kalıyor.
- Kabuk algılaması kapalı seyreltik bölge (c < 0,5) ölçütüne dayanıyor. Boşalma zamanları 50–500 birim çözünürlükte.

## 6. Önceki raporlara düzeltmeler

- **Exp11, 3B Δν = −0,004:** "tek damla, R = 19,8, kutu sınırlı" değil, içi boş bir sıvı kabuk. `summary_exp11.md` düzeltildi.
- **Exp11, 3B Δν = +0,004:** t ≈ 2000'deki bölünme 48'lik periyodik kutunun etkisini taşıyor. Yalıtık damlada bölünme çok daha geç (t ≈ 6650) ve 3 parça.

## 7. Bağlantılı deneyler (aynı istek)

- **Exp15:** faz uzayı, Vlasov–Poisson PIC; `results/exp15_phase_space/summary_exp15.md`. Plazma delikleri 1B'de yalnız birleşiyor, 2B/3B'de kıvrılıp parçalanıp sönüyor; kalıcı bölünme yok. Eşik 1B↔2B.
- **Exp16:** veriden denklem, PySINDy; `results/exp16_equations/summary_exp16.md`.

## Kaynaklar

Yukarıdaki §4 bağlantıları; ayrıca [Hutchinson 2018](https://arxiv.org/pdf/1804.08594) (Exp15 karşılaştırması).
