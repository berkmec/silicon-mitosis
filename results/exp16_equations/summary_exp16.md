# Deney 16: Veriden denklem türetme (seyrek regresyon, PySINDy / PDE-FIND)

Tarih: 2026-10-15. Kod: `experiments/exp16_equations/` (`common.py`, `s1a…s1d`, `s1a2`, `s1b2`, `s2a…s2e2`, `fig_pareto_gs.py`, zincir betikleri `run_rest.sh`, `run_fix.sh`). Sonuçlar: bu klasördeki `s*.json`, şekiller `fig_*.png`, günlükler `log_*.txt`. Araç: [PySINDy](https://github.com/dynamicslab/pysindy) 2.1.0 (PDELibrary, WeakPDELibrary, STLSQ, SmoothedFiniteDifference). PDE-FIND yöntemi: [Rudy vd. 2017, arXiv 1609.06401](https://arxiv.org/abs/1609.06401). Ayrı ortam `env/venv_sindy` (numpy 2.4.6, scipy 1.17.1, scikit-learn 1.9.1, matplotlib 3.11.2); ana `env/venv`'e dokunulmadı. Süreçler düşük öncelikli ve tek tek çalıştı.

Not: Analizi bir alt ajan yaptı. Rapor metni ajanın teslim mesajından ana oturum tarafından kaydedildi. Üç değişiklik var: §3.3'e Exp14 düzeltmesi eklendi; §2.2'de k1/k2 değerlerinin kaynağı düzeltildi (bunları kullanıcı değil, görev tanımı verdi); §0'daki BSIM4 gövde etkisi satırı buna göre düzeltildi.

## 0. Kısa cevap

**Önce boru hattı bilinen denklemlerde sınandı. Kısa özet: temiz simülasyon verisinde tam çalışıyor, gürültüde çalışmıyor.**

| Bilinen yasa (doğrulama) | Sonuç |
|---|---|
| Gray–Scott PDE (110 aday terim) | İntegratörle tutarlı türevde iki denklem tam bulundu, katsayı hatası ≤ 2·10⁻¹⁴. Merkezi fark (s = 1, 2): tam destek, hata ≤ %0,19, 200 adım tahminde v hatası %0,16. s = 5'te yalnız u bulundu. %1 gürültüde güçlü form çöktü. Zayıf form bu kurulumda gürültüsüz veride bile başarısız |
| BSIM4 darbe iyonlaşması (90 eğri) | 90/90 eğride tam olarak {1, ln x, 1/x} seçildi. β₀ = 13,8499 V (kart 13,8500), ln a = 0,08421 (kart 0,08419), ln x katsayısı 1,0000002. İlk deneme 0/90 idi; nedeni V₀ ızgarasıydı (§2.2) |
| BSIM4 gövde etkisi | Yalnız yaklaşık: φ = 0,852 V (beklenen 0,819), √ katsayısı 0,663 (bin-62 beklentisi 0,601). Görev tanımındaki k1 = 0,5409, k2 = −0,0267 başka bir bine (L ≥ 8 µm) ait |
| 2B HF sıkıştırılabilirlik (exp10 B) | {1, n^−1/2} bulundu. Katsayı −1,65335 (analitik −1,65334, bağıl hata 6·10⁻⁶). Dışkestirimde R² = 1,000000 |
| Sonlu bariyerli kuyu E1(L) (exp1) | E1 = A_sd/(L + L₀)². L₀ boya göre 0,596–0,629 nm. BenDaniel–Duke beklentisi 0,603 nm. Parametresiz yasa L ≥ 2 nm'de %0,5 içinde |
| C ∝ a (exp1) | **Test edilemez.** FEM yalnız a = 1'de çözülmüş, C(L) tablosu ölçeklemeyle üretilmiş; eğim 1 tanım gereği çıkıyor |
| FEM küre yakınsaması C(h) | **Başarısız.** 6 noktadan üs bulunamadı (C∞ serbestken p = 2,9 ± 0,4). C_exact verilirse p = 1,91 |

**Kapalı formu bilinmeyen verilerden çıkan aday denklemler** (hepsi ayrı tutulmuş veriyle doğrulandı; ayrıntı §3):

1. **exp11 2B aktif damla:** R* = 0,233·ε/ν₋, ν₋ = Δν + k(1−ε). Birini dışarıda bırakma hatası %3,0 (14 durağan damla). Bu bir akı dengesi: ara yüzün birim uzunluğuna düşen akı j ≈ 0,117·ε. Çözünme sınırı R* ≈ 2,6 civarında. Ama aynı öngörülen R* = 2,59'da ε = 0,08 koşusu çözünüyor, ε = 0,14 koşusu yaşıyor; kritik yarıçap ε'a da bağlı.
2. **exp13 filaman:** w_toplam = I / j_f, j_f = 2,55 ± 0,03 µA/µm. Uçtan beslemeli ve dağıtık balastlı 12 koşuda seyrek regresyon yalnız {I} terimini seçti; aşağı tarama doğrulamasında R² ≥ 0,9988. Dağıtık balastta V_ds 0,803–0,823 V'ta kenetli; exp12'nin Maxwell değeri 0,8018 V. **Bölünme genişliği λ'ya orantılı değil.** λ = 40 µm'de bölünme w ≈ 156–164 µm'de (≈ 4λ) oldu. λ = 5–20 µm'de filaman 203 µm'ye bölünmeden ulaştı.
3. **exp11 3B bölünme:** R_kız ∝ ν₋^(−0,87 ± 0,05), R_ebeveyn ∝ ν₋^(−0,77). Yalnız 4 nokta var. Üs, akı dengesinin iki uç değeri (−1/2 ve −1) arasında.
4. **exp8 fisyon eşiği:** Eşik, kazanç eğrisinin sıfır geçişinden sürekli okununca g_c ∝ N^p bulundu: p = −1,23 (d = ∞), −0,98 (d = 4), −1,03 (d = 2). Exp8'deki −0,73 ve −0,90 ızgara nicemlemesinin ürünü. Kapının üssü düzleştirdiği verilerle desteklenmiyor.
5. **Gray–Scott ölçeklemesi:** leke sayısı N ∝ D^(−1,03), leke alanı ∝ D^0,95. Boyut analizinin beklentisi −1 ve +1. N(t) çoğu koşuda lojistik, ama bütün koşular için tek bir büyüme yasası seçilemedi.

## 1. Yöntem

- **Seyrek regresyon:** PySINDy STLSQ, `unbias=True`. Kütüphane sütunları ve hedef RMS'leri 1 olacak biçimde ölçeklendi. Böylece eşik τ boyutsuz ve örnek sayısından bağımsız: |ξ_j|·rms(θ_j) < τ·rms(y) olan terim atılıyor. Her veri setinde τ log-ızgarada tarandı (eşik yolu).
- **Model seçimi yalnız ayrı tutulmuş veriyle yapıldı; doğru model yalnız sonradan puanlamada kullanıldı.** İki kural var:
  - Kural A: doğrulama bağıl RMSE'si en iyinin 1,05 katı içinde kalan en seyrek model.
  - Kural B: doğrulama R²'si en iyiden en çok 10⁻³ düşük olan en seyrek model.
  - Gürültüsüz BSIM4 tablolarında kural B'nin toleransı 10⁻⁶ alındı.
  - Neredeyse sabit hedeflerde R² anlamsız olduğu için orada kural A kullanıldı.
  - GS'de iki kuralın sonucu da rapordadır. **Kural A, merkezi fark türevindeki sistematik hatayı ek terimlerle emdi ve 47/50 terimli modeller seçti.**
- **Doğrulama türleri:** ayrı simülasyon koşusu (GS), dışkestirim (BSIM4: V_ds'nin üst %30'u; exp10: n < 10¹⁰), histerezisin aşağı tarama dalı (exp13), birini/bir koşulu dışarıda bırakma (küçük veri), bootstrap ve eşik duyarlılığı. Eşik duyarlılığı, doğru desteği veren τ aralığı olarak raporlandı.
- **Fizik bilgisi açıkça kullanılan yerler:**
  - BSIM4 değişkenleri (V_dseff yumuşatması, V_bseff kenetlemesi);
  - exp11 damlasında ν₋ = Δν + k(1−ε);
  - E1(L)'de sert duvar katsayısı A_sd.
  Bunlar "keşif" değil, aday kütüphanenin seçimidir.

## 2. Doğrulama: bilinen denklemleri geri bulma

### 2.1 Gray–Scott PDE-FIND (`s1a_grayscott.json`, `s1a2_weak_reduced.json`, `fig_s1a_grayscott.png`, `fig_s1a_pareto.png`)

Veri `qcell.grayscott` ile yeniden üretildi: n = 128, F = 0,026, k = 0,059, Du = 2·10⁻⁵, Dv = 10⁻⁵, h = 2,5/256, dt = 1, her adım kaydedildi.
- Eğitim: tohum 1, t = 1500–1580, 8 leke.
- Doğrulama: tohum 2, t = 2500–2580, 24 leke.
- Kütüphane: (u, v) polinomları (3. dereceye kadar), bütün 1. ve 2. uzay türevleri (u₁₂ dahil) ve bunların çarpımları; denklem başına 110 aday.
- Gerçek destek: u_t için 5 terim, v_t için 4 terim.

| Varyant | Seçilen terim (u, v), kural B | Tam destek | En büyük bağıl katsayı hatası (u, v) | Doğru desteği veren τ aralığı | 200 adım v hatası |
|---|---|---|---|---|---|
| ileri fark, s = 1 (integratörle tutarlı) | 5, 4 | ✓ ✓ | 4·10⁻¹⁵, 2·10⁻¹⁴ | 10⁻⁴ – 0,63 | 4·10⁻¹⁴ |
| merkezi fark (PySINDy FD), s = 1 | 5, 4 | ✓ ✓ | %0,099, %0,18 | 0,35 – 0,63 | %0,16 |
| merkezi fark, s = 2 | 5, 4 | ✓ ✓ | %0,11, %0,19 | 0,11 – 0,63 (u), 0,35 – 0,63 (v) | %0,16 |
| merkezi fark, s = 5 | 5, 32 | ✓ ✗ | %0,13, — | u: 0,35 – 0,63; v: yok | %0,28 |
| merkezi fark, s = 1, %1 gürültü | 90, 99 | ✗ ✗ | — (R²_val = 0,03) | yok | ıraksadı |
| zayıf form, 110 terim, %1 gürültü | 101, 73 | ✗ ✗ | — (R²_val = 0,988) | yok | %8,6 |
| zayıf form, 110 terim, s = 2, %5 gürültü | 107, 110 | ✗ ✗ | — | yok | ıraksadı |
| güçlü form, 20 terim (çarpımsız), gürültüsüz | 5, 4 | ✓ ✓ | %0,096, %0,17 | 0,034 – 0,63 (u), 0,11 – 0,63 (v) | — |
| zayıf form, 20 terim, gürültüsüz / %1 / %5 | 12/11/6, 7/7/12 | ✗ | v'nin gerçek terimleri %2 içinde, 3 küçük sahte terim | yok | — |
| güçlü form, 20 terim, %1 / %5 gürültü | 13 / 10 | ✗ | R²_val ≈ 0,01 | yok | — |

Merkezi fark s = 1'de bulunan denklemler (kural B):
`u_t = 0,02599 − 0,02599·u − 0,99952·u v² + 1,998·10⁻⁵·(u₁₁ + u₂₂)` ve
`v_t = −0,08495·v + 0,99937·u v² + 0,9982·10⁻⁵·(v₁₁ + v₂₂)` (gerçek: 0,026; −0,026; −1; 2·10⁻⁵ / −0,085; +1; 10⁻⁵).

Okuma:

1. **Hata sınırını türev tahmini belirliyor, regresyon değil.** İleri-Euler verisinin merkezi farkı her zaman t − ½'de ortalanıyor. Bu yüzden s = 1 ile s = 2 aynı ~%0,1–0,2 sapmayı veriyor. Difüzyon katsayıları sistematik olarak −%0,1 / −%0,18 çıkıyor.
2. **Kural A bu sistematik hatada aşırı uyuyor.** Seçtiği 47/50 terimli modelin katsayı hataları %56–72. Gürültüsüz veride bile seçim kuralı sonucu belirliyor; Pareto eğrisinde (`fig_s1a_pareto.png`) doğru model belirgin bir "dirsek" üzerinde.
3. **Gürültüye dayanıklılık yok.** %1 gürültüde ikinci türev gürültüyü ~1/h² büyütüyor ve güçlü form çöküyor. PySINDy'nin zayıf formu (p = 4, H = 8–12 ızgara aralığı, K = 1000–2000) bu kurulumda gürültüsüz veride bile u denklemini bulamadı. Olası neden: u ≈ 1 arka planında {1, u, u², u³} alt-alan integrallerinde neredeyse eşdoğrusal. Projenin verileri gürültüsüz simülasyon çıktısı olduğu için bu sınır aşağıdaki uygulamaları etkilemiyor. Ölçüm verisine uygulanacaksa önce çözülmeli.

### 2.2 BSIM4 tabloları (`s1b_bsim4.json` ilk deneme, `s1b2_bsim4.json` düzeltilmiş, `fig_s1b_bsim4.png` ilk deneme)

Cihaz W = 1 µm, L = 0,15 µm; kartta bin `model.62`. Beklenen değerler kartın bu bininden hesaplandı (binunit = 2, uzunluklar metre; L_eff = 0,12614 µm, W_eff = 0,95628 µm):
- β₀,eff = 13,850 V;
- a = α₀/L_eff + α₁ = 1,0878 V⁻¹;
- Φs = 0,819 V;
- √ katsayısı K1·√(1 + LPEB/L_eff) = 0,601;
- V_bs katsayısı −K2,eff + K3B·TOXE·Φs/W_eff = 0,130.

**(i) Darbe iyonlaşması, y = ln(I_sub/I_d), kütüphane {1, ln x, 1/x, x, x², 1/x², √x}**

- **İlk deneme başarısız.** x = V_ds − V₀ alındı ve V₀ 0,02 V'luk ızgarada tarandı. 90 eğrinin hiçbirinde {1, ln x, 1/x} seçilmedi; seçilen modeller 4–7 terimliydi.
- **Teşhis:** V₀ ızgarası x'e ≤ 10 mV sistematik hata katıyor. Gürültüsüz veride yoğun modeller bu hatayı emiyor.
- **Düzeltme:** x = V_ds − V_dseff alındı. V_dseff, BSIM4'ün kendi yumuşatması (δ = 0,01). V_dsat her eğride sürekli uyduruldu.
- **Sonuç:**
  - 90/90 eğride tam olarak {1, ln x, 1/x} seçildi. Eğrilerin %100'ünde eşik yolu bu desteği içeriyor.
  - β₀ medyanı 13,8499 V (p10–p90: 13,8480–13,8502).
  - ln a = 0,08421 (beklenen 0,08419).
  - Serbest ln x katsayısı 1,0000002.
  - Dışkestirim RMSE'si (ln) 1,4·10⁻⁶.
  - **Chynoweth/BSIM4 biçimi ln(I_sub/I_d) = ln(a·x) − β₀/x geri bulundu.**
- **Ders:** Doğrusal olmayan bir kaydırma parametresi ızgarada taranırsa seyrek regresyon yanlış yapı "keşfeder".

**(ii) Gövde etkisi vth(V_bs, V_ds), V_bs ≤ 0,7 V**

- **Ön kontrol:** vth'nin V_gs'e bağlılığı ≤ 0,05 mV.
- **İlk deneme:** φ ızgarada tarandı. Seyrek model φ = 1,19'da {1, s, V_bs, V_bs², s·V_bs} seçti; katsayılar birbirini götüren büyük sayılar (−14,9 + 14,4·s + …).
  - Dışkestirimi çok iyi (0,02 mV).
  - Ama yorumlanamaz: tipik bir eşdoğrusallık durumu.
- **Düzeltme denemesi:** BSIM4'ün V_bseff üst kenetlemesi (0,95·Φs) eklendi ve φ sürekli optimize edildi.
  - BSIM4 biçimi {1, √(φ − V_bseff), V_bseff, V_ds, V_ds·V_bseff} şu değerleri verdi: φ = 0,852 V, √ katsayısı 0,663, V_bs katsayısı 0,158, DIBL −8,6·10⁻⁵ V/V.
  - Eğitim RMSE'si 0,04 mV, dışkestirim RMSE'si (V_bs = 0,6–0,7) 0,88 mV.
  - Doğru destek eşik yolunda hiç çıkmadı.
  - V_bs ≥ 0,8 V'ta vth düzleşiyor (0,472 / 0,457 / 0,454 V); bu kenetlemenin kendisi.
- **Sonuç: yalnız yaklaşık geri bulundu.** φ %4, √ katsayısı %10 sapıyor. Eksik olan muhtemelen BSIM4'ün dar kanal terimi: DVT0W = −3,58, X_dep üzerinden V_bs'ye bağlı, mV mertebesinde.
- **Görev tanımındaki k1 = 0,5409, k2 = −0,0267:** bunlar uzun kanal binlerine (`model.0/1`, L ≥ 8 µm) ait. Simüle edilen 0,15 µm cihazda k1 = 0,907 ve bin ölçeklemeli k2,eff = −0,128. LPEB ile etkin √ katsayısı 0,601, veriden çıkan 0,663.

### 2.3 exp10 2B elektron gazı (`s1c_exp10.json`, `fig_s1c_exp10.png`)

Kütüphane: z^p, z = n/10¹² cm⁻², p ∈ {−2 … 2} (9 terim). Eğitim n ≥ 10¹⁰, doğrulama n < 10¹⁰ (dışkestirim).

| Varyant | Seçilen | Doğrulama R² | Not |
|---|---|---|---|
| B: Si, g = 4, HF, kapı ekranı yok | 1,0000 − 1,65335·z^−1/2 | 1,000000 | Analitik 1 − (√g/π)·r_s = 1 − 1,65334·z^−1/2. Doğru destek τ = 10⁻⁷ – 0,11 aralığında. Bootstrap'ta iki terim de her seferinde var. log-log üs −0,5000; (1 − K/K_kin)/r_s = 0,63662 = 2/π |
| C: HF, kapı ekranlı değişim | 2 terim | −384 | kuvvet yasası değil (düşük n'de doyuyor) |
| D, E: korelasyon vekili | 1 + c·z^−1/2 | 0,993 / 0,992 | yaklaşık; sabit 0,82 / 1,28 (bootstrap'ta kararsız) |
| A: tek vadi, DMC korelasyonu | 7 terim | 0,993 | seyrek yasa yok (Attaccalite parametrizasyonu) |

### 2.4 exp1 FEM (`s1d_exp1.json`, `fig_s1d_exp1.png`)

**C ∝ a: test edilemez.**
- exp1'de FEM yalnız a = 1'de, ε₀a biriminde çözülmüş.
- `transport.json`'daki C(L) ölçeklemeyle üretilmiş; log-log eğimi 1,000000 tanım gereği çıkıyor. Bu bir keşif ya da test değildir.

**Küre FEM yakınsaması: başarısız.**
- Veri: C(h), 6 mesh.
- C∞ + A·h^p serbest uydurulunca p = 2,92 ± 0,36 ve C∞ +%0,6 çıktı. Asimptot öncesi h = 0,3 noktası (%12 hata) atılırsa p = 1,19 ve C∞ −%0,5.
- {1, h, h², h³} seyrek kütüphanesi dört terimin hepsini tuttu; C∞ −%2,1.
- C_exact verilirse p = 1,91 (exp1'in değeri).
- Sonuç: bu 6 noktadan yakınsama yasası "keşfedilemez".

**Sonlu bariyerli Si/SiO₂ kübik kuyu E1(L): başarılı.**
- Serbest 3 parametreli A/(L + L₀)^p: A, L₀ ve p yozlaşık; A = 5,10 (sert duvarın 1,17 katı), L₀ = 0,73, p = 2,07.
- A = A_sd = 3π²ħ²/(2m_w) = 4,339 eV·nm² ve p = 2 sabitlenince her boy için L₀ = 0,629 / 0,615 / 0,610 / 0,604 / 0,602 / 0,597 / 0,596 nm (L = 1 … 8 nm).
- Kuram: BenDaniel–Duke sınır koşuluyla derin kuyu limitinde L₀ = 2m_b/(m_w κ), κ = √(2m_b V₀)/ħ. Değerler m_w = 0,26, m_b = 0,5, V₀ = 3,1 eV ile **L₀ = 0,603 nm**.
- **Parametresiz yasanın bağıl hatası:** +%3,2 (1 nm), +%1,1 (1,5 nm), +%0,5 (2 nm), |·| ≤ %0,2 (3–8 nm).
- Küçük L'de L₀'ın büyümesi beklenen etki: E büyüdükçe κ küçülüyor.
- Tek parametreli uydurma L₀ = 0,618 nm verdi; LOO hatası ≤ %2,4.
- Yasa: E1(L) = 3π²ħ²/(2m_w (L + 2m_b/(m_w κ))²).

### 2.5 Kütle korunumu: exp11 yükleme koşuları (`s2e_charging.json`, `s2e2_growth.json`, `fig_s2e_charging.png`)

- Çözücüde kütle korunumu kesin: c_ortalama = c₀ + J·t, dört basamak.
- **dR/dt = f(R) tek değişkenli seyrek modeli başarısız** (dışkestirim R² ≈ −50). R² = z için dz/dt = sabit beklentisi de seçilmedi.
- **Neden:** uzak alan aşırı doygunluğu c_far önce 0,011'den 0,078'e çıkıyor, sonra 0,033'e iniyor. Büyüme yalnız R'ye bağlı değil.
- R²(t) doğrusal uydurmanın R²'si 0,983–0,994; eğim 0,114. Kuram J L²/(π Δc) = 0,108–0,113 veriyor (oran 1,02–1,06).
- c_far ikinci değişken olarak eklenince seçilen model (bir koşu dışarıda doğrulama R² = 0,945, 7 koşu):
  dR/dt = 0,0039 − 6,2·10⁻⁵ R − 0,192/R² + 0,187·c_far/R.
  - Difüzyon sınırlı büyüme (c_far/R) ve Gibbs–Thomson (1/R²) terimleri seçildi.
  - Ama oranları 1,02; basit kuramın γ = 1/12'sinin 12 katı.
  - **Nicel uyum yok: aday model.**

## 3. Uygulama: kapalı formu bilinmeyen veriler

### 3.1 exp11 2B aktif damla (`s2a_droplets.json`, `s2a_dense_runs.json`, `fig_s2a_droplets.png`)

**A. Mevcut koşular (exp11 sınıfı "tek, kararlı" ya da durağan yuvarlak sonlanan "bölündü")**

- 14 durağan yuvarlak damla var. "Solucan" koşuları dışarıda: R_eş ≈ 40–51 ve en-boy ≈ 1 olan kıvrımlı yapılarla bitiyorlar ve otomatik ölçütü yanıltıyorlar.
- **Genel kütüphane 14 noktada kararlı bir seyrek yol vermedi.** Kütüphane {1, ε, ν₋, ε/ν₋, 1/ν₋, ε²/ν₋, √(ε/ν₋), ε/ν₋²}.
  - Yoğun modeller (7–8 terim): LOO %4,2–4,4.
  - Tek terimli {ε/ν₋}: kat içi yeniden seçimle LOO %7,9.
- **Akı dengesi yasası en iyi ve en seyrek model:**
  - Yasa: R* = a·ε/ν₋, a = 0,2333. Birini dışarıda bırakma hatası %3,0, R² = 0,982.
  - Bu, 2πR·j = ν₋πR² dengesidir. 2j/ε = 0,217–0,242 (14 koşu, ±%5) olduğundan birim ara yüz uzunluğuna düşen akı j ≈ 0,117·ε, R'den neredeyse bağımsız.
  - 2B'de R* ~ ℓ = 10 bölgesinde bu, R* ∝ 1/ν₋ verir.
- **Çözünme:**
  - Çözünen 6 koşunun bu yasayla öngörülen R*'ı 1,41–2,59. Durağan 14 koşununki 2,59–5,83; gözlenen en küçük durağan R* 2,65.
  - Ayrım bir çiftte bozuluyor: ε = 0,08, Δν = −0,002 (çözündü) ile ε = 0,14, Δν = +0,004 (yaşadı) aynı R* = 2,59'u veriyor.
  - Kritik yarıçap ≈ 2,6 w ama ε'a da bağlı.
- **Uzama:** Yuvarlak kalmayan 10 koşunun öngörülen R*'ı ≥ 5,5 (7 terimli modelle). Ama ε = 0,08'de 5,86 yuvarlak kaldı. Şekil kararsızlığı eşiği ε'a bağlı.
- **Gevşeme hızı:** λ/(ν₋/2) medyanı 1,25 (p10–p90 0,39–1,45, n = 32). Örnekleme 250 t₀ aralıklı ve gevşeme 1–3 örnekte bitiyor, bu yüzden kaba.

**B. Yeni yoğun koşular** (aynı çözücü, L = 64 w, R(t) her 5 t₀'da; 6 (ε, Δν) koşulu × R₀ = 3 ve 8)

- **Mevcut veri yeniden üretildi.** R₀ = 3 ve 8 aynı R*'a gidiyor; R* her iki yönden de çekici.
- Eşik yarıçapları exp11 ile 0,01 içinde: 5,015 / 3,909 / 3,090 / 2,646 / 2,985 / 4,687.
- **Seyrek dR/dt modeli başarısız.**
  - 15 terimli {1, R, 1/R, R², 1/R²} × {1, ε, ν₋} kütüphanesinde STLSQ hiçbir τ'da terim atmadı.
  - Bir koşul dışarıda doğrulama R² = −22.
  - Neden: eşdoğrusallık, ve her koşulda (R, dR/dt) tek bir monoton eğri.
- Zorlanmış en küçük model: ν₋R katsayısı −0,63; akı dengesi beklentisi −1/(2Δc) ≈ −0,5. İşaret ve mertebe uyuyor, sınama keskin değil.

### 3.2 exp13 akım filamanları (`s2c_filaments.json`, `fig_s2c_filaments.png`)

Kütüphane {1, I, I², √I}. Eğitim yukarı tarama, doğrulama aşağı tarama (histerezis dalı, bağımsız durumlar).

| Koşu grubu | Seçilen | j_f = 1/eğim (µA/µm) | Doğrulama R² | V_ds |
|---|---|---|---|---|
| B, uçtan beslemeli (8 koşu, r_D 30–1000 Ω/µm) | 7/8'inde yalnız {I} | 2,530 – 2,610 | ≥ 0,9988 ({I} seçilen 7 koşu); B_W100_rD100'de 4 terim, 0,958 | kenetli değil: medyan 0,86–1,71 V, hat düşümüyle artıyor |
| C, dağıtık balast, R_L = 5 kΩ (λ = 5/10/20/40, ince adımlı tekrar dahil 5 koşu) | hep yalnız {I} | 2,546 – 2,579 | ≥ 0,9997 | kenetli: medyan 0,8027 / 0,8037 / 0,8076 / 0,823 V, aralık boyunca ≤ 8 mV değişim |
| D, R□ = 10¹⁰ (dy 0,05 / 0,1) | 4 terim | 2,57 / 3,01 | 0,999 / 0,968 (araya serpiştirilmiş) | 0,79–0,81 V |
| A, global bağlaşım (5 set) | yetersiz (filaman bölgesinde 1–18 nokta) | 1,7 – 2,9 | negatif (A4: 0,73) | sete göre 0,61–0,72 V |

- **Yasa: w_toplam = I / j_f, j_f ≈ 2,55 µA/µm.** Kesişim ≈ 0, filaman sayısından ve geometriden bağımsız (B ve C). Bu S-NDC'nin "eşit alan" resmidir: akım yalnız filamanın genişlemesiyle artıyor.
- **Kenetlenme gerilimi:** C'de exp12'nin analitik Maxwell değeriyle (0,8018 V) %0,1–2,7 içinde, λ büyüdükçe hafif artıyor.
- **Bölünme ve λ:**
  - Bölünme yalnız λ = 40 µm, R_L = 5 kΩ'da oldu. Ebeveyn genişliği 163,5 µm (kaba adım) ve 155,5 µm (ince adım), yani 3,9–4,1·λ. Kızlar 88 + 78 ve 83 + 72,5 µm.
  - İki ayrı tohumda (W = 300) bölünme V_DD = 3,02 / 3,04 V'ta oldu; kızlar 90 + 83 ve 77 + 98 µm.
  - **λ = 5, 10, 20 µm'de filaman 203 µm'ye (10–40·λ) bölünmeden ulaştı.** Exp13'ün "bölünme ≈ 4λ'da" gözlemi tek bir λ'ya dayanıyor; λ'ya orantılı bir yasa değil.
  - W = 600'de bölünme yok. Bölünme genişliğini λ'nın fonksiyonu olarak türetmeye veri yetmiyor.

### 3.3 exp11 3B bölünme yarıçapı (`s2d_scaling.json`, `fig_s2d_scaling.png`)

ε = 0,2, ν₋ = 0,006 / 0,008 / 0,010 / 0,012 (Δν = −0,002 … +0,004), bölünme t ≈ 2000'de.

| Büyüklük | Değerler | Kuvvet yasası |
|---|---|---|
| R_kız (doğumda) | 11,13 / 8,80 / 7,15 / 6,14 | ∝ ν₋^(−0,865 ± 0,015), bootstrap %5–95: −0,93 … −0,82; LOO ≤ %1,8 |
| R_ebeveyn (bölünmeden hemen önce) | 12,96 / 10,50 / 8,85 / 7,56 | ∝ ν₋^(−0,774 ± 0,020), bootstrap −0,86 … −0,75 |

- Akı dengesine göre üs R ≪ ℓ'de −1/2, R ≫ ℓ'de −1 olmalı (ℓ = √(D/k) = 10 w). Ölçülen değerler arada ve R ≈ ℓ ile tutarlı.
- Sınırlar: yalnız 4 nokta, tek tohum, 64³ kutu. Örnekleme 250 t₀ olduğundan doğumda kızlar zaten büyümüş; 2R_kız³/R_eb³ = 1,07–1,27.
- Δν = −0,004'te bölünme yok. Exp11 bunu "R = 19,8, kutu sınırlı tek damla" diye kaydetmişti. **Exp14 düzeltmesi:** bu yapı içi boş bir sıvı kabuk (merkez konsantrasyonu t = 1000'de 0,78, t = 2000'de 0,03). Yani bu nokta bölünme verisinin dışında kalıyor.

### 3.4 exp8 küme fisyonu (`s2d_scaling.json`)

- **Yöntem:** g_c, E_tek − E_parça kazanç eğrisinin sıfır geçişinden alındı (log g'de doğrusal ara değer). ×1,45'lik g ızgarasının nicemlemesi böylece kalkıyor.
- **Sonuç:**
  - d = ∞: g_c = 0,738 / 0,563 / 0,421 / 0,369 / 0,301 / 0,219 (N = 10 … 26), p = −1,23 ± 0,07 (bootstrap −1,37 … −1,10).
  - d = 4: p = −0,98 ± 0,08 (−1,21 … −0,82).
  - d = 2: p = −1,03 ± 0,14 (−1,33 … −0,87).
  - Exp8'in ızgaradan okuduğu eğimler: −1,08 / −0,73 / −0,90.
- **Okuma:**
  - Üç kapı uzaklığı da sıvı damla (Rayleigh) beklentisi g_c ∝ 1/N ile uyumlu; d = ∞ biraz daha dik.
  - Exp8'deki "kapı eğimi düzleştiriyor" yorumu verilerle desteklenmiyor.
  - LOO hataları %15–24: küçük kümelerin sihirli sayı (LJ13, LJ19) etkileri büyük.

### 3.5 Gray–Scott: leke sayısı ve ölçek (`s2b_gs_count.json`, `s2d_scaling.json`, `fig_s2b_gs_count.png`)

**Ölçek (gs_ramp, D yarıya indirilerek):**
- Plato N = 82,8 / 139,3 / 282,8 / 580,8 (D = 1 / 0,5 / 0,25 / 0,125). N ∝ D^(−0,945 ± 0,05); ilk evre atılınca D^(−1,03).
- Ortalama leke alanı ∝ D^(0,952 ± 0,04).
- Bu, difüzyon boyu ∝ √D boyut analiziyle uyumlu (beklenen −1 ve +1).

**N(t) büyüme yasası:**
- Dört 3 parametreli aday (lojistik, Gompertz, √N cephe, doğrusal doyum) yörünge uydurması, BIC ve 5 katlı çapraz doğrulamayla karşılaştırıldı.

| Seri | En iyi (BIC) | Rakiplerin ΔBIC'i | Lojistik ln2/r |
|---|---|---|---|
| pilotGS (256², L = 2,5) | lojistik, R² = 0,975 | ≥ 33,8 | 830 t |
| GS F0,026 k0,059 (128²) | lojistik, R² = 0,87 | 4,6 / 13,5 / 43,2 | 415 t (exp11 sınıflandırması T_çift = 703, tanımı farklı) |
| gs_track (ölçek değişimi sonrası) | lojistik, R² = 0,77 | 6,2 / 6,6 / 12,7 | 431 t |
| F0,042 k0,063 | lojistik | 1,5 / 1,5 / 7,7 (ayırt edilemiyor) | 1275 t |
| F0,018 k0,055 | √N cephe | 0,1 / 1,4 (ayırt edilemiyor) | — |
| F0,034 k0,063 | doğrusal doyum | ≥ 11,4 | — |
| F0,010 k0,050 | hepsi eşdeğer (R² = 0,54) | ≤ 0,5 | — |

- Sonuç: en yaygın aday lojistik (dN/dt = rN(1 − N/K)), ama bütün (F, k) için tek bir yasa seçilemiyor.
- PySINDy'nin doğrudan kullanımı başarısız. {N, N², √N, N^1,5} kütüphanesi [0, 1]'de eşdoğrusal; hiçbir eşikte terim atılmadı ve sonuç bilgi vermiyor.

## 4. Özet tablo

| Veri | Denklem | Katsayılar (belirsizlik) | Doğrulama | Bilinen yasa? | Uyarı |
|---|---|---|---|---|---|
| GS (yeniden üretildi) | u_t = Du∇²u − uv² + F(1−u); v_t = Dv∇²v + uv² − (F+k)v | s = 1 merkezi farkta hepsi ≤ %0,18; doğru destek τ ∈ [0,35; 0,63] | ayrı koşu R² = 0,99996; 200 adım tahmin %0,16 | evet, geri bulundu | gürültüde başarısız; seçim kuralına duyarlı |
| BSIM4 I_sub | ln(I_sub/I_d) = ln(a·x) − β₀/x, x = V_ds − V_dseff | β₀ = 13,8499 (p10–p90 ±0,002), ln a = 0,0842 | dışkestirim RMSE 1,4·10⁻⁶ | evet (90/90 eğri) | V_dseff tanımı BSIM4'ten alındı |
| BSIM4 vth | c₀ + c₁√(φ − V_bseff) + c₂V_bseff + DIBL | φ = 0,852, c₁ = 0,663, c₂ = 0,158 | dışkestirim 0,88 mV | kısmen (φ %4, c₁ %10 sapma) | dar kanal terimi eksik; k1 = 0,5409 farklı bin |
| exp10 HF (B) | K/K_kin = 1 − (2/π)·r_s | −1,65335·z^−1/2 (hata 6·10⁻⁶) | dışkestirim R² = 1 | evet | — |
| exp1 E1(L) | E1 = 3π²ħ²/(2m_w (L + L₀)²) | L₀ = 0,596–0,629 nm (kuram 0,603) | parametresiz hata ≤ %0,5 (L ≥ 2 nm) | evet (BenDaniel–Duke) | 1 nm'de +%3 |
| exp1 C ∝ a; C(h) | — | — | — | — | dairesel / veri yetersiz |
| exp11 yükleme | dR/dt ≈ 0,187 c_far/R − 0,192/R² + … | 4 terim | koşu dışarıda R² = 0,945 | kısmen (terimler doğru, oran 12×) | aday |
| exp11 2B damla | R* = 0,233·ε/ν₋ | a = 0,2333; 2j/ε = 0,22–0,24 | LOO %3,0, R² = 0,982 | akı dengesi biçimi | j(ε) ampirik |
| exp11 2B çözünme | öngörülen R* ≲ 2,6 → çözünür | — | 19/20 net, 1 çift aynı R* (2,59) | — | kritik yarıçap ε'a bağlı |
| exp13 filaman | w_toplam = I/j_f | j_f = 2,55 ± 0,03 µA/µm | aşağı tarama R² ≥ 0,9988 (12 koşu) | S-NDC eşit alan resmi | A setlerinde veri yetersiz |
| exp13 kenetlenme | V_ds ≈ V_M | 0,803–0,823 V (C) | aralık boyunca ≤ 8 mV değişim | exp12 Maxwell 0,8018 V | B'de hat düşümü yüzünden yok |
| exp13 bölünme | w_div ≈ 4λ (yalnız λ = 40) | 155–164 µm | 3 tohum | — | λ = 5–20'de bölünme yok: yasa değil |
| exp11 3B | R_kız ∝ ν₋^−0,87; R_eb ∝ ν₋^−0,77 | ±0,05 / ±0,06 (bootstrap) | LOO ≤ %2,6 | üs akı dengesi sınırları arasında | 4 nokta |
| exp8 fisyon | g_c ∝ N^p | p = −1,23 / −0,98 / −1,03 (d = ∞/4/2) | LOO ≤ %24 | Rayleigh −1 | 6 nokta, sihirli sayılar |
| GS ölçek | N ∝ D^−1, alan ∝ D | −1,03 ± 0,01 / 0,95 ± 0,04 | 4 nokta | boyut analizi | — |
| GS N(t) | çoğunlukla lojistik | ln2/r = 415–1275 t | BIC + 5 katlı CV | — | koşular arasında tutarsız |

## 5. Kullanılmayan veya yetersiz veriler

- **exp5 `mos1d.json`:** yalnız iki çalışma noktası (nominal ve kalibre), tarama yok; regresyon yapılamaz.
- **exp6 RDF:** yalnız iki katkı yoğunluğu var. İki noktalı üs σ_VT ∝ N_A^0,16 (MC) ve N_A^0,27 (Stolk) çıkıyor; bu bir yasa değil. Nokta oranları geri çekildi ve kullanılmadı.
- **exp8 `pf_validation.json` ve `selection.json`:** analitik kontroller zaten yapılmış; λ(d) için koşul başına 5 nokta var. Kalite için alınmadı.
- **exp12 indirgenmiş model ve uzaysal taramalar:** çoğunlukla kategorik sonuçlar (bistabil mi değil mi). `maxwell.json` tek bir analitik değer.
- **exp11 A koşulları:** büyük ölçüde bölünme yok; yalnız A2 yükleme koşuları kullanıldı (§2.5). `b3d.json` (makale hızları) parametre kaydı içermediği için dışarıda.

## 6. Sınırlar ve dürüstlük notları

1. **Gürültüye dayanıklılık yok.** Gray–Scott'ta %1 gürültü güçlü formu çökertti. PySINDy'nin zayıf formu bu kurulumda (p = 4, H = 8–12, K = 1000–2000) gürültüsüz veride bile u denklemini bulamadı. Ölçüm verisine uygulanmadan önce test fonksiyonu ve alt-alan seçimi ayarlanmalı.
2. **Seçim kuralı sonucu değiştiriyor.** Kural A, sistematik türev hatasını emen yoğun modeller seçti. Kural B (ΔR² ≤ 10⁻³) bu çalışmada GS'nin ilk sonuçları görüldükten sonra eklendi; iki kuralın sonucu da raporda. BSIM4'te kullanılan 10⁻⁶ toleransı ve s1b2'deki V_dsat/φ düzeltmeleri de ilk başarısızlıktan sonra geldi. Teşhis gerekçeleri §2.2'de.
3. **Fizik bilgisi kütüphaneye girdi.** BSIM4 değişkenleri, ν₋, A_sd ve c_far gibi büyüklükler kütüphaneye konuldu. Bu yerlerde sonuç "yapı keşfi" değil, "verilen değişkenlerle doğru yapının seçilmesi"dir.
4. **Uygulama verileri küçük.** 4–14 nokta ya da tek tohum. Kuvvet yasası üslerinin belirsizlikleri bootstrap ve LOO ile verildi; ±0,05–0,2 mertebesinde. Bunlar aday denklemler; yeni koşullarda sınanmadan yasa sayılmamalı.
5. **Yeni simülasyon yalnız kendi çözücüyle yapıldı.** Bu çalışmada yeni simülasyon yalnızca §2.1 Gray–Scott ve §3.1-B aktif damla için yapıldı (`qcell.grayscott`, `qcell.active`, L = 64 kutu). Diğer bütün sonuçlar projenin mevcut verisinden.
6. **BSIM4 300 K kartı.** Darbe iyonlaşması ve gövde etkisi bulguları kartın kendi denklemlerine aittir; 4 K fiziği hakkında bir şey söylemez.

## Kaynaklar

- PySINDy: [github.com/dynamicslab/pysindy](https://github.com/dynamicslab/pysindy); de Silva vd., JOSS 5, 2104 (2020); Kaptanoglu vd., JOSS 7, 3994 (2022).
- PDE-FIND: Rudy, Brunton, Proctor, Kutz, Sci. Adv. 3, e1602614 (2017), [arXiv 1609.06401](https://arxiv.org/abs/1609.06401).
- SINDy: Brunton, Proctor, Kutz, PNAS 113, 3932 (2016), [arXiv 1509.03580](https://arxiv.org/abs/1509.03580).
- Zayıf form: Messenger & Bortz, [arXiv 2007.02848](https://arxiv.org/abs/2007.02848); Reinbold vd., [arXiv 1910.08463](https://arxiv.org/abs/1910.08463).
- Gray–Scott: Pearson, Science 261, 189 (1993), [arXiv patt-sol/9304003](https://arxiv.org/abs/patt-sol/9304003).
- Aktif damla: Zwicker vd., Nat. Phys. 13, 408 (2017), [arXiv 1603.01571](https://arxiv.org/abs/1603.01571).
- BSIM4 4.5 kılavuzu (UC Berkeley); SKY130 kartı [skywater-pdk-libs-sky130_fd_pr](https://github.com/google/skywater-pdk-libs-sky130_fd_pr).
- BenDaniel & Duke, Phys. Rev. 152, 683 (1966) (etkin kütle sınır koşulu).
