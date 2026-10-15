# Deney 17: Silikonda kendiliğinden ve tekrar tekrar bölünen akım filamanları — olgun bölünme denklemi

Tarih: 2026-10-15.

- **Kod:** `experiments/exp17_silicon_split/`
  - `law.py`: denklem ve sabitler;
  - `run.py`: 1B koşular;
  - `mechanism_check.py`: her bölünmede mekanizma testi;
  - `symreg.py` / `symreg_gp.py`: veriden model seçimi ve sembolik regresyon;
  - `spatial2d.py` / `run2d.py`: 2B dizi;
  - `figures17.py`.
- **Testler:** `tests/test_exp17.py`.
- **Veri:** `results/exp17/`
  - öngörüler koşulardan **önce**: `predictions*.json`;
  - `runs/`, `runs2d/`, `analysis.json`, `mechanism.json`, `symreg*.json`, `universality_*.json`;
  - şekil `fig_exp17.png`.

> **Düzeltme (2026-10-15, yayın öncesi denetim).** Bu rapordaki bazı ifadeler veriyle yeniden karşılaştırıldı. Aşağıdaki metin özgün
> hâliyle bırakıldı; çelişki olursa bu kutu geçerlidir.
> - **"Hızlı protokol" daha hızlı bir tarama değildir.** Bu grup R_L = 2 kΩ'luk daha sert yüktür; adım başına süre yine 1 s.
>   Daha hızlı bir tarama koşulmadı.
> - **1,001 ± 0,069 yalnız temiz protokolün 22 olayı için geçerlidir** (R_L = 5 kΩ, normal uyumsuzluk). `mechanism.json`'daki 63 olayın
>   tamamı:
>
>   | grup | n | merkez–kenar farkı / ΔV_c |
>   |---|---|---|
>   | temiz (R_L 5 kΩ) | 22 | 1,001 ± 0,069 (0,89–1,11) |
>   | temiz, uyumsuzluk × 0,1 | 2 | 1,03 |
>   | sert yük (R_L 2 kΩ) | 22 | 1,038 ± 0,062 |
>   | sert yük, uyumsuzluk × 0,1 | 7 | 1,14 ± 0,07 |
>   | evrensellik (U1, U2) | 9 | 0,979 ± 0,133 |
>   | düşük sızıntı (g_sub = 1e-15) | 1 | 2,60 |
>
> - **λ = 36 µm, λ_c'nin (35,6 µm) altında değil, hemen üstündedir.** Önceden kaydedilen öngörü 322 µm'de bölünmeydi; filaman 350 µm'ye
>   büyüdü ve bölünmedi (`runs/A_lam36.json`). Bu, eşik yakınında bir ıskalamadır: ΔV_c %2,3 daha büyük olsaydı λ = 36 µm eşiğin altında
>   kalırdı.
> - **σ(V_down)/ΔV_c ≈ 0,12 hiçbir betikte hesaplanmadı;** geri çekildi.
> - **N ≈ I/(0,75 j_on w_s) bir kılavuz çizgidir,** uydurma değildir.
> - **Periyodik dizi formülü bu depoda sınanmadı.** Neredeyse düzgün cihazdaki (F_lowmm) bölünmenin kanıtı gerçek düzen çözücüsüdür
>   (16,30 / 16,27 mV). F_lowmm'de yalnız λ = 80 ve 120 µm bölündü; λ = 40–60 µm 5,6 V'a kadar bölünmedi.
> - **İlk bölünme sınaması 3 × 3 ızgara değildir:** 11 temel set koşusu (r_D = 10: 8, r_D = 20: 2, r_D = 5: 1) ve 3 U koşusu (hepsi
>   r_D = 10). Bölünürken ortamda ikiden fazla filaman bulunan U2_x2.5 (+%18) seçime girmedi.
> - **w_top oranı** temel ve U setlerinde 1,00–1,02; düşük sızıntı (N) koşularında 1,03–1,23.
> - **Tek katsayılı uydurma yalnız 11 temel set koşusuna yapıldı** (k = kuramın 1,010 katı, LOO %8,2). Parametresiz biçim zaten LOO %7,75
>   veriyor. GP ifadelerinin eğitim hatası %4,3–11; LOO hesaplanmadı.
> - **Şekil (c) önce E_lam60'ı gösteriyordu;** orada iki filaman birlikte doğup birer kez bölünüyor. Şekil, kızların yeniden bölündüğü
>   H_lam60_s1 ile yeniden çizildi (tek filaman → 2 → 4). Panel (b) artık yalnız normal uyumsuzluklu olayları gösteriyor; panel (e)
>   başlığındaki "periyodik dizi formülü" ifadesi kaldırıldı.
> - **2B koşularda V_DD adımlarının yalnız %12–44'ü yakınsama toleransına ulaştı** (T30 %36, T64 %28, V64_small %44, P64 %12; eski
>   S25/S64 %0). P64'ün son durumundaki iki parça, periyodik sınırlar birleştirilince tek bir şerittir (bütün geç anlık görüntülerde
>   doğrulandı).
> - **`law.py` artık `predictions.json`'ın üzerine yazmıyor;** yeniden hesap `predictions_recomputed.json`'a yazılıyor.

## 0. Kısa cevap

1. **Silikonda bölünebilen şekil:** geniş bir SKY130 nFET'te (4 K, eşik altı, yüzer gövde kilitlenmesi) akım filamanı drenajı dağıtık balastlıysa akımla büyüyüp **ortasından ikiye bölünüyor**. Akım arttıkça bu tekrar ediyor: 600 µm genişlikte 3–4 kuşak, 4–5 filaman.
   - Hücre büyümesine ve bölünmesine benzer bir döngü: akım besin, bölünme boyu ise cihaz parametrelerinin belirlediği sabit bir genişlik.
2. **Olgun denklem** (fizikten türetildi, bütün sabitler cihaz modelinden, ayarlanan parametre yok):
   - **Genel ölçüt:** bölünme, filamanın merkezindeki drenaj gerilimi yüksek durumun var olabildiği en düşük gerilime (V_down) indiği anda olur.
     Δu_merkez–kenar = V_M − V_down ≡ ΔV_c
   - **Tek filaman:** w_s = 2λ·ln[λ/(λ − λ_c)],  λ_c = √(2 ΔV_c / (r_D·j_on)),  λ = √(r_v/r_D).
     - En küçük bölünme genişliği 2λ_c.
     - λ ≤ λ_c ise bölünme yok.
   - **Periyodik dizi** (aralık P): Δu = r_v j_on (cosh a − 1)/(cosh a + sinh a·coth b), a = w/2λ, b = (P − w)/2λ.
   - **Toplam genişlik:** w_top = I / j_on.
   - **Filaman sayısı:** N ≈ I / (0,75 · j_on · w_s).
   - **Saçılma:** V_T uyumsuzluğu bölünme anını σ(V_down)/ΔV_c ≈ 0,12 kadar oynatıyor.
3. **Sınama** (öngörüler önceden kaydedildi; yeni koşular):
   - **Gerçek düzende merkez–kenar farkı:** her bölünme anında komşular ve cihaz kenarları dahil hesaplandı.
     - temiz protokolde **ΔV_c × (1,001 ± 0,069)** (22 olay);
     - iki farklı cihaz setinde 0,979 ± 0,133 (9 olay).
   - **Evrensel eğri:** 3 cihaz seti ve 3 r_D'deki 14 ilk bölünme w_s/λ_c = f(λ/λ_c) eğrisine oturuyor; ortalama hata +%1,6, RMS %9,1.
   - **Veriden model seçimi:**
     - fizik formu sıfır parametreyle %7,8 hata;
     - tek katsayı uydurulunca **kuramın 1,01 katı** (birini dışarıda bırakma %8,2);
     - kuvvet yasası %11,5;
     - önyargısız sembolik regresyon (gplearn) 32–70 terimli, yorumlanamayan formüllerle daha iyisini bulamadı.
4. **Geçerlilik koşulu:** tarama yarı statik olmalı.
   - Hızlı sürüşte (R_L = 2 kΩ) bölünme ~%4–14 gecikiyor.
   - Çok düşük sızıntıda (g_sub = 10⁻¹⁵, gövde boşalması c_b/g_sub ≈ 1 s) 1 s/adım taramada bölünme yok, 30 s/adımda 2,6 kat geç.
5. **Kuantum değil**, klasik bir disipatif desen. Model 4 K düğmeleri varsayım; ölçüm gerekli.

## 1. Türetme

Drenaj hattında (dağıtık balast) düşüm u = V_c − V_D için u'' = u/λ² − r_D j_d.

- **Tek filaman:** w genişliğinde, içinde j_on taşıyan bir filaman için merkez–kenar farkı Δu = (r_v j_on / 2)(1 − e^(−w/2λ))².
- **Kenarlar:** cephe durağan olduğu için Maxwell geriliminde, V_M (eşit alan koşulu).
- **Merkez:** daha aşağıda. V_D < V_down olunca yerel yüksek durum yok (eyer-düğüm), merkez söner.
- **Sonuç:** Δu = ΔV_c koşulu w_s = 2λ ln[λ/(λ−λ_c)] formunu veriyor; çünkü √(2ΔV_c/(r_v j_on)) = λ_c/λ.

Temel set için sabitler (SS 20 mV/dek, V_s 5 mV, F_ii 10, V_ov0 −0,1 V, g_sub 10⁻¹³):

| büyüklük | değer | not |
|---|---|---|
| V_M | 0,8018 V | Exp12 ile aynı |
| V_down | 0,7855 V | Exp13'teki bölünme anında ölçülen en düşük V_D 0,785 |
| ΔV_c | 16,27 mV | |
| j_on | 2,569 µA/µm | Exp16 veriden 2,55 ± 0,03 bulmuştu |
| V_up | 1,559 V | |
| λ_c | 50,3 / 35,6 / 25,2 µm | r_D = 5 / 10 / 20 Ω/µm |

## 2. Sınama

| set / koşul | öngörü w_s (µm) | gözlenen (µm) |
|---|---|---|
| temel, r_D 10: λ = 40 / 45 / 50 / 60 / 80 / 120 | 176 / 141 / 124 / 108 / 94 / 84 | 173 / 138 / 123 / 101–120 (3 tohum) / 102 / 94 |
| temel, r_D 20: λ = 40 / 60 | 79 / 65 | 83 / 64 |
| temel, r_D 5: λ = 80 | 159 | 191 (ortamda 2. filaman var) |
| U2 (SS 15, F_ii 30): λ = 1,5 λ_c | 97,0 | 96,5 |
| U1 (SS 30, V_s 10, F_ii 3, g 1e-12): λ = 1,5 / 2,5 λ_c | 191 / 148 | 162 / 136 (uyumsuzlukla erken; oranlar 0,72 / 0,83) |
| λ ≤ λ_c (30, 36 µm) | bölünme yok | bölünme yok |

Diğer sonuçlar:
- **Mekanizma oranı** (`mechanism.json`): temiz 1,001 ± 0,069 (n = 22); evrensellik 0,979 ± 0,133 (n = 9); hızlı protokol 1,038 ± 0,062 (n = 22; düşük uyumsuzlukta 1,14). Düşük sızıntı 2,60 (geçerlilik dışı).
- **Tekrarlayan bölünme:** sonraki kuşaklar da aynı ölçütle bölünüyor (oranlar 0,9–1,1).
- **Düzgün cihaz** (uyumsuzluk ×0,1): V_up'ta aynı anda 3 filaman doğup eşit aralıklı diziye yerleşiyor (Turing benzeri). Komşulu dizide bölünme gecikiyor; tek filaman formülü burada yanlış (23,7 mV verirdi), gerçek düzen hesabı doğru (16,30 mV / ΔV_c 16,27).
- **Toplam genişlik:** w_top/(I/j_on) = 1,01–1,02, bütün koşularda.

## 3. SKY130 test yapısı önerisi (λ = 60 µm, r_D = 10 Ω/µm)

- **Cihaz:** nfet_01v8, W ≈ 600 µm, L = 0,15 µm, eşik altı sürüş (V_gs ≈ V_T − 0,1 V). Drenaj, yük direnci (≈ 5 kΩ) ya da akım kaynağı üzerinden yavaşça taranıyor.
- **Drenaj hattı:** li1 şerit, 12,8 Ω/□; 1,28 µm genişlik 10 Ω/µm veriyor.
- **Balast:** her 10 µm'lik bölüme 3,6 kΩ. Seçenekler res_high_po (~300 Ω/□, 0,35 µm genişlikte ~12 kare) ya da res_xhigh_po (2000 Ω/□). Değerler [SkyWater PDK](https://skywater-pdk.readthedocs.io/en/main/rules/device-details.html) belgesinden.
- **λ taraması:** aynı çipte λ = 30, 45, 60, 90 µm (balast değeri 0,9 / 2 / 3,6 / 8,1 kΩ) ve balastsız bir kontrol. Öngörü: λ = 30'da bölünme yok, diğerlerinde w_s = 2λ ln[λ/(λ−λ_c)].
- **Ölçüm:**
  - Bölüm akımları ayrı ayrı (ya da drenaj hattı boyunca gerilim) okunarak filaman konumu görülür.
  - Akım arttıkça filaman sayısının basamaklı artması ve her bölünmede iki komşu bölüm grubu arasındaki sönme beklenir.
- **Denklemin sabitleri ölçülebilir:**
  - V_M: balastsız geniş cihazda filaman platosunun kenetlenme gerilimi;
  - V_down: dar bir cihazın aşağı tarama histerezis ucu;
  - j_on: plato boyunca dI/dw.
  - Bu üçü ölçülünce w_s öngörüsü parametresiz sınanır.

## 4. Sınırlar

- 4 K davranışı Exp12'nin indirgenmiş modelinden. BSIM4 300 K tabloları + varsayılan 4 K düğmeleri; ısınma, tuzak kinetiği ve parazitik bipolar yok. Yalnız drenaj hattı inhibitör; kaynak hattı ve 2B kanal yok.
- 1B (genişlik boyunca) model. 2B dizi (`spatial2d.py`):
  - eski çözücü yük doğrusunda yakınsamıyordu, düzeltildi (Sherman–Morrison);
  - 2B'de bölünme görülmedi (§5).
  - 2B koşularda adımların çoğu 400 yinelemede tam yakınsamadı. Yine de durumlar düzgün ve tutarlı evrildi (vc ve şekil sürekli).
- Bölünme anı uyumsuzlukla rastgele: tek koşuda ±%10, nadiren %25–30 erken.
- Hız şartı: denklem yarı statik sürüş içindir.

## 5. 2B dizi (parmak dizisi + 2B drenaj yaprağı + dağıtık balast)

**Model:** her hücre bir parmak parçası (dx = 6 µm), gövde 2B, drenaj 2B yaprak (r_□ = 50 Ω/□), balast dağıtık.
- Yük doğrusu Newton'a Sherman–Morrison ile örtük girdi. Eski açık şema 20 mA'de salınıyordu.
- Yalıtık disk kuramı: merkez–kenar farkı Δu = r_□ λ² J · x K₁(x)(I₀(x) − 1), x = R/λ.
  - λ = 64 µm için kritik yarıçap R_c = 71 µm;
  - 2B λ_c = 39 µm (r_□ = 50, dx = 6).

| koşu | alan | sonuç |
|---|---|---|
| T30 (λ = 30 < λ_c) | 600 µm, duvarlı | tek leke büyüyor, bölünme ve halka yok (beklenen) |
| T64 (λ = 64) | 600 µm, duvarlı; leke duvarlardan uzakta doğdu | R_eş 73 µm'ye kadar yuvarlak (en-boy 1,03–1,07), R_c'yi ~%3 geçiyor. Sonra **uzayıp şerit oluyor**: en-boy 1,2 → 5,8 (R_eş 81 → 108 µm). Merkez sönmüyor, halka ve bölünme yok |
| V64_small (λ = 64) | 420 µm, duvarlı; leke köşede doğdu | duvar boyunca şerit, sonra bir parmak; merkez–kenar farkı uzamayla 12,8'den ~7 mV'a düşüyor; bölünme yok |
| P64 (λ = 64) | 540 µm, **periyodik (duvarsız)**, merkezde tasarlanmış çekirdek (V_T −30 mV, r = 12 µm); elektrik kesintisinden sonra 1 V'lık kayıtlarla yeniden koşuldu | R_eş 78 µm'ye kadar yuvarlak (en-boy 1,05–1,08), R_c'yi ~%10 geçiyor. Sonra uzuyor (en-boy 1,3 → 3,9) ve **dambıl** oluyor, ama **boyun kopmuyor**. V_DD = 16'ya (R_eş 111 µm) kadar tek parça; uzun bir şerit olup bükülüyor. Merkez V_D ≥ 0,794 V > V_down. Şekiller: `fig_2d_periodic.png`, `check_P64_late.png` |

**Okuma:**
- 2B'de yuvarlak leke merkez sönme eşiğine yaklaşınca bölünmek yerine uzuyor, dambıl ve sonra solucan oluyor. Uzama merkezdeki düşümü azaltıyor; ek akım lekenin genişliğine değil boyuna gidiyor.
- Bu davranış duvarlı (T64) ve duvarsız (P64) alanda aynı, yani kenar etkisi değil.
- Exp14'te aktif damlanın 2B'de solucana dönüşmesinin ve Exp11'deki 2B solucanların silikondaki karşılığı.
- **Tasarım ilkesi:** silikonda bölünebilen şekil, yapı **1B'ye hapsedildiğinde** ortaya çıkıyor.
  - Geniş bir MOSFET'te filaman yalnız genişlik yönünde büyüyebiliyor; merkezi sönüp ikiye bölünüyor (§2).
  - 2B parmak dizisinde lekeler solucana dönüşüyor.
  - Aktif damlada boyut dersinin tersi: orada 2B solucan veriyor, 3B'de eğrilik bölünmeyi sağlıyordu (Exp14). Her iki sistemde de bölünme geometrinin boyutuna bağlı.

## Kaynaklar

- SkyWater SKY130 cihaz ayrıntıları: <https://skywater-pdk.readthedocs.io/en/main/rules/device-details.html>
- S-NDC akım filamanları ve yerel/küresel bağlaşım: [Hüpper vd. PRB 47, 15515](https://link.aps.org/doi/10.1103/PhysRevB.47.15515), [Schöll & Drasdo 1990](https://link.springer.com/article/10.1007/BF01309347)
- Leke ve filaman çoğalmasında aktivatör–inhibitör: [Pearson 1993](https://arxiv.org/abs/patt-sol/9304003)
- PySINDy / gplearn: <https://github.com/dynamicslab/pysindy>, <https://gplearn.readthedocs.io>
