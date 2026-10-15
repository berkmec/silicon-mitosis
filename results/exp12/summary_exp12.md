# Deney 12: 4 K'de SKY130 nFET — yüzer gövde kilitlenmesi ve akım filamanı (disipatif desen)

Tarih: 2026-10-11. Kod: `experiments/exp12_filament/` (spice_tools, calibrate, single_kink, local_bistability, reduced_model,
spatial, figures), `src/qcell/fet4k.py`, testler `tests/test_exp12.py`. Şekiller: `fig_exp12.png`, `fig_calibration.png`.
Veri: `spice/*.cir|dat|log` (ngspice), `calibration.json`, `single_kink.json`, `local_bistability.json`,
`reduced_bistability.json`, `spatial.json`, `spatial_convergence.json`, `maxwell.json`.

Soru (H22): 4 K'de SKY130 nFET'te darbe iyonlaşması ve donmuş gövde bir disipatif desen, yani akım filamanı üretir mi?
Üretirse bölünür veya çoğalır mı, boyu nedir, kuantumla ilgisi var mı?

## 0. Kısa cevap

1. **Eşik üstünde hayır.** Resmî SKY130 BSIM4 kartıyla gövde 0'dan 1 V'a çıkınca darbe iyonlaşma akımı yalnız 1,3–1,5 kat artıyor.
   Yüzer gövdeyle kink var, ama histerezis ve sağlam bir S-eğrisi yok.
2. **Eşik altında evet: kilitlenme.** Burada geri besleme üstel. BSIM4 gövde katsayısı γ' = 0,26–0,49 ve BSIM4 çoğaltma oranı M ile
   SOI'den bilinen "tek transistör kilitlenmesi" ortaya çıkıyor. Bu, kink ile birlikte gelen bir histerezis: düğme kombinasyonlarının
   %24–97'sinde görülüyor, eşik altı eğimine bağlı.
3. **Filaman çıkıyor ama yalnız belirli koşulda.** Kanal genişliği boyunca ve bir yük direnciyle, donmuş gövdenin potansiyel yayılma
   boyu ℓ kanal genişliğinden çok küçükse akım filamanı oluşuyor (R_□ ≳ 10¹⁰–10¹² Ω/□). V_ds analitik Maxwell geriliminde (0,80 V)
   kenetleniyor; filaman **bölünmüyor, akımla doğrusal genişliyor** (47 → 94 µm).
4. **Kuantumla ilgisi yok.** Filaman onlarca µm genişliğinde, klasik ve sıcak taşıyıcılı. "Küçülerek çoğalma" yok.

## 1. Araçlar ve 4 K varsayımları

- **Gerçek araç.** ngspice 44.2 (WSL Debian) ve resmî SKY130 `nfet_01v8` tt BSIM4 4.5 kartı, 300 K. Kanal akımı `id`, darbe iyonlaşma
  akımı `isub` ve `vth`, cihaz içi büyüklükler olarak çıkarıldı (V_gs 0,3–1,8; V_ds 0–2,6; V_bs 0–1,0 V; W = 1 µm, L = 0,15 µm).
- **BSIM4 4 K'de geçerli değil.** 4 K etkileri açık ve taranan düğmeler olarak verildi:
  - eşik altı eğimi SS (10–80 mV/dek);
  - darbe iyonlaşması çarpanı F_ii (1–100);
  - gövde diyotu V_on = 1,05 V ve eğimi V_s (5–25 mV);
  - donmuş kuyu sızıntısı g_sub (10⁻¹⁵–10⁻¹¹ A/(V·µm));
  - yanal gövde tabaka direnci R_□ (10⁶–10¹² Ω/□);
  - işletme noktası V_ov0.
- **Kompakt model kalibrasyonu (`calibrate.py`).** I_d'de medyan hata %9 (p90 %30). I_sub için β₀,eff = 10,9 V gerekti (p90 %87).
  Uyum orta düzeyde olduğu için ana motor olarak kullanılmadı; yerine BSIM4 tabloları kullanıldı.

## 2. Doğrulama ve bulunan hatalar

| Kontrol | Sonuç |
|---|---|
| ngspice zinciri | cihaz içi `isub` = gövde kaynağı akımı (V_ds ≥ 1,5 V) |
| Uzaysal çözücü, düzgün cihaz | 1 ve 10 hücre, indirgenmiş modelin dal izlemesini birebir veriyor: V_up = 1,56 V, V_down = 0,78 V |
| Filaman kenetlenmesi | analitik Maxwell (eşit alan) V_ds = 0,8018 V; simülasyon 0,800–0,803 V |
| Izgara ve gevşeme | dy = 0,5 / 0,25 / 0,1 µm ve 12 / 48 alt adımda R_□ = 10¹²'de aynı sonuç |
| Testler | 36/36 (3 yeni regresyon testi) |

Yol boyunca bulup düzelttiğim hatalar:
1. BSIM4'ün kendi 300 K eklem diyodu, eklem akımı parametreleri sıfırlanınca bile gövdeyi 0,72 V'ta kenetliyordu. Bu yüzden 4 K diyot
   düğmesini kullanabilmek için tablo yaklaşımına geçtim.
2. Kök sayımına dayalı "bistabilite" ölçütü aşırı genişti. Yerine soğuk cihazdan yukarı ve yüksek daldan aşağı tarama ile dal izleme
   koydum.
3. Örtük adım, pozitif geri besleme bölgesinde ters yöne gidiyordu. Etkin adımı büyüme süresiyle sınırladım.
4. Tek hücrede Laplace köşe hatası vardı.

## 3. Sonuçlar

| Adım | Sonuç |
|---|---|
| Eşik üstü (ngspice, yüzer gövde, V_gs = 1,2 V) | V_b yaklaşık 0,7–0,8 V'a çıkıyor; kink var, histerezis yok (yukarı ve aşağı tarama üst üste) |
| Eşik üstü bistabilite (BSIM4 tablosu, 486 düğme kombinasyonu) | 14'ünde, o da tek bir ızgara adımında (0,05 V) ve aşırı koşullarda (V_s = 5 mV, sızıntı ≈ 0). Sağlam değil |
| Geri besleme eğimi d ln I_ii / dV_b | eşik altında 4,6–9 V⁻¹, eşik üstünde 0,8–1,8 V⁻¹ (300 K, BSIM4) |
| Eşik altı kilitlenme (indirgenmiş model, 108 kombinasyon/SS) | kink + histerezis: SS 10 / 15 / 20 / 30 / 80 mV/dek'te 26 / 39 / 48 / 69 / 105. Medyan pencere 0,04 / 0,04 / 0,04 / 0,28 / 0,16 V; bazı durumlarda cihaz kilitli kalıyor |
| Uzaysal model (bütün satırlar) | **eşik altı** işletme: V_ov0 = −0,1 V (V_gs ≈ V_T − 0,1 V), SS 20 mV/dek, V_s 5 mV, F_ii 10, g_sub 10⁻¹³; W = 100 µm, 200–1000 hücre, her hücrenin kendi gövde potansiyeli, R_L = 5 kΩ, SKY130 Pelgrom eşleşmesi (σ ≈ 12 mV / 0,5 µm hücre) |
| Uzaysal, R_□ = 10¹² Ω/□ (ℓ ≈ 2 µm) | tek filaman. V_ds = 0,80 V'ta kenetli; genişlik V_DD 1,4 → 2,0 V'ta 47 → 94 µm (akımla doğrusal); aşağı taramada daralıp kayboluyor |
| Uzaysal, R_□ = 10¹⁰ (ℓ ≈ 22 µm) | ince ızgarada filaman, kaba ızgarada bütün kanal birlikte geçiyor. **Exp13'te çözüldü:** dy = 0,1 ve 0,05 µm aynı sonucu veriyor (tek filaman); kaba ızgara artefakttı |
| Uzaysal, R_□ ≤ 10⁸ (ℓ ≳ 224 µm > W) | bütün kanal birlikte geçiyor. Yük doğrusunda testere dişi: kararlı bir durağan durum bulunamadı (olası gevşeme salınımı, çözülmedi) |

Filaman ölçütü: ℓ = √(L_gövde / (R_□ · g_sub)), yani donmuş gövde potansiyelinin yanal yayılma boyu. Filaman yalnız W ≫ ℓ olduğunda
oluşuyor. Yük direnciyle (global bağlaşım) tek bir domen kararlı; çoğalma değil genişleme görülüyor. Bu, S-NDC teorisinin (Schöll)
klasik beklentisi.

## 4. Ölçülebilir öngörüler (SKY130, 4 K)

- **Kullanılacak cihaz:** geniş nFET (W = 20–100 µm, L = 0,15 µm), V_gs ≈ V_T − 0,1 V, drenaj seri bir direnç ya da akım kaynağı
  üzerinden sürülüyor.
- **Filaman imzası:** bir akım aralığında V_ds sabit kalıyor (Maxwell gerilimi, modelde ≈ 0,8 V). Bu aralık W ile orantılı olmalı.
- **Gerilim sürüşünde:** histerezis (kink V_up, dönüş V_down).
- **Dar cihazda (W ≲ ℓ):** sabit-V_ds bölümü yok, yalnız ani geçiş var. Çok parmaklı cihazda parmak parmak basamaklar beklenir.
- **Bu ölçüm donmuş kuyunun R_□ ve g_sub değerlerini de verir.** Bu ikisi şu an bilinmiyor ve sonucu belirliyor.

## 5. Geçerlilik sınırları

- BSIM4 300 K'de çıkarıldı. 4 K'ye özgü her şey varsayım ve taranan düğme olarak verildi.
- İndirgenmiş modelde eşik üstü akım kaba (doğrusal büyüme). Kilitlenme analizi eşik altına dayanıyor.
- Modelde olmayanlar: sıcak taşıyıcı ve kafes ısınması, yavaş tuzak/donma kinetiği (literatürdeki dinamik histerezisin bir adayı),
  parazitik bipolar, donmuş kuyuda sığ katkı darbe iyonlaşma kırılımı (Schöll tipi). Si:B için 4 K kırılım parametresi bulunamadı.
- Uzaysal sonuç tek bir temsilci parametre setine dayanıyor. Kanal uzunluğu yönü ve 2B/3B gövde dahil değil. R_□ ≤ 10¹⁰ durumları
  tam yakınsamadı.

## 6. Hipotez durumu

H22 model düzeyinde ve koşullu olarak destekleniyor. Koşullar: eşik altı işletme ve donmuş kuyu R_□ ≳ 10¹⁰–10¹² Ω/□.
- Desen bir hücre bölünmesi değil, genişleyen tek bir akım filamanı.
- Kuantum değil, klasik.
- Ölçülmeden kesinleşmez.

## Kaynaklar

- Kriyojenik CMOS kink ve histerezis: [arXiv 2003.00207](https://arxiv.org/pdf/2003.00207), [arXiv 1811.09987](https://arxiv.org/pdf/1811.09987),
  [arXiv 2104.08467](https://arxiv.org/pdf/2104.08467), [arXiv 2104.03094](https://arxiv.org/pdf/2104.03094), [arXiv 1811.12632](https://arxiv.org/pdf/1811.12632),
  [arXiv 1811.11497](https://arxiv.org/pdf/1811.11497), [arXiv 1908.02656](https://arxiv.org/pdf/1908.02656); FDSOI'de kink yok: [arXiv 1809.09013](https://arxiv.org/pdf/1809.09013)
- S-NDC akım filamanları (Schöll): [Hüpper vd. PRB 47, 15515](https://link.aps.org/doi/10.1103/PhysRevB.47.15515), [Schöll & Drasdo 1990](https://link.springer.com/article/10.1007/BF01309347),
  iki seviyeli model: [Albuquerque vd. 2006](https://www.redalyc.org/pdf/464/46436303.pdf), [arXiv 1409.7283](https://arxiv.org/pdf/1409.7283)
- ESD/snapback filamanları (oda sıcaklığı, karşılaştırma için): [IISc yayınları](https://faculty.dese.iisc.ac.in/mayank/wp-content/uploads/sites/5/2022/12/55.pdf)
- SKY130 model kartı: [skywater-pdk-libs-sky130_fd_pr](https://github.com/google/skywater-pdk-libs-sky130_fd_pr)
