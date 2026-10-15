# Aşama 4 test planı: SKY130 minimum aygıtlarında kazara nokta / Coulomb blokajı araması

> **Geçerlilik uyarısı (exp9, 2026-10-10):** Bu planın 4 K öngörüsü (%4–21 kazara nokta) iki nedenle güvenilir değil.
> (1) Görüntü modelinin düzeltmesi z<10 nm'de ekstrapole; bu kısım kesilince katı oran %10 → %6,25 (efektif), %3,75 → %2,5 (nominal).
> (2) Modelde iyonlaşma T'den bağımsız; NA = 1,7e17 cm⁻³ ise Mott eşiğinin (~3,7e18, bellekten) çok altında, 4 K'de tam iyonlaşma varsayımı geçersiz.
> Ayrıntı: `results/exp9/COMPARISON.md`. Bu plan, dondurma (freeze-out) bilinçli bir model çıkana kadar karar kapısı olarak kullanılmamalı.
> 300 K negatif sonucu (ölçek ~30× kısa) ise bu sorunlardan etkilenmez.

Bu plan, Aşama 2b tahminini (H10) çürütülebilir hâle getirir. Kaynak hesaplar: `results/exp5..exp7`. Hiçbir üretim veya ölçüm yapılmadı;
donanım ve bütçe kararı kullanıcıdadır. Aşağıdaki tüm sayılar **model tahminidir**, ölçüm değildir.

## 1. Test edilen tahmin (H10)
126 nm × 420 nm etkin kanallı SKY130 `nfet_01v8` minimum aygıtların yaklaşık **%4–21'inde** (katı ölçüt %4–10), eşik civarında gate taramasında, 4 K'de
Coulomb blokajı salınımı görülür (nokta kenarı ≈ 15–23 nm, E_C ≈ 2–4 meV, yani E_C/k_B ≈ 25–40 K).

## 2. Öngörülen, ölçülebilir imzalar
| İmza | Tahmin | Nasıl ölçülür |
|---|---|---|
| Tepe aralığı ΔV_g | e/C_g; 20 nm'lik nokta için C_g ≈ 3,3 aF (paralel plaka, 4,148 nm oksit) → ΔV_g ≈ 50 mV (30–100 mV aralığı) | I_d–V_g, V_sd ≈ 20–100 µV |
| Elmas yüksekliği | E_C/e ≈ 2–4 mV | I_d–(V_g, V_sd) haritası |
| Sıcaklık | E_C/kT = 8,5 (4 K) → görünürlük 0,99; 10 K'de 0,40; 20 K'de 0,03 (E_C = 34 K için) | 1,5–20 K taraması |
| Akım seviyesi | tünel bariyeri 20–30 meV, genişlik 5–20 nm: R_T ≈ 0,6 MΩ – 140 GΩ; 100 µV'ta 0,4 fA – 80 pA | fA çözünürlüklü elektrometre gerekir (R_T ≳ 1 GΩ için) |
| Tek cihaz ölçekleme | pik yüksekliği ∝ 1/(R_s + R_d), pik genişliği ∝ kT | sıcaklık bağımlılığı |

## 3. Kaç aygıt ölçülmeli (binom güç analizi, `results/exp7/testplan.json`)
≥ 3 pozitif aygıtı %95 olasılıkla görmek için gerekli aygıt sayısı:

| gerçek oran p | gerekli n |
|---|---|
| %12 | 51 |
| %8 | 77 |
| %4 | 156 |
| %2 | 313 |

n = 192 aygıt: p ≥ %4 ise ≥ 1 pozitif görme olasılığı 0,9996, ≥ 3 pozitif olasılığı 0,984. Hiç pozitif çıkmazsa (0 / 192) %95 güvenle p < ≈ %1,6.

## 4. Karar ölçütleri
- **GO:** ≥ 3 aygıtta (a) periyodik tepeler (ΔV_g sabit ±%20), (b) çok-kapılı veya V_sd haritasında elmas yapısı, (c) tepe genişliği ve yüksekliğinin
  kT ile tutarlı sıcaklık bağımlılığı. Üçü birden gerekir.
- **NO-GO (130 nm hattı kuantum için):** ≥ 150 aygıtta hiçbir elmas yapısı yok. Bu durumda SKY130 yalnızca klasik kontrol/okuma elektroniği olarak kalır.
- **Belirsiz:** Yalnızca düzensiz tepeler veya RTN benzeri basamaklar var, elmas yok: kusur kaynaklı rastgele telegraf gürültüsü ile ayırt edilemiyor, ölçüm
  koşullarını (sıcaklık, filtre, V_sd) değiştirerek tekrar edilir.

## 5. Yanlış pozitif riskleri
- **RTN / tuzak kaynaklı basamaklar:** ΔV_g periyodikliği ve elmas yapısı olmadan "nokta" sayılmaz.
- **Ölçüm devresi:** fA akımlarda kablo kaçakları ve 1/f gürültü; ölçüm hattı doğrulaması için bilinen bir direnç/diyot referansı kullan.
- **Tahmin kusurları:** tünel direnci modellenmedi (üst sınır olabilir), birikinti istatistiği 80 örneğe dayanıyor, E_C sabiti MOSFET geometrisine tam uymuyor.

## 6. Yapılması gerekenler (kullanıcı kararı + doğrulanmamış bağımlılıklar)
1. Üretim yolu: açık MPW servislerinin güncel erişilebilirliği ve fiyatı **doğrulanmadı**.
2. Yerleşim: resmî PDK katman haritası ve DRC kuralları ile KLayout/Magic doğrulaması gerekir; bu oturumda **yerleşim (GDS) üretilmedi**,
   doğrulamasız bir GDS üretmek hatalı çip riski taşır.
3. Ölçüm: ≤ 4 K kriyostat, fA elektrometre, çok kanallı anahtarlama; yoksa harici laboratuvar ortaklığı.
4. Kriyojenik model verisi: NIST SKY130 veri seti dosyaları (Aşama 3) sunucu hatasıyla indirilemedi; yeniden denenecek.
