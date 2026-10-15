# Deney 8 — "Bölünme" mekanizmalarının modellenmesi

Üç aday mekanizma ve kendi birleştirmem. Tüm çözücüler analitik sonuçlarla doğrulandı (aşağıda).

## 1. Doğrulama

- Faz alanı çözücüsü: duvar gerilimi σ = 0.9389 (analitik 0.9428, hata %0.41); lineer büyüme hızı σ(q) en çok %0.32 sapma; kütle sürüklenmesi 2.6e-10.
- Bifurkasyon: A_c analitik 0.38490, sayısal 0.38490; eşik civarı genlik ∝ δ^0.502 (süperkritik çatal bifurkasyonu için beklenen 0,5).
- Küme minimizasyonu: LJ7 −16,505384, LJ13 −44,326801, LJ19 −72,659782 (literatür değerleri, hata < 5e-7).

## 2. Mekanizma 2+3 (faz ayrışması + bifurkasyon): Coulomb-frustrasyonlu Cahn-Hilliard, gate ile ekranlı

Düzgün durum, A_c eşiğinde sonlu dalga boylu desene bölünüyor (doğrulanmış süperkritik bifurkasyon). Eşik gate uzaklığına bağlı:

| d | 1.5 | 3.0 | 6.0 | 12.0 | ∞ |
|---|---|---|---|---|---|
| A_c (φ̄=0) | 0.4846 | 0.3983 | 0.3853 | 0.3849 | 0.3849 |

Gate yaklaştıkça (d küçülünce) Coulomb itmesi ekranlanıyor ve A_c **artıyor**: düzgün durum daha geniş A aralığında kararsız, yani desen/bölünme daha geniş aralıkta sürüyor. Ölçülen dalga boyu d=1,5'te d=∞'a göre 1.58/1.73/1.53 kat büyük (φ̄=0, A=0,1/0,2/0,3); d ≳ 12'de gate'in etkisi yok.
**Gate'in bölünmeyi başlatması:** φ̄ = −0,4, A = 0,15'te düzgün durum d ≥ 6'da kararlı (σ_max ≈ 0, desen yok), d ≤ 3'te kararsız (σ_max > 0, desen oluşuyor). Yani aynı yük/frustrasyon gücünde gate'i yaklaştırmak, bölünmeyi (düzgün durumun kararsızlaşıp hücrelere ayrılmasını) başlatıyor.

| φ̄ | A | d | λ_lineer (en hızlı büyüyen) | λ_ölçülen (t=3000) | σ_max | desen? |
|---|---|---|---|---|---|---|
| -0.4 | 0.05 | 1.5 | 13.1 | 23.9 | 0.048 | evet |
| -0.4 | 0.05 | 3.0 | 13.1 | 21.3 | 0.044 | evet |
| -0.4 | 0.05 | 6.0 | 13.0 | 18.3 | 0.043 | evet |
| -0.4 | 0.05 | 12.0 | 13.0 | 17.3 | 0.043 | evet |
| -0.4 | 0.05 | inf | 13.0 | 17.2 | 0.043 | evet |
| -0.4 | 0.15 | 1.5 | 15.7 | 21.3 | 0.016 | evet |
| -0.4 | 0.15 | 3.0 | 16.2 | 18.0 | 0.003 | evet |
| -0.4 | 0.15 | 6.0 | 6283.2 | 15.5 | -0.000 | hayır |
| -0.4 | 0.15 | 12.0 | 6283.2 | 15.3 | -0.000 | hayır |
| -0.4 | 0.15 | inf | 6283.2 | 15.3 | -0.000 | hayır |
| 0.0 | 0.1 | 1.5 | 9.3 | 26.0 | 0.189 | evet |
| 0.0 | 0.1 | 3.0 | 9.3 | 21.4 | 0.182 | evet |
| 0.0 | 0.1 | 6.0 | 9.2 | 17.6 | 0.181 | evet |
| 0.0 | 0.1 | 12.0 | 9.2 | 16.5 | 0.181 | evet |
| 0.0 | 0.1 | inf | 9.2 | 16.5 | 0.181 | evet |
| 0.0 | 0.2 | 1.5 | 9.8 | 22.9 | 0.133 | evet |
| 0.0 | 0.2 | 3.0 | 9.7 | 14.9 | 0.117 | evet |
| 0.0 | 0.2 | 6.0 | 9.7 | 13.4 | 0.114 | evet |
| 0.0 | 0.2 | 12.0 | 9.7 | 13.3 | 0.114 | evet |
| 0.0 | 0.2 | inf | 9.7 | 13.2 | 0.114 | evet |
| 0.0 | 0.3 | 1.5 | 10.5 | 18.2 | 0.080 | evet |
| 0.0 | 0.3 | 3.0 | 10.4 | 12.6 | 0.055 | evet |
| 0.0 | 0.3 | 6.0 | 10.2 | 11.9 | 0.051 | evet |
| 0.0 | 0.3 | 12.0 | 10.2 | 11.9 | 0.051 | evet |
| 0.0 | 0.3 | inf | 10.2 | 11.9 | 0.051 | evet |

NEGATİF SONUÇ (kayıt: `division_diag.json`): hazır tek bir domain Coulomb frustrasyonunda **bölünmedi**. t=3000'de hiçbir başlangıç ve A için n=1:

| başlangıç | A | n (t=3000) | eşdeğer yarıçap | pozitif faz oranı |
|---|---|---|---|---|
| disk R0=18 | 0.0 | 1 | 17.4 | 0.095 |
| disk R0=18 | 0.004 | 1 | 16.3 | 0.084 |
| disk R0=18 | 0.016 | 1 | 12.6 | 0.050 |
| disk R0=18 | 0.05 | 1 | 3.9 | 0.005 |
| rod 60x12 | 0.0 | 1 | 20.0 | 0.126 |
| rod 60x12 | 0.004 | 1 | 19.1 | 0.115 |
| rod 60x12 | 0.01 | 1 | 17.7 | 0.098 |
| rod 60x12 | 0.02 | 1 | 15.3 | 0.073 |
| ellipse 30x9 | 0.0 | 1 | 15.7 | 0.078 |
| ellipse 30x9 | 0.004 | 1 | 14.6 | 0.067 |
| ellipse 30x9 | 0.01 | 1 | 12.7 | 0.050 |
| ellipse 30x9 | 0.02 | 1 | 9.2 | 0.027 |

Disk A ≲ 0,02'de yavaşça küçülüyor, A = 0,05'te çözünüyor (φ̄ = −0,8 homojen durumu kararlı). Çubuk (60×12) A=0'da bile daireye dönüyor: 2B'de Rayleigh-Plateau kopması yok. Yani bu modelde 'bölünme' tek hücrenin şekil kararsızlığı değil, düzgün durumun bifurkasyonla çok sayıda hücreye ayrılmasıdır.

## 3. Mekanizma 1 (elektrostatik/Coulomb): yüklü parçacık kümesinin fisyonu

Kısa menzilli çekim (Lennard-Jones) + gate ile ekranlı Coulomb itmesi; tek küme enerjisi E(N) ile iki parça enerjisi E(N₁)+E(N₂) karşılaştırıldı (parçalar sonsuz uzakta).

- d=inf: g_c(N) = N=10: 0.827, N=13: 0.570, N=16: 0.570, N=19: 0.393, N=22: 0.393, N=26: 0.271; log-log eğim = -1.08
- d=4.0: g_c(N) = N=10: 0.827, N=13: 0.827, N=16: 0.570, N=19: 0.570, N=22: 0.570, N=26: 0.393; log-log eğim = -0.73
- d=2.0: g_c(N) = N=10: 1.200, N=13: 1.200, N=16: 0.827, N=19: 0.827, N=22: 0.570, N=26: 0.570; log-log eğim = -0.90

Rayleigh/sıvı-damla ölçeklemesi için g_c ∝ 1/N beklenir (eğim −1); gate ekranlaması eşiği büyütmeli (daha fazla yük gerekir).

## 4. Sınırlar

- Faz alanı modeli jenerik: φ'nin fiziksel karşılığı (elektron/boşluk birikintisi, ferroelektrik polarizasyon) bu çalışmada bağlanmadı; parametreler birimsiz.
- Tek hazır domainin bölünmemesi 2B ve bu parametre aralığı için geçerli; daha büyük sistem, gürültülü başlangıç veya 3B'de farklı olabilir.
- Küme fisyonu klasik (kuantum etkisi yok); basin-hopping kısa koşulu, enerjiler üst sınırdır, eşikler biraz kayabilir.
- g ızgarası geometrik (×1,45 aralık): g_c değerleri ±1 adım (~%45) belirsiz; eğimler bu ızgara çözünürlüğüyle okunmalı (d=∞ için −1,08, Rayleigh beklentisi −1).
- Yüksek g'de (örn. N=19, d=∞, g=1,2) E_tek − E_parça kazancı düşer: bu fiziksel değil, kısıtsız tek-küme aramasının kendisinin parçalı bir duruma inmesidir. Eşikler ilk geçişten okunduğu için etkilenmez; yüksek-g kazanç eğrisi nicel okunmamalı.
- Gate uzaklığı d birimsiz (kapasitör mesafesi ölçeği); SKY130'a çeviri bu çalışmada yapılmadı.