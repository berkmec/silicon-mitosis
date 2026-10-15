# İnceleme: Exp6 (Aşama 2b), Exp10 ve Exp11 — tutarlılık, denge ve mantık denetimi

Tarih: 2026-10-11. İstek: "son 2 raporu karşılaştır, dengesizlik veya mantık hatası var mı; Exp2b Monte Carlo'yu da ekle".

Yöntem:
1. Kendi kontrolüm ve yeni doğrulama hesapları:
   - 3B Δν taraması (−0,004 … +0,004),
   - A için şarj şiddeti taraması,
   - A'nın dambıl koşusunun 3 kat uzatılması,
   - Gray–Scott soy takibi,
   - 2B tabaka güç/boy tahmini,
   - Exp6'nın kutu kenarı (Campbell) hesabı ve düzeltilmiş yeniden koşusu.
2. Üç bağımsız denetçi. Bağlamım olmadan yalnız dosyaları okudular ve dosya değiştirmediler:
   - (i) Exp10 ↔ Exp11 fiziği,
   - (ii) Exp11'in iç mantığı ve dengesi,
   - (iii) Aşama 2b Monte Carlo.

   Kritik bulguların hepsini kendim yeniden hesapladım veya yeniden koştum. "Doğrulandı" sütunu bunu gösterir.

Şiddet: **kritik** = ana sonucu değiştirir, **orta** = bir iddiayı daraltır, **düşük** = sayı veya ifade düzeltmesi.

## 0. Kısa sonuç: ne değişti?

| Rapor | Önceki iddia | Düzeltilmiş hâli |
|---|---|---|
| Exp11 | "Bölünmede B (disipatif) açık ara güçlü" | **Adil değil.** Aynı 2B çekirdekte büyüme→ikili bölünme A 0/28, B 0/30; kalıcı bölünme 0/28'e karşı 2/30 (Fisher p = 0,49). 3B'de bölünme Δν = 0 dahil (sızıntılı şarj ≡ Ohta–Kawasaki ≡ ekransız 3B Coulomb) Δν'nin her değerinde sürekli. Yani disipasyona özgü değil. Belirleyici etkenler: boyut (3B kopar, 2B kopmaz), uzun menzilli itme (kapı onu ekranlar) ve otokatalitik aktivatör–inhibitör yapısı (Gray–Scott 2B'de de çoğalır) |
| Exp10 | "Yakın kapı (süperkapasitör) bölünmeyi açıyor" | **Model içi tutarsızlıktan doğan yapay sonuç.** Kapı Hartree ve değişimi ekranlıyor, korelasyon vekilini ekranlamıyor. Tutarlı (kapı ekranlı) HF'de tabaka her kapıda kararlı: S/K_kin ≥ +1,006 (C 0,05–3000 µF/cm², n 10⁹–10¹³ cm⁻²) |
| Aşama 2b (Exp6) | "Minimum aygıtların %4–21'inde a ≈ 15–23 nm kazara nokta, E_C/k_B 25–40 K" | **Geri çekildi.** Bu noktaların hepsi kutu kenarı artefaktı: dört kenar dışlanınca 0/80. Yalnız rastgele RDF potansiyeli (kubbe çıkarılmış) ele alınınca çukurlar sık ama büyük (medyan a = 35–41 nm, E_C/kT 4 K'de ≈ 3,5–4). Ekranlama olmadığı için bu sayılar da güvenilmez |
| Kuantum hedefi | "R ≲ 18–20 nm (izole küre)" | Geometriye bağlı: Si içinde FEM tabanlı R ≈ 6–9 nm, izole küre 18 nm, kapılı disk (yalnız plaka) ≲ 30 nm |
| B'nin silikon karşılığı | "en iyimser ≈ 66 nm" | Kink literatüründeki mekanizma bant-bant darbe iyonlaşması (E ≈ E_g), 45 meV değil. En iyimser durumda 2B tabakada ≈ 123 nm, 3B'de ≈ 1,6 µm. Sonuç (kuantum yok) aynı, sayılar düzeldi |

## 1. Bulgular

### 1.1 Exp11 (süperkapasitör A ile disipatif B)

| # | Şiddet | Bulgu | Kanıt | Doğrulandı |
|---|---|---|---|---|
| 1 | kritik | "B açık ara güçlü" aynı koşullarda tutmuyor. A'nın sürüş şiddeti (J) taranmamıştı. A'nın 3B karşılığı yoktu. §4'teki "hız, gürültü, rampa, eşik" satırları test edilmemiş bir A'ya karşı yazılmıştı | 2B ortak çekirdek: büyüme→ikili 0/28 ile 0/30; kalıcı 0/28 ile 2/30 | ben + denetçi (ii) |
| 2 | kritik | 3B bölünme disipasyona özgü değil | 3B, ε = 0,20: Δν = +0,004 / +0,002 / 0 / −0,001 / −0,002 hepsi t ≈ 2000'de iki eşit kıza bölündü (R_kız 6,1 / 7,2 / 8,8 / 9,9 / 11,1). Δν = −0,004 t = 3000'de tek ve yuvarlak (R = 19,8). Bölünme yarıçapı Δν boyunca sürekli; Δν = 0'da özel bir şey yok. Sıvı damla (fisyon) tahmini R ≈ 6–7 w, gözlenenle aynı mertebede | ben (koşular) + denetçi (ii) |
| 3 | kritik | 2B B sınıflandırması yanlıştı. Kovaryans en-boy oranı kıvrımlı solucanlarda ≈ 1 veriyor; ilk 11 koşuda bu alan hiç yoktu | İç teğet yarıçap ölçütüyle (R_eq / r_iç): 12 yuvarlak, 9 solucan/labirent, 1 solucan kopması, 2 büyük başlangıç bölünmesi, 6 çözündü (önce "17 tek damla, 15'i yuvarlak" deniyordu) | ben (yeniden sınıflama) + denetçi |
| 4 | orta | "L = 64 pilotu yuvarlak R ≈ 14,9" yanlış: kutuyu saran bir solucan. Makale hızlarındaki bütün 2B koşular solucan; solucan boyu kutuyla ölçekleniyor (en-boy 23 → 36) | pilotB_t30000: n = 0, süzülme var | denetçi (ii) |
| 5 | orta | A'nın en elverişli koşusu (A2_eq_A0.01) dambıldı ve T'de kesilmişti | Koşuyu t = 60 000'e uzattım: bölünmedi, solucana döndü (R_eq/r_iç 2,0 → 2,57 doygun, en-boy 10) | ben |
| 6 | orta | A'da güçlü şarj hücre sayısını artırıyor (12/12, 1 → 13–41), ama bölünmeyle değil: uzak alan spinodali aşıyor ve ortamdan doğum oluyor. Eski doğum kuralı bunları yanlışlıkla "bölünme" sayıyordu | Kural düzeltildi (önceki/sonraki uzak alan, NaN, N_yeni > 3·N_önceki). Eski sonuçların hiçbiri değişmedi | ben + denetçi |
| 7 | orta | Gray–Scott "rampası" aslında bir uzunluk ölçeği değişimi (D/2 ≡ kutu √2 büyük). Bu yüzden N'nin artması kaçınılmaz; asıl soru artışın bölünmeyle mi doğumla mı olduğu | Soy takibi: 3000 zaman biriminde doğum 0, net +44 artışın hepsi bölünmeden (eşik titremesi büyük: 325 bölünme / 234 birleşme) | ben |
| 8 | orta | Gray–Scott çoğalan 6/42 değil 5/42: (0,034; 0,059) labirent (süzülen) | Kural düzeltildi | ben + denetçi |
| 9 | orta | Fizik tahmini: olay başına enerji 45 meV değil ≈ E_g; 2B solucan ρ = 2,43 bir hücre boyu değil; kuantum hedefi tutarsız | Bkz. §0 ve §2 | ben + denetçi |
| 10 | düşük | §1 denklemi ve active.py açıklamasında Coulomb teriminin işareti "−∇²" yazılmış; kod doğru şekilde +∇²[G∗(c−⟨c⟩)] kullanıyor (itici) | — | ben |
| 11 | düşük | "703 adım ≈ Pearson'ın 660'ı" elmayla armut: 660, F = 0,02'de Lyapunov süresi; aynı uyum Pearson noktasında 1522 adım veriyor | — | denetçi |
| 12 | düşük | 3B kızlar "yuvarlak" değil: doğumda en-boy 2,6, sonra 5,1. "3,3 katı" 17,9 nm'ye göre 3,7 kat. Anlamsız ikilenme süreleri (10¹⁸–10¹⁹) koruma ile kaldırıldı. float32 ve %5 gürültü kontrolleri artık `checks_review.py` ile üretiliyor | — | ben + denetçi |

### 1.2 Exp10 (φ'nin fiziği, kapı haritası)

| # | Şiddet | Bulgu | Kanıt | Doğrulandı |
|---|---|---|---|---|
| 1 | kritik | Kapı haritası (`check.py`, with_corr=True) kapı ekranlı olmayan korelasyon vekiline dayanıyor; özette bu yazmıyor. Tutarlı HF'de kararsızlık hiçbir yerde yok. "Yakın kapı bölünmeyi açıyor / süperkapasitör sezgisi doğrulandı" cümlesi geri alınmalı | S/K_kin (ε_g = 3,9): HF'de en küçük +1,006 (C = 3000). Korelasyonla örneğin n = 10⁹, C = 0,83'te −22,4. Ana tablodaki kararsız varyantların hepsinde (A, B, D, E) değişim-korelasyon Hartree'den az ekranlanmış | ben (ızgara) + denetçi (i) |
| 2 | orta | "SKY130 C_ox'ta kararsızlık yok" ifadesi yalnız n ≥ 10¹⁰ taramasında geçerli. Vekille n_max: SiO₂ 4,148 nm 3,4·10⁹; 2 nm 9,8·10⁹; ε = 20 4,148 nm 7,2·10⁹ cm⁻². Haritadaki "yok" hücreleri tarama tabanının ürünü | — | denetçi (i), ben (ızgarada aynı işaretler) |
| 3 | orta | Aynı "kapı" iki zıt rol oynuyor, raporlarda ayrımı yazılmamış: (i) düzgün tabakanın domenlere ayrılması: Hartree dengeleyicidir, ekranlanması ayrışmayı kolaylaştırır (exp8 §2); (ii) var olan yüklü hücrenin ikiye bölünmesi: Coulomb iticidir, ekranlanması bölünmeyi engeller (exp8 §3) | — | denetçi (i) |
| 4 | düşük | "Hartree |K|'nın 11–27 katı": E varyantında Fang–Howard çarpanı unutulmuş, doğrusu 11–21 | F(k_F) = 0,78 | denetçi (i) |
| 5 | düşük | ε_c = (11,7 + ε)/2 seçimi belgelenmemiş; ε = 20 sonuçlarını belirliyor (7,8 seçilse n_max ≈ 4 kat artar). Tabakanın kendi kalınlığı (ağırlık merkezi 2,45 nm) seri ≈ 4,2 µF/cm² ekliyor; etkin kapı bağlaşımı bunu geçemez | — | denetçi (i) |

### 1.3 Aşama 2b Monte Carlo (Exp6) ve Exp9'un bu konudaki hükümleri

| # | Şiddet | Bulgu | Kanıt | Doğrulandı |
|---|---|---|---|---|
| 1 | kritik | Katkılar yalnız 126×420 nm kutuda, kapı ve taban düzlemleri sonsuz, S/D yok. Bu, ortalama potansiyelde sabit bir kubbe yaratıyor ve σ_ψ²'nin yarısından fazlası bu kubbe. Nokta yalıtım testi yalnız x kenarlarını dışlıyor, STI kenarı şeritleri "nokta" sayılıyordu | Campbell (aynı çekirdek): kubbe varyansın %55 / %59'u. Rastgele σ 13,8 / 20,2 mV (raporlanan 20,6 / 31,5). Merkez−kenar 48,6 / 75,8 meV, köşe 70,7 / 112,9 meV | ben (`review_campbell.py`) + denetçi (iii) |
| 2 | kritik | Düzeltilmiş yeniden koşu (aynı örnekler) | Orijinal kural 3/80. Dört kenar dışlanınca **0/80** (her iki katkı düzeyinde). Kubbe çıkarılmış (yalnız rastgele RDF): 31/80 (%39) ve 76/80 (%95); katı noktaların medyan boyu a = 41 / 35 nm | ben (`review_rerun.py`) |
| 3 | kritik | FEM düzeltmesi fiziksel sınırın dışına taşınmış. Arayüzde gerçek oran 2ε_Si/(ε_Si+ε_ox) = 1,5'i aşamaz (dielektrik görüntü yükü), ama uydurulan 1 + 1,8/z₀ z₀ = 1 nm'de 2,8 veriyor. Exp9'un "üst sınır" dediği %10–21 aşırı düzeltme | Denetçinin tam q-uzayı Green fonksiyonu FEM'i %0,3–2 içinde tutturuyor; rastgele σ 17,6 / 25,9 mV | denetçi (iii), ben (sınır argümanı) |
| 4 | orta | Exp9'un "σ_VT MC, Stolk ile %30 uyumsuz" hükmü yanlış. MC kutu ortalamalı ψ_s'nin dağılımını ölçüyor; beklenen W/(W+D_G) ve sonlu kutu çarpanlarıyla 3,72 / 4,71 mV, MC 3,94 / 4,52 mV, yani uyumlu. RDF, dökümhane eşleşme bütçesinin (13,4 mV) yaklaşık %42'si | — | denetçi (iii) |
| 5 | orta | 4 K uyarısı yanlış mekanizmaya işaret ediyordu. Tükenme bölgesindeki alıcılar 4 K'de de iyonize. Donan nötr gövde, ve bu gövde yüzer hâle gelir. Exp6 gövdeyi topraklı varsayıyor, Exp11/Exp12'nin B yolu ise aynı sıcaklıkta donmuş gövdeyi varsayıyor; çelişki | — | denetçi (iii), ben (fizik) |
| 6 | orta | Ekranlama çelişkisi. n_s ≈ 10¹² cm⁻²'de q_TF ≈ 1,84 nm⁻¹ ve 20 nm'lik dalgalanmalar yaklaşık 6,6 kat sönümlenir. Noktaların doluluk düzeyleri (15–45 meV) tam bu yoğunluğa karşılık geliyor. "Eşik altı, ekranlama yok" ile "eşik civarında Coulomb blokajı" aynı anda doğru olamaz | — | denetçi (iii) |
| 7 | orta | İstatistik. 3/80 için %95 aralığı %0,8–10,6 (raporlanan "± %2,1" 1σ Wald). "%4–21" nokta tahminlerinin yayılımı. "n" değerleri aygıt değil kayıt sayısı. "≥ 100 aygıt" %4'te üç pozitif için yetmez: 156 gerekir | — | denetçi (iii) |
| 8 | düşük | "nimg 8 ile 12 farkı ≤ 0,26 mV" hiçbir sonuç dosyasında yok; kayıtlı kontrol 4 ile 8 arası ve 1,20 mV. check1 kendini kendisiyle karşılaştırıyor | `rdf.json` checks | ben |

### 1.4 Raporlar arası

| # | Şiddet | Bulgu |
|---|---|---|
| 1 | orta | Aynı nokta için E_C farklı raporlarda 3,7 kat farklı. a = 20 nm için: FEM hücresi ε = 11,7 → 19,5 K; Exp6 sabiti → 29 K; FEM hücresi ε = 3,9 → 58,5 K; Exp11 izole küre → 63–71 K. Düzlemsel bir çukurun kapı kapasitansı C ≥ C_ox·a² = 3,3 aF olduğundan E_C üst sınırı ≈ 280 K. Projede tek bir FEM tabanlı "düzlemsel çukur" kapasitansı kullanılmalı |
| 2 | orta | 4 K gövde sınır koşulu: Exp6 topraklı ve iletken gövde, Exp11 ve Exp12 donmuş ve yüzer gövde. Exp12 yüzer gövdeyi esas alacak; Exp6 yeniden koşulursa aynı koşul kullanılmalı |

## 2. Yeni doğrulama hesapları

| Hesap | Sonuç | Dosya |
|---|---|---|
| 3B Δν taraması (ε = 0,20, 64³) | Δν −0,002 … +0,004: hepsi t ≈ 2000'de ikiye bölündü; bölünme yarıçapı Δν arttıkça 11,8 → 6,6 w azalıyor. Δν = −0,004: t = 3000'de tek, R = 19,8 | `b3d_*.json` |
| A şarj şiddeti (J = 10⁻⁴, 5·10⁻⁴, 2·10⁻³; 4 Coulomb/kapı durumu) | 0/12 bölünme, 12/12 ortamdan doğum; Coulomb'suz kontrolde de aynı | `runs/A3_*` |
| A2_eq_A0.01 uzatma (t = 20k → 60k) | bölünme yok; dambıl solucana dönüştü (R_eq/r_iç 2,57, en-boy 10) | `runs/A2_eq_A0.01_dinf__ext.*` |
| GS soy takibi (D/2 sonrası 3000 zaman birimi) | N 83 → 127; doğum 0, ölüm 44; artış bölünmeden | `gs_track.json` |
| 2B tabaka güç/boy (en iyimser köşe, N_min = 100) | E = 45 meV: 55 nm; E = E_g: 123 nm; N_min = 10 ile 31 / 69 nm. Tipik (D = 1, 1 W/cm²): 175 / 390 nm | `phys.json` |
| Exp6 Campbell + düzeltilmiş yeniden koşu | Bkz. §1.3 | `results/exp6/review_*.json` |

## 3. Düzeltilmiş karşılaştırma (hücre üreten mekanizmalar)

| Mekanizma | 2B'de (MOSFET kanalı gibi) | 3B'de | Hücre boyu ve kuantum ilgisi (4 K) | Geçerlilik |
|---|---|---|---|---|
| A: süperkapasitör (yüklü hücre + kapı + şarj) | 0/28 kalıcı bölünme; uzar/solucan olur; güçlü şarjda ortamdan doğum | ekransız 3B Coulomb (≡ Δν = 0) **bölünür** | silikon karşılığı yok: tutarlı HF'de tabaka hiçbir kapıda ayrışmıyor | 2B kapı ekranlaması itmeyi öldürüyor |
| B1: aktif damla (Zwicker, disipatif) | 0/30 ikili bölünme; 9 solucan, 1 kopma | bölünür (Δν'nin geniş aralığında) | — | 2B'de hücre değil filaman/solucan |
| B2: Gray–Scott (otokatalitik, disipatif) | **çoğalır** (5/42 (F,k), 6/6 dayanıklı, doğum 0) | — | silikon karşılığı en iyimser 123 nm (E_g), tipik ≈ 0,4 µm; hedef 6–20 nm | sınır tünel bariyeri değil, taşıyıcı sıcak |
| C: kazara nokta (RDF, Aşama 2b; bölünme değil) | orijinal tahmin geri çekildi (kenar artefaktı, 0/80); yalnız rastgele RDF: %39–95, a ≈ 35–41 nm | — | E_C/kT ≈ 3,5–4 (Exp6 sabiti), plaka sınırıyla ≲ 17 | ekranlama yok; kesin sayı için öz-uyumlu model gerekir |

Okuma: "Disipatif mi, süperkapasitör mü?" sorusunun cevabı tek bir kazanan değil. Bölünmeyi belirleyen üç şey var: (1) boyut, (2) uzun menzilli itme ve (3) otokatalitik aktivatör–inhibitör yapısı. MOSFET kanalı 2B ve kapılı olduğu için orada yalnız (3) sınıfı çalışıyor. Kuantum ölçeğine ise hiçbiri inmiyor.

## 4. Yapılan düzeltmeler

- `results/exp11/summary_exp11.md`: §0, §3, §4, §5, §6 düzeltildi; başa düzeltme notu eklendi; Coulomb işareti düzeltildi.
- `src/qcell/active.py`: açıklamadaki işaret düzeltildi (kod değişmedi).
- `experiments/exp11_compare_division/compare.py`: şekil ölçütü (iç teğet yarıçap), doğum kuralı, ikilenme koruması, GS labirent kuralı; yeniden üretilen `classification.json`, `stats.json`, `comparison.csv`, şekiller.
- `results/exp10/summary_exp10.md`: başa düzeltme kutusu (kapı haritası vekile dayanıyor; tutarlı HF kararlı; 11–21; n_max).
- `research/roadmap/ROADMAP.md`: §3c (Aşama 2b) düzeltme notu ve geri çekme, §3d/§3e düzeltmeleri.
- `research/hypotheses/hypotheses.md`: H10, H16, H17–H21 güncellendi.
- `results/exp9/COMPARISON.md`: başa not (G7, H, FEM "üst sınır" hükümleri düzeltildi).

## 5. Açık kalanlar

- Aşama 2b için öz-uyumlu (Thomas–Fermi) ekranlama, gerçek geometri (S/D elektrotları, STI) ve tam q-uzayı çekirdeğiyle yeniden yapılmış bir model. Bu yapılmadan kazara nokta oranı için sayı verilmemeli.
- Projede tek bir FEM tabanlı düzlemsel çukur kapasitansı.
- 2B Zwicker solucanlarının kopması kutu boyutunun bir ürünü olabilir. Ayırt etmek için L ≥ 256 ve T ≥ 5·10⁴ koşusu gerekir (yapılmadı).
- Exp12: 4 K yüzer gövde + darbe iyonlaşması (E ≈ E_g) ile akım filamanı; aşağıdaki Exp12 raporunda.
