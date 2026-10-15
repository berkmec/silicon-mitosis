# Hipotezler ve durumları

Fikir: elektrikle uyarılan transistör+süperkapasitör hücreleri, içlerinde daha küçük hücreler doğurur;
en küçük nesil kuantum rejimine geçer.

Not: Deney 0 (`results/exp0`) yazdığım analitik bir ön hesaptı, fiziksel simülasyon değildi.
Aşağıdaki durumlar Deney 1'e (`results/exp1`, 3B FEM, doğrulanmış) dayanır; kapsam dışı olanlar açıkça belirtildi.

| # | Hipotez | Durum (Deney 1) |
|---|---|---|
| H1 | Nesil küçüldükçe E_C ∝ 1/L büyür ve bir nesilde E_C ≥ 10 kT olur | FEM'de C ∝ a (ölçek değişmezliği), r=0.5 için c_Σ ≈ 23 ε0·a (±~%1,2). 4 K: nesil 2 (25 nm, ε=3,9) veya nesil 4 (6 nm, ε=11,7). 77 K: yalnız <1 nm'de (model geçersiz). 300 K: yok. **Kısmen destekleniyor (yalnız kriyojenik).** |
| H2 | Yüksek kapasitanslı "süperkapasitör" hapsetmeyi güçlendirir | FEM'de ayrıca test edilmedi. E_C = e²/2C ilişkisi gereği C arttırmak E_C'yi düşürür; hücrede c_g toplamın ≈%76'sı. Deney 0'daki sonuç (ters etki) fizik kimliğiyle tutarlı. **Desteklenmiyor; doğrudan FEM testi yapılmadı.** |
| H3 | Çoğalma kuralı derinlemesine devam eder (enerji/kırılma penceresi) | Çoğalma kuralı FEM'e taşınmadı; Deney 0'ın analitik kuralı yalnız bir ön hesap. **Sınanmadı.** |
| H4 | Seviye ayrışması (hapsolma) Coulomb'dan önce gelir | FEM (Si/SiO2, V0=3,1 eV): ΔE sert duvar formülünden küçük (3 nm'de 0,69×, 1 nm'de 0,34×). 300 K'de ΔE ≥ 10 kT nesil 5'te (3,1 nm), marj küçük. **Destekleniyor, daha zayıf biçimde.** |
| H5 | İç içe iletken kabuk adayı "ekranlar" / bağlaşımı değiştirir | Yüzer kabuk eş-potansiyel; adanın C'si pencereden yalnız ≈%1 etkileniyor, kabuk ise C'yi %12 artırıyor. **Desteklenmedi.** |

## Deney 2 ve 3 (zincir, çift nokta, oda sıcaklığı qubit taraması) sonuçları

| # | Hipotez | Durum |
|---|---|---|
| H6 | Boyutu kademeli küçülen nokta zinciri, küçük hücreleri sırayla yükler ve en küçük hücre kuantum rejimine geçer | FEM kapasitans + ana denklem (12/6/3 nm): 4 K'de ada 0 ve 1 yük basamakları gösteriyor, **3 nm'lik en küçük ada hiç yüklenmiyor**; 77 K ve 300 K'de yük sürekli (kuantizasyon ≈ 0), akım salınımı yok. **Desteklenmedi (oda sıcaklığında).** |
| H7 | Çift hücrede (iki nokta) 300 K'de başlatılabilir ve koherant qubit oluşturulabilir | FEM Schrödinger + tek-fonon akustik altın kural, 9 geometri: 300 K polarizasyon ≥ 0,9 olan yok (en iyi 0,74, L=1 nm, d=0,5 nm, Δ=50 meV). En iyi geometride bile X-kapı doğruluğu 0,78 (τφ=100 fs, iyimser hızlı kapı); elektriksel 10 ps kapıda 0,5 (ölü). **Gösterilemedi.** |
| H8 | Fonon darboğazı büyük Δ'da koheransı korur | Δ=50 meV'de <F²> ≈ 8e-5 (q₀ ≈ 9 nm⁻¹ nokta boyutundan çok küçük dalga boyu), Δ≈10 meV'de ≈ 0,3-0,4. Yani akustik fonon etkisi büyük Δ'da gerçekten bastırılıyor. **Destekleniyor, yalnız akustik tek-fonon kanalı için.** |

Önemli: En iyi geometri taramanın kenarında (en küçük L, en ince d). 1 nm altında sürekli ortam modeli geçersiz; atomistik yöntem gerekli.
Optik fonon, Raman, yük gürültüsü ve vadi/orbital kanalları hesaplanmadı: T2 değerleri ÜST sınırdır. 300 K'de τφ = 10 fs/100 fs/1 ps senaryoları tahmin değil, varsayımdır.

## SKY130 yığını (Aşama 2a-2b), yol haritası: `research/roadmap/ROADMAP.md`

| # | Hipotez | Durum |
|---|---|---|
| H9 | SKY130 nFET kanalında dikey olarak 2B elektron gazı oluşur (kuantize alt bantlar) | 1B Poisson-Schrödinger (toxe 4,148 nm, ndep 1,7e17): tabaka ağırlık merkezi 2,45 nm, E₀ ≈ 50 meV (n_s = 10¹² cm⁻²). **Destekleniyor (model).** |
| H10 | Rastgele katkı potansiyeli minimum aygıtta kazara kuantum nokta oluşturur | MC (4 varyant, 3B FEM kontrollü): cihazların ≈ %4–21'inde bağlı durumlu nokta (a ≈ 15–23 nm, E_C/k_B ≈ 25–40 K). 3B FEM görüntü modelinden %18'e kadar büyük potansiyel verdi (düzeltme eklendi). Tünel bariyeri modellenmedi, istatistik küçük. **Geçerlilik sorunlu (exp9, bkz. `results/exp9/COMPARISON.md`):** düzeltmenin z<10 nm kısmı kesilince σ_ψ %23 düşüyor, katı oran %2,5–10'a iniyor; model oranı T'den bağımsız sayıyor, ama 4 K'de katkılar tam iyonlaşmıyor (NA = 1,7e17 ≪ N_c ≈ 3,7e18, bellekten). Doğrulanmadı; 4 K tahmini yeni model çıkana kadar kullanılmamalı.** **İnceleme (2026-10-11): oran tahmini geri çekildi.** Noktalar kutu kenarı artefaktı: kubbe σ_ψ²'nin %55–59'u, dört kenar dışlanınca 0/80. Yalnız rastgele RDF: 31/80–76/80, a ≈ 35–41 nm, E_C/kT ≈ 3,5–4. FEM düzeltmesi fiziksel 1,5 sınırını aşıyor; ekranlama yok. 4 K'de tükenme bölgesi alıcıları iyonize, asıl sorun yüzer gövde. **Belirsiz; öz-uyumlu modelle yeniden yapılmalı.** |
| H11 | V_T tanımı BSIM4 vth0 ile uyumlu | Sabit yük tanımı (n_s=10¹² cm⁻²) V_g = 0,513 V vs vth0 = 0,519 V; doğrusal ekstrapolasyon 0,34 V. **Tanıma bağlı; belirsiz.** |

## "Bölünme" mekanizmaları (Deney 8, soyut model)

| # | Hipotez | Durum |
|---|---|---|
| H12 | Yüklü bir küme, yük belirli bir eşiği aşınca iki parçaya bölünmek enerjice avantajlıdır (Coulomb fisyonu) ve eşik N^−1 ile ölçeklenir | LJ kümesi + Coulomb: g_c ∝ N^−1,08 (d=∞). g ızgarası ×1,45, eşikler ±%45 belirsiz. **Destekleniyor.** |
| H13 | Gate (kapasitör) uzaklığı bölünmeyi/hücre boyutunu belirler | Domain dalga boyu d=1,5'te 1,5–1,7× büyük; A_c(d) gate yaklaştıkça artıyor; φ̄=−0,4, A=0,15'te gate yaklaşınca desen başlıyor. **Destekleniyor (soyut model).** |
| H14 | Simetrik durum kararsızlaşınca sonlu boyutlu hücrelere ayrılır (bifurkasyon) | Süperkritik çatal bifurkasyonu: A_c analitik = sayısal, genlik üssü 0,502. **Destekleniyor (doğrulanmış).** |
| H15 | Hazır tek bir hücre (domain) Coulomb frustrasyonuyla ikiye bölünür | Disk, çubuk, elips: hiçbiri bölünmedi (küçüldü/çözündü). **Çürütüldü (bu 2B modelde, yalnız dengede).** Exp11: sürülen (disipatif) sistemde tek hücre bölünebiliyor (H17). Exp11'deki Coulomb + şarj koşularında da (28 koşu) bölünme yok (H18). |
| H16 | Bu mekanizmalar SKY130 evirme tabakasında geçerlidir | Exp10: φ ≈ (n − n₀)/n₀. Yerel ters kompresibilite n_s ≲ 1,5–3·10¹² cm⁻²'de negatif (çift kuyu fiziksel olarak var), ama kapı ekranlı Hartree k_F'ye kadar |K|'nın 11–21 katı (inceleme: 27 değil, Fang–Howard çarpanı unutulmuştu): düzgün tabaka kararlı, bölünme yok. Kararsızlık yalnız n ≲ 3·10⁹–2·10¹⁰ (r_s ≳ 20–50), model geçerliliği dışında. **Çürütüldü (SKY130 çalışma aralığında, dengede).** Yeniden kontrol: kod Thomas–Fermi analitiğiyle birebir; bölünme yalnız kapı ≤1 nm (SiO₂) ya da ≤2 nm (ε=20) ve n ≲ 1–5·10¹⁰ cm⁻² (r_s ≈ 12–26) için çıkıyor; hücre boyutu kesme sınırına bağlı, güvenilmez. Exp11: tutarlı HF (korelasyon olmadan) hiçbir kapasitansta (0,05–3000 µF/cm²) kararsızlık vermiyor; korelasyonlu sonuç, kapı ekranlı olmayan vekile dayanıyor (model içi tutarsızlık). **İnceleme (2026-10-11): "yakın kapı bölünmeyi açar" sonucu bu tutarsızlığın ürünü.** Tutarlı HF'de S/K_kin ≥ +1,006, n_s = 10⁹–10¹³ cm⁻² ve her ε_kapı için. **Çürütüldü (her kapıda, tutarlı modelde).** |

## Sürülen (disipatif) bölünme: Deney 11 (`results/exp11/summary_exp11.md`)

| # | Hipotez | Durum |
|---|---|---|
| H17 | Disipatif sürüş (Prigogine) tek bir hücreyi çoğaltır | Gray–Scott: 5/42 (F,k) noktasında çoğalma (inceleme: 6 değil, biri labirent), %5 gürültü ve 3 tohumla 6/6. Aktif damla 2B'de ikili bölünme 0/30 (9 solucan, 1 kopma); 3B'de Δν −0,002…+0,004 aralığının hepsinde ikiye bölünüyor, **Δν = 0 dahil** (denge eşdeğeri). **Gray–Scott için destekleniyor; aktif damlada bölünme disipasyona özgü değil (sürüm 2).** |
| H18 | Süperkapasitör fikri (yüklü hücre + kapı ekranlı Coulomb + elektrikle şarj) tek hücreyi böler | 2B'de 0/28 kalıcı bölünme. Dambıl koşusu t = 60 000'e uzatıldı, bölünmedi. Güçlü şarjda 12/12 ortamdan doğum var, bölünme yok. 3B'de ekransız Coulomb'lu yüklü damla (≡ Δν = 0) bölünüyor. **2B kapılı tabakada çürütüldü; 3B ekransız durumda destekleniyor (sıvı damla fisyonu).** |
| H19 | Sızıntılı şarj (s = J − λc), disipatif modelin Δν = 0 çizgisine (Ohta–Kawasaki λ/q²) eşdeğerdir | Cebirsel olarak gösterildi ve testle doğrulandı. 2B'de bu çizgide 0/5 (ε = 0,20 noktası T'de hâlâ değişiyordu); **3B'de bölünüyor** (iki eşit kız, R = 8,8). **Destekleniyor (eşdeğerlik); 3B'de bölünme için yeterli.** |
| H20 | Sürüş arttıkça hücreler bölünerek küçülür ("küçülerek çoğalma") | Gray–Scott'ta D/2 ≡ uzunluk ölçeğinin √2 küçülmesi; N'nin artması bu yüzden kaçınılmaz. Gösterilen şey artışın **bölünmeyle** olması: soy takibinde doğum 0, N 83 → 127. Taban ℓ. **Destekleniyor (Gray–Scott; mekanizma olarak bölünme).** |
| H21 | Disipatif hücreler silikonda 4 K kuantum boyutuna iner (R ≈ 6–20 nm, geometriye bağlı) | Olay enerjisi ≈ E_g (kink mekanizması) ile en iyimser 123 nm (2B tabaka), ≈ 1,6 µm (3B), tipik ≈ 0,4 µm. 45 meV varsayımı kink literatürüyle desteklenmiyor. **Desteklenmiyor (mertebe tahmini).** |
| H22 | 4 K'de geniş bir SKY130 nFET'te darbe iyonlaşması akım filamanları (disipatif desen) oluşur | Exp12 (ngspice + resmî BSIM4 + indirgenmiş 4 K modeli). Eşik üstünde geri besleme zayıf, histerezis yok. Eşik altında yüzer gövde kilitlenmesi var: histerezis 26–105/108 kombinasyonda. Donmuş kuyu R_□ ≳ 10¹⁰–10¹² Ω/□ (ℓ ≪ W) iken yük direnciyle tek bir akım filamanı oluşuyor; V_ds analitik Maxwell geriliminde (0,80 V) kenetli, filaman akımla genişliyor (47 → 94 µm) ve bölünmüyor. **Model düzeyinde ve koşullu destekleniyor; klasik, kuantum değil; ölçülmeli.** Exp13: 5 parametre setinin 5'inde tek filaman; R_□ = 10¹⁰ dy ≤ 0,1 µm'de yakınsadı. |
| H23 | Eşik altı kilitlenme ortamında filaman, yerel bir inhibitör (dirençli drenaj) varsa Gray–Scott gibi bölünür | Exp13: yalnız global bağlaşımla bölünme yok. Uçlardan beslenen hatta filamanlar kontaklarda doğuyor; aşağı taramada ortadan bölünme var. Dağıtık balastta (λ = 40 µm, 5 kΩ, W = 300 µm) büyüme → ikiye bölünme 3/3 tohumda (≈ 4λ'da; merkezde V_D en düşük). Tekrarlayan çoğalma yok; kızlar 70–100 µm. **Model düzeyinde destekleniyor (tek bölünme); tekrarlayan çoğalma ve küçülme yok.** Exp16: λ = 5–20 µm'de filaman 203 µm'ye bölünmeden ulaştı, yani "≈ 4λ" genel bir yasa değil. |
| H24 | Aktif damlada bölünme boyuta bağlıdır (6B hiperküpte Cahn–Hilliard ile gösterilebilir) | Exp14: boyut sürekli değişken yapıldı (O(d−1) indirgeme; 3B Kartezyen ve gerçek 6B hiperküp 16⁶ ile doğrulandı). Bölünme yalnız bir boyut penceresinde oluyor. Altında şerit var (2B), üstünde damla bölünmeden içten boşalıp sıvı kabuk oluyor. Pencere Δν ile yukarı kayıyor: −0,004'te 2,25–2,5; −0,001'de 2,3–3,5; +0,004'te 2,5–4; Δν ≥ +0,012'de 6B de bölünüyor. Pasif RP bandı kR < √(d−2) d = 3–6'da %2,5–4 içinde doğrulandı. **Destekleniyor, ama tekdüze değil: yüksek boyut daha güçlü sürüş istiyor.** |
| H25 | Gerçek bir fiziksel 6B faz uzayı sisteminde (Vlasov–Poisson) yapılar boyuta göre bölünür | Exp15 (doğrulanmış PIC): elektron delikleri 1B'de yalnız birleşiyor. 2B/3B'de kıvrılıp parçalanıp sönüyor; kalıcı bölünme 0–1. Eşik 1B↔2B (Morse & Nielson 1969 ile uyumlu); kuvvetli B kırılmayı bastırıyor. **Boyuta bağlılık var ama bölünme yok (yok olma); silikonda 4 K'de model geçersiz.** |
| H26 | Proje verilerinden seyrek regresyonla (PySINDy) denklem türetilebilir | Exp16: bilinen yasalar geri bulundu (Gray–Scott PDE ≤ %0,19; BSIM4 Chynoweth 90/90; HF −1,65335; BenDaniel–Duke). %1 gürültüde başarısız. Aday yasalar: R*₂B = 0,233 ε/ν₋ (LOO %3); filaman w = I/j_f, j_f = 2,55 µA/µm; fisyon g_c ∝ N^≈−1 (Exp8'deki "kapı üssü düzleştirir" yorumu desteklenmiyor). **Temiz simülasyon verisinde destekleniyor; gürültülü veride değil.** |
| H27 | Silikonda (SKY130 nFET, 4 K) dağıtık balastlı drenajlı geniş cihazda akım filamanı büyüyüp tekrar tekrar bölünür; bölünme genişliği cihaz sabitlerinden öngörülebilir | Exp17: genel ölçüt merkez V_D = V_down (merkez–kenar farkı = V_M − V_down). Tek filaman için w_s = 2λ ln[λ/(λ−λ_c)], λ_c = √(2(V_M−V_down)/(r_D j_on)); dizi için kapalı form (sınanmadı); w_top = I/j_on. Öngörüler önceden kaydedildi. Temiz protokolün 22 bölünmesinde oran 1,001 ± 0,069 (63 olayın tamamı Exp17 düzeltme kutusunda), iki farklı cihaz setinde 0,979 ± 0,133 (n = 9); 14 ilk bölünmede RMS %9; λ = 36 µm (eşiğin %1 üstü) ıskalandı. Fizik formu veriden tek katsayıyla uydurulunca kuramın 1,01 katı. Sınır: yarı statik sürüş gerekli (düşük sızıntıda gecikme). 2B parmak dizisinde (duvarlı ve duvarsız) leke bölünmüyor; yuvarlak → dambıl → solucan. Bölünme için 1B'ye hapsetmek gerekiyor. **Model düzeyinde güçlü destek (1B); tekrarlayan bölünme var; ölçülmeli. Kuantum değil.** |


## Açık / sonraki deney adayları
- Hücreler arası tünel bağlaşımıyla koherans (hiçbir deneyde yok).
- Çoğalma fiziği Exp11'de modellendi (disipatif: Gray–Scott, aktif damla). Silikon karşılığı açık: 4 K darbe iyonlaşması filamanları (H22), e-h damlaları.
- Poisson-Schrödinger öz-uyum döngüsü ve çok-elektron etkileri (şu an tek-parçacık, bağımsız).
- <1 nm için atomistik yöntem (tight-binding / DFT); sürekli ortam burada geçersiz.
