# Yol haritası: silikon üzerinde kuantum cihaza giden yol (SKY130 hedefiyle)

Bu belge `prompt.txt` isteğine cevaptır. Her bölümde "ölçüldü / literatür / varsayım / bilinmiyor" ayrımı yapıldı.
Sayısal kaynaklar: `results/exp1..exp4`. Model sınırları: `research/notes/model_assumptions.md`.

## 1. Kısa cevap

Hedefin mevcut hâli, yani **oda sıcaklığında, kriyojeniksiz, standart 130 nm SkyWater CMOS'ta kendiliğinden çoğalan hücrelerle kuantum rejimi**,
ölçtüğüm fizikle uyuşmuyor. Üç bağımsız sebep var:

1. **Ölçek:** Oda sıcaklığında Coulomb blokajı (E_C ≥ 10 kT) için ada kenarı ≈ 0,4 nm (ε=3,9) veya ≈ 0,1 nm (ε=11,7) gerekiyor.
   SKY130'daki en küçük poly genişliği 150 nm. Fark ≈ 400–1500×. 77 K için bile 1,5 nm / 0,5 nm gerekiyor (aşağıda tablo).
2. **Bağlaşım:** Hücreler arası tünel bağlaşımı için nm ölçeğinde aralık gerekiyor. SKY130'da poly-poly aralığı 210 nm; silikon içinde 20 meV'lik
   elektrostatik bariyer için tünel genliği 210 nm'de 10⁻²⁹ çarpanına iniyor (WKB). Yani iki kapı-tanımlı hücre arasında bağlaşım yok.
3. **Koherans:** Oda sıcaklığında en iyi geometride bile (L=1 nm, d=0,5 nm) başlangıç polarizasyonu 0,74, kapı doğruluğu ≤ 0,78 (iyimser) ve
   10 ps'lik elektriksel kapıda 0,5 (ölü). Yalnızca tek-fonon akustik kanal hesaplandı, yani gerçek T2 bunun altında.

Ulaşılabilir, doğrulanabilir bir alt hedef var: **SKY130'da en küçük aygıtlarda (150 nm kapı) ≤ 1 K'de tek elektron / Coulomb blokajı davranışı
aramak** (aşağıda Aşama 3-4). Bu, "qubit" değil, "SET davranışı var mı?" sorusudur.

## 2. Bu projede ölçülenler (kendi simülasyonlarım)

| Konu | Sonuç | Doğrulama |
|---|---|---|
| Kapasitans (ada+S/D+kapı kabuğu) | c_Σ ≈ 23 ε0·a (r=0,5), belirsizlik ~%1,2 | küre analitiği (eğim 1,91), dolfinx ile 8e-16 uyum |
| Coulomb eşiği | E_C/k_B ≈ 1171 K·nm / a (ε=3,9) | FEM'den çıkan sabit, ölçek değişmezliği |
| Seviye ayrışması | Sert duvar formülünden küçük (3 nm'de 0,69×, 1 nm'de 0,34×) | sert duvar ve ayrılabilir sonlu kuyu analitik |
| Tünel yarılması | Δ ∝ exp(−κ d), uyum: FEM/WKB = 1,11 / 1,06 / 1,03 (L=1/1,5/2 nm) | WKB analitik |
| Landau-Zener | sayısal ≡ analitik, hata ≤ 1,2e-3 | γ = 0,1…4 |
| Zincir (12/6/3 nm) | 4 K'de yük basamakları; 77 K ve 300 K'de yük sürekli | FEM+ana denklem, tek-ada limiti testi |
| Çift nokta, 300 K | P ≤ 0,74, T1(akustik) fs–µs, kapı ≤ 0,78 | Aer ≡ kanal hesabı |
| Fonon darboğazı | Δ=50 meV'de ⟨F²⟩≈8e-5, Δ≈10 meV'de ≈0,3–0,4 | mesh ve kareleme yakınsaması |

Hipotez durumu: `research/hypotheses/hypotheses.md` (desteklenen, çürüyen ve sınanmayanlar ayrı).

## 3. Ölçek boşluğu (SKY130)

E_C ≥ 10 kT için gereken en büyük ada kenarı (r=0,5 hücre geometrisi, FEM sabiti):

| T (K) | ε=3,9 | ε=11,7 |
|---|---|---|
| 0,05 | 2342 nm | 781 nm |
| 0,1 | 1171 nm | 390 nm |
| 0,3 | 390 nm | 130 nm |
| 1 | 117 nm | 39 nm |
| 4 | 29 nm | 9,8 nm |
| 20 | 5,9 nm | 2,0 nm |
| 77 | 1,5 nm | 0,5 nm |
| 300 | 0,4 nm | 0,1 nm |

- 150 nm'lik bir ada için E_C/k_B ≈ 7,8 K, yani 10 kT ölçütü ≈ 0,78 K'de sağlanır. Bu, seyreltme soğutucusu (dilution) gerektirir.
- Bu, kapı-tanımlı **tek bir** büyük ada için. Kanalda kusur/pürüz kaynaklı **kazara (accidental) nokta** gate uzunluğundan çok küçük olabilir
  ve E_C daha büyüktür. Bu bir tahmin değil, literatürde 28/40 nm ticari aygıtlarda gözlenen davranış (bkz. §4).
- Silikonda elektrostatik bariyer için WKB (m*=0,19): bariyer 5/20/50 meV → bozunma uzunluğu 6,3/3,2/2,0 nm; 30 nm'de t oranı 9e-3 / 8e-5 / 3e-7;
  210 nm'de 4e-15 / 2e-29 / 3e-46. SKY130'un poly aralığı bu yüzden hücre bağlaşımı için kullanılamaz.

## 3b. SKY130 yığını ve 1B Poisson-Schrödinger (Aşama 2a)

Parametreler resmî model dosyasından (`sky130_fd_pr__nfet_01v8__tt.pm3.spice`, [skywater-pdk-libs-sky130_fd_pr](https://github.com/google/skywater-pdk-libs-sky130_fd_pr)):
toxe = 4,148 nm, epsrox = 3,9, ndep = 1,7×10¹⁷ cm⁻³, xj = 150 nm, vth0 = 0,519 V, lint = 11,9 nm (çizilen 150 nm kapıdan etkin kanal ≈ 126 nm).
Rastgele-katkı eşleşmesi: vth0_slope = 3,356e-3 (birimi dosyadan kesinleştiremedim; V·µm kabul edilirse σ(150 nm × 420 nm) ≈ 13 mV).

Çözücü doğrulaması: Airy üçgen kuyu (hata < %0,5), derinlik yaklaşımı (hata ≲ kT mertebesi), sabit/monoton yük (3 birim test).

| Büyüklük | Sonuç (ndep = 1,7e17) | Not |
|---|---|---|
| Fermi potansiyeli, derinlik genişliği | φ_F = 0,430 V, W_dep = 81 nm | analitik |
| C_ox | 0,833 µF/cm² | analitik; C_inv(PS) = 0,75 µF/cm² |
| V_T (klasik derinlik formülü) | 0,135 V | n+ poly, ideal |
| V_g, n_s = 10¹² cm⁻² (sabit yük) | 0,513 V | BSIM4 vth0 = 0,519 V'a %1 yakın |
| V_T (doğrusal ekstrapolasyon) | 0,338 V | vth0'dan 0,18 V düşük |
| Tabaka ağırlık merkezi (n_s = 10¹²) | 2,45 nm | alt bant: unprimed E₀ = 50 meV, primed 101 meV (E_F üstü) |
| Yüzey alanı | 0,38 MV/cm | |

Okuma: V_T tanımına bağlı olarak vth0 ile uyum **tanım seçimine bağlı** (sabit yük: uyumlu; doğrusal ekstrapolasyon: 0,18 V sapma). Hangi tanımın
BSIM4'ün vth0'ı olduğunu burada kesinleştiremedim, bu yüzden "etkin katkı 3,9e17" sonucu (ekstrapolasyon tanımına kalibre) bir bulgu değil, tanıma bağlı
bir kalibrasyondur. Dikey kuantizasyon (tabaka ≈ 2–2,5 nm, E₀ ≈ 50 meV) elektron gazının 2B olduğunu doğruluyor; dolayısıyla SKY130 kanalında
kuantum nokta oluşumu için gereken yanal sınırlama ek olarak kusur/katkı potansiyelinden gelmelidir.

Kaba tahmin (doğrudan hesap değil): FEM sabitinin ε ölçeklemesiyle, ε_eff ≈ 7,8 için E_C/k_B ≈ 585 K·nm / a. a = 10 / 20 / 30 nm → 59 / 29 / 20 K;
4 K'de E_C ≥ 10 kT için a ≲ 15 nm. Dikkat: bu sabit, ada + kaynak/savak + saran kapı boşluğu geometrisinden geldi; bir MOSFET kanalındaki kazara nokta
için geometri uymuyor, yalnızca büyüklük mertebesi gösterir. Aşama 2b bunu gerçek geometri ile sınar.

## 3c. Rastgele katkı potansiyeli ve kazara nokta (Aşama 2b)

> **Geri çekme (2026-10-11 incelemesi, `results/exp11/review_exp6_exp10_exp11.md`).** Aşağıdaki kazara nokta oranları ve
> "%4–21" tahmini **geçersiz**.
> - **Kenar artefaktı.** Katkılar yalnız kutu içinde, S/D elektrotu yok. Bu, ortalama potansiyelde sabit bir kubbe yaratıyor:
>   σ_ψ²'nin %55/%59'u, kenarlar 49/76 meV aşağıda. Yalıtım testi STI kenarlarını dışlamıyordu. Aynı örneklerde dört kenar
>   dışlanınca **0/80**.
> - **Yalnız rastgele RDF potansiyeli (kubbe çıkarılmış).** Çukurlar sık ama büyük: 31/80 ve 76/80, medyan a = 41/35 nm,
>   E_C/kT ≈ 3,5–4 (4 K).
> - **FEM düzeltmesi aşırı.** Fiziksel sınır 1,5, uydurulan formül 2,8'e kadar çıkıyor.
> - **Ekranlama çelişkisi.** n_s ≈ 10¹²'de 20 nm'lik dalgalanmalar yaklaşık 6,6 kat sönümlenir; model ekranlama içermiyor.
> - **4 K gövde.** Topraklı değil, yüzer.
>
> Öz-uyumlu ekranlama ve gerçek geometri olmadan sayı verilmemeli. σ_VT (MC) aslında Stolk ile uyumlu; exp9'un "%30 uyumsuz"
> hükmü yanlıştı.

Model: iyonlaşmış akseptörler (−q) derinlik bölgesinde Poisson dağılımlı; arayüzdeki potansiyel, topraklı kapı düzlemi (silikon-eşdeğer mesafe
4,148 × 11,7/3,9 = 12,4 nm) ve nötr yığın düzlemi arasında görüntü yöntemiyle (nimg = 8; kayıtlı yakınsama kontrolü nimg 4 / 8 = 1,20 mV) toplandı. Taşıyıcı ekranlaması yok
(eşik altı için makul), kuantum hapsolma için ayrı bağlı-durum ölçütü kullanıldı. Doğrulama: tek yük + görüntü analitik (fark < 1e-9 V).
80 örnek / katkı seti; nokta sayıları küçük, bu yüzden oranların belirsizliği büyük (3/80 = %3,8 ± %2,1).

| Büyüklük | ndep = 1,7e17 (W_dep 81 nm) | N_eff = 3,9e17 (W_dep 55 nm) |
|---|---|---|
| Ortalama katkı sayısı (kanal hacmi) | 728 | 1133 |
| σ (arayüz potansiyeli) | 20,6 meV | 31,5 meV |
| σ_VT (MC) | 3,9 mV | 4,5 mV |
| σ_VT Stolk formülü (literatür, burada türetilmedi) | 5,7 mV | 7,1 mV |
| σ_VT BSIM4 eşleşme parametresi | ≈ 13 mV (birimi kesin değil) | |
| Klasik izole birikinti içeren cihaz oranı | %48 | %61 |
| Bağlı durum tutan nokta, katı ölçüt | %3,8 (n=3), a ≈ 23 nm, E_C/k_B ≈ 25 K, E₁ ≈ 7 meV | %3,8 (n=4), a ≈ 17 nm, E_C/k_B ≈ 34 K, E₁ ≈ 13 meV |
| Bağlı durum tutan nokta, iyimser ölçüt | %7,5 (n=6), a ≈ 18 nm | %12,5 (n=11), a ≈ 15 nm |

Düzeltme notu: İlk hesapta ≈ 7 nm'lik klasik birikintileri nokta saymıştım; yanal kuantum düzeyi (≈ 80 meV) birikinti derinliğinden (15–45 meV)
büyük olduğundan bunlarda bağlı durum yoktur. Bağlı-durum ölçütüyle (E₁ < derinlik) oran %48'den %4–12'ye düştü.

**3B FEM kontrolü ve duyarlılık (`experiments/exp6_rdf/fem_check.py`):** Tek bir akseptör için gerçek 3B FEM (oksit 4,148 nm ε=3,9, Si ε=11,7,
topraklı kapı ve nötr düzlem, mesh yakınsaması <%1) görüntü modelimden **sistematik olarak büyük** potansiyel verdi: r=0'da derinlik 10 / 30 / 60 nm için
−%18 / −%7 / −%3,5 (görüntü modeli küçük kalıyor). FEM'e uydurulan derinlik bağımlı düzeltme (f = 1 + 1,8/z₀ · exp(−r²/2(z₀+5)²)) eklendiğinde:

| Varyant | σ_ψ | σ_VT (MC / Stolk) | bağlı nokta, katı | bağlı nokta, iyimser |
|---|---|---|---|---|
| ndep 1,7e17, görüntü | 20,6 meV | 3,9 / 5,7 mV | %3,7 (n=3) | %7,5 (n=6) |
| ndep 1,7e17, FEM düzeltmeli | 28,6 meV | 4,4 / 5,7 mV | %3,7 (n=3) | %7,5 (n=6) |
| N_eff 3,9e17, görüntü | 31,5 meV | 4,5 / 7,1 mV | %3,7 (n=4) | %12,5 (n=11) |
| N_eff 3,9e17, FEM düzeltmeli | 43,9 meV | 5,1 / 7,1 mV | %10 (n=8) | %21 (n=18) |

Düzeltme faktörü r=0 civarında uyduruldu, r büyüdükçe FEM görüntüden küçük kalıyor (r=40 nm'de +%6 ters yön), yani düzeltme yaklaşık. 80 örnek ve 3–18 nokta:
oranların belirsizliği büyük.

**[GERİ ÇEKİLDİ, bkz. kutu] Sınanabilir tahmin (Aşama 4 için):** 126 nm × 420 nm minimum aygıtların yaklaşık **%4–21'inde** (katı ölçüt %4–10, iyimser %7,5–21; dört varyantın birleşimi)
4 K'de eşik civarında Coulomb blokajı salınımı beklenir (E_C/k_B ≈ 25–40 K, yani E_C/kT ≈ 6–10 at 4 K).

**Güncelleme (exp9, `results/exp9/COMPARISON.md`):** Bu tahmin yukarıdaki gibi ayrı bir kontrol olmadan kullanılmamalı. Kıskaç duyarlılığı
σ_ψ'yi %23 değiştiriyor; modelde iyonlaşma T'den bağımsız, oysa NA ≪ N_c olduğu için 4 K'de katkılar tam iyonlaşmaz. Tahmin
"önceki" olarak okunmalı; yeni aralık katı ölçütte %2,5–10 (düzeltmeli üst sınır). Bu, tek bir aygıt değil **çok sayıda özdeş aygıt dizisi** ölçümü gerektirir (≥ 100 aygıt,
%4 ± 2 oran için). Tahmin yalnızca klasik birikinti + kaba bağlı-durum ölçütüne ve E_C sabitine (ε_eff = 7,8 ölçeklemesi) dayanır; tünel bariyer
dirençleri (R_T > h/e² koşulu) modellenmedi, dolayısıyla oran bir ÜST sınır olabilir.

## 3d. "Bölünme" (hücre doğurma) mekanizmalarının modeli (Aşama 2c, Deney 8)

Kullanıcının önerdiği üç aday mekanizma modellendi, hepsi analitik sonuçlarla doğrulandı; ayrıca kendi birleştirmem eklendi:
Coulomb-frustrasyonlu faz ayrışması (Ohta-Kawasaki tipi Cahn-Hilliard), Coulomb çekirdeği **geçit (gate) ekranlı**: ĝ(q) = A(1 − e^{−2qd})/q.

| Mekanizma | Model | Doğrulama | Sonuç |
|---|---|---|---|
| 1. Elektrostatik/Coulomb fisyonu | Lennard-Jones kümesi + gate-ekranlı Coulomb, tek küme vs iki parça | LJ7/13/19 global minimumları (hata < 5e-7) | Eşik yük g_c ∝ N^−1,08 (d=∞; Rayleigh/sıvı-damla beklentisi −1); gate yaklaştıkça g_c yükseliyor (d=4: eğim −0,73, d=2: −0,90) |
| 2. Faz ayrışması (Landau-Ginzburg) | Cahn-Hilliard + ĝ(q), yarı-örtük spektral | duvar gerilimi %0,41, σ(q) ≤ %0,32, kütle korunumu | Domainler doğuyor; boyut, A ve gate'e bağlı |
| 3. Simetri kırılması/bifurkasyon | Düzgün durumun lineer kararlılığı + genlik üssü | A_c analitik = sayısal (0,38490); genlik ∝ δ^0,502 | Süperkritik çatal bifurkasyonu (beklenen 0,5); A_c(d) gate yaklaştıkça artıyor (d=1,5: 0,4846) |

Bulgular (ayrıntı: `results/exp8/summary_exp8.md`):
- **Gate (kapasitör) bölünmeyi kontrol ediyor.** φ̄ = 0'da desen dalga boyu d=1,5'te d=∞'a göre 1,5–1,7× büyük; d ≳ 12'de etki yok. φ̄ = −0,4, A = 0,15'te düzgün
  durum d ≥ 6'da kararlı (desen yok), d ≤ 3'te kararsız: aynı frustrasyonda gate'i yaklaştırmak bölünmeyi başlatıyor.
- **Negatif sonuç:** Hazır tek bir domain (disk, çubuk, elips) Coulomb frustrasyonunda **bölünmedi**; küçüldü, A=0,05'te çözündü (`division_diag.json`).
  Bu modelde bölünme, tek hücrenin şekil kararsızlığı değil, düzgün durumun bifurkasyonla çok sayıda hücreye ayrılmasıdır.

**Güncelleme (exp10, `results/exp10/summary_exp10.md`):** φ ≈ (n − n₀)/n₀ olarak SKY130 inversiyon tabakasına bağlandı.
Değişim (Hartree–Fock, kapı ekranlı) + korelasyon (Attaccalite 2002) ile yerel ters kompresibilite n_s ≲ 1,5–3·10¹² cm⁻²'de negatif:
çift kuyu fiziksel olarak var. Ama kapı ekranlı Hartree k_F'ye kadar |K|'nın 11–21 katı; düzgün tabaka kararlı, SKY130 çalışma
aralığında bölünme yok. Aşağıdaki paragraf bu güncellemeden öncedir.
**İnceleme (2026-10-11):** exp10'un "kapı çok yaklaşınca bölünme açılıyor" haritası, kapı ekranlı olmayan korelasyon vekilinin ürünü.
Tutarlı kapı ekranlı HF'de tabaka hiçbir kapıda (0,05–3000 µF/cm²) ayrışmıyor. Yukarıdaki "gate bölünmeyi başlatıyor" bulgusu soyut
modelde geçerli, silikon karşılığı yok.

Cihaza çeviri ve açık noktalar: Model birimsiz; φ'nin fiziksel karşılığı (2B elektron sıvısı faz ayrışması, ferroelektrik polarizasyon vb.) bu çalışmada bağlanmadı.
Düşük yoğunluklu 2B elektron gazında negatif sıkıştırılabilirlik kaynaklı faz ayrışması literatürde tartışılır (bellekten, doğrulanmadı); SKY130 evirme tabakasında
böyle bir kararsızlık olup olmadığı için değişim-korelasyon içeren bir model gerekir (1B Poisson-Schrödinger çözücüm bunu içermiyor).

## 3e. Sürülen (disipatif) bölünme: süperkapasitör fikri ile disipatif yapı karşılaştırması (Deney 11)

Kullanıcının itirazı üzerine (Prigogine, disipatif yapılar) iki model aynı Cahn–Hilliard çekirdeğinde, aynı çözücüyle kıyaslandı:
A = süperkapasitör fikri (yüklü hücre + kapı ekranlı Coulomb + elektrikle şarj), B = disipatif yapı (Zwicker 2017 aktif damla + Gray–Scott).
Ayrıntı: `results/exp11/summary_exp11.md`.

Sürüm 2 (incelemeden sonra; `results/exp11/review_exp6_exp10_exp11.md`):

| Sonuç | A: süperkapasitör | B1: aktif damla | B2: Gray–Scott | C: kazara nokta (2b) |
|---|---|---|---|---|
| 2B kalıcı bölünme | 0/28 (güçlü şarjda 12/12 ortamdan doğum) | 2/30 (ikili 0; 9 solucan) | 5/42 çoğalıyor, 6/6 dayanıklı | bölünme yok; tahmin geri çekildi |
| 3B | ekransız Coulomb ≡ Δν = 0: bölünür | Δν −0,002…+0,004 hepsi bölünür | — | — |
| Silikon | tutarlı HF'de ayrışma yok | 2B'de filaman | en iyimser 123 nm (E_g), tipik ≈ 0,4 µm | a ≈ 35–41 nm, ekransız |

Okuma: Tek kazanan yok. Bölünmeyi boyut (3B kopar, 2B kopmaz), uzun menzilli itme (kapı ekranlar) ve otokatalitik
aktivatör–inhibitör yapısı belirliyor. "Disipatif" etiketi tek başına belirleyici değil: 3B'de denge eşdeğeri Δν = 0 da bölünüyor.
MOSFET kanalında (2B, kapılı) yalnız Gray–Scott türü otokatalitik sistemler aday. 4 K kuantum hedefine (R ≈ 6–20 nm) hiçbiri inmiyor.
Silikonda gerçek disipatif aday: 4 K bulk CMOS'ta darbe iyonlaşması + yüzer (donmuş) gövde. Bu, Exp12'nin konusu.

## 3f. 4 K yüzer gövde kilitlenmesi ve akım filamanı (Deney 12)

Araç: ngspice 44.2 ve resmî SKY130 BSIM4 kartı (300 K; kanal akımı, darbe iyonlaşması, V_th tabloları). 4 K etkileri açık düğmeler.
Ayrıntı: `results/exp12/summary_exp12.md`.

| Sonuç | Değer |
|---|---|
| Eşik üstü (V_gs 1,2–1,8) | geri besleme ×1,3–1,5; yüzer gövdeyle kink var, histerezis yok |
| Eşik altı | yüzer gövde kilitlenmesi (SOI "tek transistör kilitlenmesi" benzeri); kink + histerezis düğmelerin %24–97'sinde |
| Kanal genişliği boyunca (W = 100 µm, yük direnci) | donmuş kuyu R_□ ≳ 10¹⁰–10¹² Ω/□ iken tek akım filamanı; V_ds Maxwell geriliminde (0,80 V, analitik = simülasyon) kenetli; filaman akımla genişliyor, bölünmüyor |
| Ölçüm önerisi | geniş nFET, V_gs ≈ V_T − 0,1 V, seri dirençli drenaj: sabit-V_ds bölümü (aralık ∝ W) ve histerezis |

Okuma: SKY130'da 4 K'de gerçekçi bir disipatif desen adayı var, ama çoğalan ya da küçülen hücre değil; genişleyen tek bir klasik filaman.
Kuantum hedefine katkısı yok. Ölçülürse donmuş kuyunun bilinmeyen parametrelerini verir.

## 3g. Yerel inhibitörlü filaman: bölünme (Deney 13)

Dış inceleme (`Desktop/inceleme.txt`) "eşik altında çok hücreli model yapılmadı" diyordu. Bu doğru değildi: Exp12'nin uzaysal modeli
zaten buydu. Ama inceleme gerçek bir boşluğa işaret ediyordu: yerel inhibitör test edilmemişti. Exp13'te dirençli drenaj eklendi.
Sonuçlar:
- Yalnız global bağlaşımla 5/5 parametre setinde tek filaman.
- Uçtan beslenen hatta filamanlar kontaklarda doğuyor.
- Dağıtık balastta (λ = 40 µm) filaman ≈ 4λ'da ikiye bölünüyor (3/3 tohum; mekanizma: merkezde drenaj gerilimi en düşük).
- Tekrarlayan çoğalma ve küçülme yok; kızlar 70–100 µm, klasik.

Tasarım önerisi: bölümlenmiş drenaj + poly balast, bölüm akımları ayrı ölçülür. Ayrıntı: `results/exp13/summary_exp13.md`.

## 3h. Boyut ve bölünme (Deney 14–16)

Kullanıcının isteği iki parçalıydı: "6B hiperküpte Cahn–Hilliard; bölünmenin boyuta bağlı olduğunu göster" ve "fiziksel 6B faz uzayı (Vlasov–Poisson) ile veriden denklem".
- **Exp14 (aktif CH, d sürekli):** bölünme bir boyut penceresinde oluyor. Pencerenin altında 2B şerit, üstünde içi boş sıvı kabuk var.
  - Pencere iç tüketim Δν ile yukarı kayıyor; 6B ancak Δν ≥ +0,012'de bölünüyor.
  - Gerçek 6B hiperküp (16⁶) O(5) indirgemeyi ≤ %2,5 farkla doğruladı.
  - Tam 6B ızgarada bölünme simülasyonu imkânsız: ≥ 10¹¹ nokta gerekir. Pasif tüpler küçük hiperküpte çözülüyor ("boyut laneti").
- **Exp15 (PIC):** plazma delikleri 1B'de birleşiyor, 2B/3B'de parçalanıp sönüyor. Eşik 1B↔2B, bölünme yok.
- **Exp16 (PySINDy):** bilinen yasalar geri bulundu, gürültüye dayanıksız. Exp8 ve Exp13 yorumlarına iki düzeltme getirdi.
- **Exp11'e düzeltmeler:** 3B Δν = −0,004 yapısı bir kabuk. 3B Δν = +0,004 bölünmesi 48'lik kutudan etkilenmiş.

Ayrıntı: `results/exp14_dimension/`, `results/exp15_phase_space/`, `results/exp16_equations/` içindeki `summary_*.md`.

## 3i. Silikonda bölünen filamanlar: olgun denklem (Deney 17)

Hedef: "verilerden olgun bir denklem ve silikonda bölünebilen şekiller".
- **Denklem:** Exp12/13 modelinden türetildi; bölünme, filaman merkezindeki drenaj gerilimi V_down'a indiğinde olur.
  - Tek filaman: w_s = 2λ ln[λ/(λ−λ_c)], λ_c = √(2(V_M−V_down)/(r_D j_on)).
  - Dizi için kapalı form (bu çalışmada sınanmadı); w_top = I/j_on. Filaman sayısı I/(j_on w_s) ile basamak basamak artıyor
    (0,75'lik çizgi kılavuzdur, uydurma değil).
- **Sınama:** öngörüler önceden kaydedildi. Temiz protokoldeki 22 bölünmede oran 1,00 ± 0,07 (63 olayın grup tablosu Exp17 raporunun düzeltme kutusunda). Üç cihaz setinin 14 temiz ilk bölünmesi tek evrensel eğride (RMS %9). λ = 36 µm (λ_c'nin %1 üstü) ıskalandı: öngörülen 322 µm'de bölünme olmadı. Veriden tek katsayı uydurulunca kuramın 1,01 katı.
- **Şekil:** filaman akımla büyüyüp tekrar tekrar ikiye bölünüyor (3–4 kuşak). Düzgün cihazda Turing benzeri eşit aralıklı dizi oluşuyor.
- **SKY130 test yapısı:** li1 drenaj hattı + poly balast, λ taraması. Sabitler tek bir referans cihazda ölçülebilir.
- **Geçerlilik koşulu:** yarı statik sürüş.

Ayrıntı: `results/exp17/summary_exp17.md`.

## 4. Literatür (arama sonuçları, doğrulama durumu ile)

| Bulgu | Kaynak | Durum |
|---|---|---|
| SKY130: poly genişliği 0,150 µm, poly-poly aralığı 0,210 µm | [SkyWater PDK kural sayfası](https://skywater-pdk.readthedocs.io/en/main/rules/periphery.html) | birincil kaynaktan okundu |
| Ticari bulk 28 nm PMOS'ta 14 K'de Coulomb blokajı tümseği | [arXiv 2206.09086](https://arxiv.org/pdf/2206.09086) | arama özeti, tam metin okunmadı |
| Ticari bulk 40 nm MOSFET'lerde 50 mK'de kanalda kuantum nokta ve Coulomb blokajı | [arXiv 2005.13863](https://arxiv.org/pdf/2005.13863) | arama özeti, tam metin okunmadı |
| Tek elektron tuzaklama, minimum boy nMOS, 50 mK, çok kanallı okuma yongası | [arXiv 2101.08295](https://arxiv.org/pdf/2101.08295) | düğüm belirtilmemiş |
| 22FDX ve 28FDSOI'de endüstriyel kuantum nokta ve qubit çalışmaları | [arXiv 2412.08422](https://arxiv.org/pdf/2412.08422), [arXiv 2504.20572](https://arxiv.org/pdf/2504.20572) | arama özeti |
| Oda sıcaklığında Coulomb elmasları sub-5 nm Si adasında (araştırma cihazı, standart CMOS değil) | [Shin vd., Nano Lett. 2011 (Korea Univ.)](https://pure.korea.ac.kr/en/publications/room-temperature-charge-stability-modulated-by-quantum-effects-in/) | arama özeti |
| Si derin safsızlık qubiti, TFET'te, 5–10 K'ye kadar çalışma | [arXiv 1804.03364](https://arxiv.org/abs/1804.03364) | özet okundu; T1/T2 özette yok |
| SKY130 kriyojenik ölçümleri: 77 K model ve NIST veri seti | [NIST veri seti](https://data.nist.gov/od/id/mds2-2997) | **sıcaklık noktaları doğrulanamadı** (sayfa içeriği yüklenmedi); 4 K verisi var mı bilinmiyor |

Aramada **130/180 nm düğümünde** kuantum nokta / Coulomb blokajı raporu bulamadım. Bu "yok" demek değil, "bu aramada bulunmadı" demektir.
Arama motoru bir kaynakta "≤ 80 nm cihazlarda 300 K'de iletkenlik tepesi" dedi; özetten doğrulayamadığım için kullanmadım.

## 5. Aşamalı yol haritası (her aşamada geçiş ölçütü)

| Aşama | İş | Çıktı / geçiş ölçütü | Durum |
|---|---|---|---|
| 0 | Elektrostatik, hapsolma, zincir, çift nokta simülasyonu | doğrulanmış modeller | **tamam** |
| 1 | Tünel bağlaşımı (WKB doğrulamalı), dinamik yükleme (LZ doğrulamalı), ölçek boşluğu | bu belge | **tamam** |
| 2a | SKY130 yığınında dikey 2B elektron gazı: 1B öz-uyumlu Poisson-Schrödinger (parametreler resmî model dosyasından, §3b) | alt bant enerjileri, tabaka kalınlığı, V_T tutarlılığı | **tamam** (1B) |
| 2b | Rastgele katkı potansiyeli ve kazara nokta oluşumu, 126 nm × 420 nm kanal (birinci mertebe model, §3c) | kazara nokta oranı, boyutu, E_C | **yeniden yapılmalı**: 2026-10-11 incelemesinde oran tahmini geri çekildi (kenar artefaktı, aşırı FEM düzeltmesi, ekranlama yok). Gerekenler: öz-uyumlu Thomas–Fermi ekranlama, S/D + STI geometrisi, tam q-uzayı çekirdeği, yüzer gövde (4 K) |
| 2c | "Bölünme" (hücre doğurma) mekanizmalarının modeli: Coulomb fisyonu, Cahn-Hilliard + gate-ekranlı Coulomb, bifurkasyon (§3d) | eşik, boyut seçimi, gate bağımlılığı | **tamam** (soyut, birimsiz model; fiziksel malzemeye bağlanmadı) |
| 2d | Sürülen (disipatif) bölünme: süperkapasitör (A) ile disipatif yapı (B) kıyası, 3B aktif damla doğrulaması, küçülerek çoğalma rampası (§3e) | hangi sürüşte, hangi bantta bölünme; silikon ölçeği | **tamam** (sürüm 2, incelemeden sonra: belirleyici olan boyut, menzil ve otokataliz; "disipatif" tek başına değil); 4 K akım filamanı modeli Exp12'de |
| 2e | 4 K yüzer gövde kilitlenmesi ve akım filamanı (Exp12, §3f) | eşik altı kilitlenme, filaman koşulu (W ≫ ℓ), Maxwell gerilimi | **tamam** (model; 4 K parametreleri varsayım; ölçüm önerisi hazır) |
| 2f | Yerel inhibitörlü filaman bölünmesi (Exp13, §3g) | çoklu filaman, bölünme mi doğum mu | **tamam** (model): dağıtık balastta tek bölünme, tekrarlayan çoğalma yok |
| 2h | Silikonda tekrarlayan filaman bölünmesi ve olgun bölünme denklemi (Exp17, §3i) | önceden kaydedilmiş öngörüyle sınama, evrensellik, mekanizma | **tamam** (model): temiz protokolün 22 bölünmesinde merkez–kenar farkı = V_M − V_down (1,00 ± 0,07); ölçüm önerisi hazır |
| 2g | Boyuta bağlı bölünme: d = 2–6 aktif CH + 6B hiperküp (Exp14), faz uzayı PIC (Exp15), veriden denklem (Exp16) (§3h) | bölünme penceresi, d_c, doğrulanmış indirgeme | **tamam** (model): pencere 2 < d ≲ 3,5–4, Δν ile kayıyor; 6B güçlü sürüşte bölünüyor; plazmada bölünme yok |
| 3 | SKY130 kriyojenik veri: NIST veri setini indirip sıcaklık noktalarını doğrula | ≥ bir 4 K veya daha düşük veri var mı? | **engellendi**: veri seti kaydı (API) okundu, dosya sunucusu (ETest Tile.zip 3,8 MB, MPW-5 Test Tile.zip 39 MB, README) Cloudflare 524 hatası verdi; sonra yeniden denenecek |
| 4 | Test yapısı tasarımı: 150 nm × 420 nm minimum aygıtlar, ince şeritler, çok kapılı birkaç yapı; açık MPW ile üretim, ≤ 1 K'de ölçüm | Coulomb blokajı salınımı var mı? E_C ne? | donanım ve bütçe gerekir; erişilebilirlik ve fiyat doğrulanmalı |
| 5 | Çift nokta / tünel bağlaşımı yalnızca Aşama 4 olumluysa | tünel yarılması ölçümü | koşullu |
| 6 | Oda sıcaklığı | sub-5 nm adalar: standart CMOS'ta mümkün değil; post-işleme (oksitleme ile daraltma, e-ışın) gerekir | bu projenin kapsamı dışında; "no-go" (standart süreçte) |

Geçiş ölçütleri:
- Aşama 3'te 4 K altı SKY130 verisi yoksa, Aşama 4'ün ilk ölçümü zaten bu veriyi üretir.
- Aşama 4'te Coulomb salınımı **görülmezse** (ör. E_C/k_B < 1 K), 130 nm hattı kuantum için kapatılır, SKY130 yalnızca klasik kontrol/okuma elektroniği olarak kalır.
- Aşama 5 yalnızca Aşama 4'te tutarlı (tekrarlanabilir) tek elektron davranışı varsa başlar.

## 6. Hâlâ bilinmeyenler

- Elektrikle "hücre doğurma" mekanizması: Exp11'de gösterildi. 2B'de Gray–Scott türü otokatalitik sistemlerde; 3B'de aktif damlada ve ekransız Coulomb'lu yüklü damlada. Ölçek değişiminde yeni lekeler bölünmeyle doğuyor. Silikonda karşılığı (4 K darbe iyonlaşması + yüzer gövde) Exp12'de; hücreler klasik olmalı (mertebe tahmini: en iyimser yaklaşık 0,1 µm).
- SKY130 minimum aygıtlarında kazara nokta oluşur mu ve E_C kaç olur? (Aşama 2b tahmini geri çekildi; öz-uyumlu modelle yeniden yapılmalı. Projede tek bir düzlemsel çukur kapasitansı tanımı gerekli: aynı 20 nm nokta için raporlar arasında E_C 3,7 kat farklıydı.)
- 4 K altı SKY130 aygıt davranışı (NIST veri setinin içeriği doğrulanamadı).
- Optik fonon, Raman, yük gürültüsü: 300 K koherans üst sınırları gerçek T2'yi aşar; τφ senaryoları fenomenolojik.
- Landau-Zener süpürme modeli: dephasing banyo ile tutarlı değil (uzun zaman sonuçları güvenilmez), F(q) yalnızca antikesişim enerjisinde hesaplandı.
- Elektron-fonon bağlaşım sabiti Ξ=9 eV ve SiO₂ efektif kütlesi 0,5 m₀ varsayım.

## 7. Verilmesi gereken karar

1. Hedefi, "SKY130'da ≤ 1 K tek elektron davranışı arayışı" olarak yeniden tanımlıyor musun? (Aşama 2-4)
2. Aşama 4 için donanım ve bütçe (MPW, ölçüm kriyostatı) mevcut mu? Yoksa yalnızca simülasyon + literatür kapsamında mı kalalım?
3. Oda sıcaklığı hedefi, standart CMOS dışı işlem (post-işleme) kapsamına bırakılsın mı?
4. ~~Exp12 yapılsın mı?~~ Kullanıcı kabul etti (2026-10-11); Exp12 başladı.
