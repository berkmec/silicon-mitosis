# Deney 2 ve 3 varsayımları

- Zincir: S - 12 nm - 6 nm - 3 nm - D, ε=3,9, ideal iletken adalar, ortak kapı kutusu (kenar 3,4·a0). FEM kapasitans (209 k düğüm) +
  ortodoks ana denklem (nmax=6, 343 durum), R_T = 1 MΩ. Seviye 1→2 mesh değişimi en çok %4,1 (tam yakınsamamış).
- Çift nokta: iki küp kuyu (kenar L), bariyer kalınlığı d, SiO2 kutu (pay 0,6 nm), V0=3,1 eV, m_Si=0,26, m_SiO2=0,5 (varsayım).
- Fonon: yalnız boyuna akustik, deformasyon potansiyeli Ξ=9 eV (varsayım, literatür ~5-10 eV), ρ=2329 kg/m³, v=8430 m/s.
  Γ_em = Ξ² q0³ <|F|²> (n+1) / (2π ħ ρ v²); T1 = 1/(Γ0 (2n+1)); T2 ≤ 2·T1 (yalnız tek-fonon, ÜST sınır).
- F(q) FEM dalga fonksiyonlarından, alt bölmeli (8 alt dörtyüzlü, 15 noktalı) kareleme ile; kareleme derecesi etkisiz, mesh yakınsaması ~%2.
- Qiskit Aer: thermal_relaxation_error ile X kapısı; ortalama kapı doğruluğu kanal olarak hesaplandı, Aer yoğunluk matrisi ile tutarlı.

# Deney 1 (FEM) varsayımları

- Araçlar: scikit-fem 12.0.2 (P1 Laplace, P2 Schrödinger), gmsh 4.15.2, dolfinx 0.9.0 (çapraz kontrol, WSL).
- Elektrostatik: ada/kaynak/savak/kapı ideal iletken; dolgu tek ε (3,9 veya 11,7); Maxwell kapasitans matrisi C_ij = φ_i^T K φ_j.
- Geometri: ada kenarı a, tünel aralığı 0,1·a, kaynak/savak küpleri 0,3·a, kapı kabuğu boşluğu b = a/r.
- Doğrulama: eş merkezli küreler (analitik), dolfinx ile aynı meshte 1e-15 uyum, satır toplamı=0 ve simetri testleri.
- Schrödinger: tek parçacık, BenDaniel-Duke, m_Si=0,26, m_SiO2=0,5 (varsayım), V0=3,1 eV, bariyer payı 0,8 nm.
- Taşıma: FEM'den gelen C_ij üstünde ortodoks (sabit etkileşim) teori; ilk ilkelerden değil.
- Bilinen sınırlar: kuantum kapasitans yok, <1 nm geçersiz, çok-elektron yok, çoğalma dinamiği yok.

# Model varsayımları (Deney 0, analitik ön hesap)

- Hücre: kenarı L = L0·r^g; öz-kapasitans küre (yarıçap L/2, εr=11.7); kaynak/savak kapasitansı öz-kapasitansın yarısı.
- Kapı kapasitansı: paralel plaka, kalınlık = 0.1·L, εox=3.9, "boost" çarpanı süperkapasitörü temsil eder (kaba).
- Hapsolma: kübik kutu, m*=0.26 m0, ΔE = 3ħ²π²/(2m*L²).
- Eşik: enerji ≥ 10·kB·T.
- Taşıma: ortodoks teori, sabit etkileşim modeli, tek nokta, R_T=1 MΩ (>> h/e²).
- Çoğalma kuralı: ½·C_g·V² ≥ E_C(çocuk). Kırılma: 1e9 V/m.
- Bilinen kaba noktalar: kapı kapasitansı çok küçük hücrelerde fazla tahmin edilir; <2 nm'de sürekli ortam yaklaşımı geçersiz;
  kuantum kapasitans, çok-seviyeli taşıma, kotünelleme ve dekoherans yok.
