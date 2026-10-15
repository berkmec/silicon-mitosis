# Deney 2 ve 3 — çoğalma zinciri, çift kuantum nokta ve oda sıcaklığı qubit taraması

## Önce sınırlar

- Oda sıcaklığında kuantum durumunu ve qubiti **zorlayamam**: bu belge hangi parametrelerde neyin mümkün olduğunu gösterir.
- Fonon hesabı yalnızca **tek-fonon akustik** kanaldır (T2 için ÜST sınır, T2 = 2·T1). Optik fonon, iki-fonon (Raman), vadi/orbital ve yük gürültüsü hesaba katılmadı. Bu yüzden 300 K'de dephasing için τφ = 10 fs / 100 fs / 1 ps **senaryoları** ayrıca verildi; bunlar tahmin değildir.
- Elektron-fonon bağlaşımı: Si boyuna akustik deformasyon potansiyeli Ξ = 9 eV (varsayım; literatür aralığı yaklaşık 5-10 eV).

## Adım A — kademeli küçülen zincir (12 / 6 / 3 nm, ε=3,9)

- FEM kapasitans matrisi: satır toplamı hatası ≤ 5e-13, simetri hatası ≤ 1.2e-14; mesh seviye 1→2 adacık bloğu en büyük göreli değişim %4.1 (tam yakınsamamış, üst sınır).
- Adaların E_C/k_B değerleri: 120 K, 232 K, 404 K.
- Yük kuantizasyon indeksi (1 = tam basamak, 0 = sürekli) ve detrend edilmiş akım salınımı:

| T (K) | salınım genliği/ort. akım | kuantizasyon (ada 12/6/3 nm) |
|---|---|---|
| 4 | 6.0748 | 0.93 / 0.98 / 1.00 |
| 77 | 0.0179 | -0.00 / -0.06 / 0.58 |
| 300 | 0.0010 | 0.00 / -0.18 / -0.29 |

Okuma: 4 K'de ada 0 ve ada 1 yük basamakları gösteriyor ama en küçük ada (3 nm) hiç yüklenmiyor; 77 K ve 300 K'de yük sürekli (kuantizasyon ≈ 0), akım salınımı yok. 77 K'de ada 2 için görünen 0,58, doluluğun zaten ≈ 0 olmasından kaynaklanıyor, kuantizasyon değil.

## Yakınsama (çift nokta, L=1 nm, d=0,5 nm)

| h_int (nm) | kareleme derecesi | dof | Δ (meV) | <|F|²> |
|---|---|---|---|---|
| 0.22 | order6/sub1 | 7145 | 49.93 | 7.568e-05 |
| 0.16 | order4/sub0 | 19682 | 49.71 | 7.961e-05 |
| 0.16 | order9/sub0 | 19682 | 49.71 | 7.960e-05 |
| 0.16 | order6/sub1 | 19682 | 49.71 | 7.960e-05 |
| 0.12 | order6/sub1 | 43610 | 49.63 | 8.110e-05 |

## Adım B/C — çift kuantum nokta tablosu (300 K)

Kapı doğruluğu: ortalama X-kapı doğruluğu; 0,5 = kanal tamamen depolarize (kapı ölü). THz = Ω=0,05·Δ/ħ (RWA güvenli), hızlı = Ω=0,5·Δ/ħ (iyimser, RWA dışı), elektriksel = 10 ps darbe. τφ = 100 fs senaryosu.

| L (nm) | d (nm) | Δ (meV) | P(300K) | q0 (1/nm) | <|F|²> | T1 akustik (s) | Q | F_X THz | F_X hızlı | F_X elektriksel | F=0,99 için gereken T2 (s) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1.0 | 0.5 | 49.6 | 0.74 | 8.9 | 8.1e-05 | 6.8e-13 | 7.0 | 0.549 | 0.784 | 0.500 | — |
| 1.0 | 0.8 | 11.3 | 0.21 | 2.0 | 4.0e-01 | 3.4e-15 | 0.1 | 0.500 | 0.500 | 0.500 | — |
| 1.0 | 1.2 | 1.6 | 0.03 | 0.3 | 3.2e-02 | 2.1e-12 | 0.2 | 0.500 | 0.549 | 0.502 | — |
| 1.5 | 0.5 | 18.9 | 0.35 | 3.4 | 8.1e-02 | 5.7e-15 | 0.3 | 0.500 | 0.500 | 0.500 | — |
| 1.5 | 0.8 | 3.5 | 0.07 | 0.6 | 1.5e-01 | 9.3e-14 | 0.3 | 0.500 | 0.500 | 0.500 | — |
| 1.5 | 1.2 | 0.4 | 0.01 | 0.1 | 2.8e-03 | 4.4e-10 | 0.1 | 0.630 | 0.663 | 0.663 | 6.2e-10 |
| 2.0 | 0.5 | 9.3 | 0.18 | 1.7 | 3.2e-01 | 6.2e-15 | 0.2 | 0.500 | 0.500 | 0.500 | — |
| 2.0 | 0.8 | 1.6 | 0.03 | 0.3 | 5.1e-02 | 1.4e-12 | 0.2 | 0.500 | 0.524 | 0.500 | — |
| 2.0 | 1.2 | 0.2 | 0.00 | 0.0 | 6.3e-04 | 1.2e-08 | 0.0 | 0.663 | 0.666 | 0.667 | 9.3e-10 |

Aer tutarlılık kontrolü (ilk konfigürasyon): popülasyon |1> Aer = 0.382262, kanal = 0.382262.

## Sonuç (dürüst okuma)

- 300 K'de başlangıç polarizasyonu ≥ 0,9 (Δ ≥ 76 meV) olan konfigürasyon: 0 / 9. En iyisi (L=1.0 nm, d=0.5 nm): Δ=50 meV, P=0.74, T1 (akustik) = 6.8e-13 s.
- Bu en iyi konfigürasyonda bile (τφ=100 fs, iyimser hızlı kapı) X-kapı doğruluğu 0.78; τφ=1 ps olsaydı 0.94. %99 doğruluğa, T2 üst sınırı (2·T1) verilse bile ulaşılamıyor (akustik T1 ve termal popülasyon sınırı).
Tablodaki '—' = bu konfigürasyonda %99'a T2 = 2·T1 ile bile ulaşılamıyor.
- **Optimum taramanın kenarında** (en küçük L, en ince d). Δ bu yönde artıyor; ancak 1 nm altında sürekli ortam ve etkin kütle modeli geçersiz, bu yüzden daha küçüğüne bu modelle çıkarım yapılamaz (atomistik yöntem gerekir).
- Fonon hesabı yalnızca tek-fonon akustik kanaldır; optik fonon, Raman, yük gürültüsü hesaplanmadı. T1 değerleri konfigürasyona göre fs'ten µs'ye kadar çok geniş değişiyor; Δ büyüdükçe fonon dalga boyu nokta boyutundan küçülüp bağlaşım bastırılıyor (fonon darboğazı), Δ ≈ 10 meV civarında bağlaşım en güçlü (λ ≈ nokta boyutu).
- Elektriksel 10 ps kapı darbesi bütün konfigürasyonlarda koherans süresinden çok uzun; bu modelde elektrikle kontrol edilen oda sıcaklığı silikon qubiti gösterilemedi. Gösterilebilen tek şey: en küçük/ince geometride koherans için bir *fs* penceresi (Q≈7), ve bunun kullanılabilmesi için fs mertebesinde THz/optik kontrol ile τφ ≳ ps gerektiği.