# Deney 13: Yerel inhibitörlü (dirençli drenaj) SKY130 nFET — 4 K eşik altı kilitlenmede çoklu filaman ve bölünme

Tarih: 2026-10-11. Kullanıcının paylaştığı dış inceleme (`Desktop/inceleme.txt`) üzerine yapıldı.
- Kod: `experiments/exp13_rail/` (run, figures) ve `experiments/exp12_filament/spatial.py` (drenaj hattı ve soy takibi eklendi).
- Testler: `tests/test_exp13.py`.
- Şekil: `fig_exp13.png`.
- Veri: `A*/B*/C*/D*.json`, `*_fine.json`, `split_robustness.json`.

## 0. İnceleme metnine cevap

| İddia | Değerlendirme |
|---|---|
| "Exp12 yakınsamıyor (dy 0,5 / 0,25 / 0,1 → %47 / %86 / %78)" | Satırlar karışmış. %47 R_□ = 10¹²'nin satırı; %86 ve %78 R_□ = 10¹⁰'un. "Filaman genişliyor" sonucunun dayandığı R_□ = 10¹²'de üç ızgara %47 / %48 / %47 veriyor, yani yakınsamış. R_□ = 10¹⁰ burada dy = 0,1 ve 0,05 ile yeniden koşuldu; ikisi aynı sonucu (tek filaman) veriyor, kaba ızgara artefakttı. Ayrıca sütun V_ds değil V_DD |
| "Eşik altında çok hücreli model yapılmadı" | Yanlış. Exp12'nin uzaysal modeli tam olarak bu: 200–1000 hücre, hücre başına gövde potansiyeli, V_ov0 = −0,1 V (V_gs ≈ V_T − 0,1 V), yük direnci, SKY130 eşleşmesi. Exp12 özetinde bunu açık yazmadığım için yanlış okunmaya açıktı; düzeltildi |
| "Otokataliz varsa 2B'de de bölünür" | Eksik. Gray–Scott tipi bölünme için otokatalizin yanında daha geniş menzilli yerel bir inhibitör gerekir. Exp12'deki tek inhibitör global yük doğrusuydu, bu da teoride tek domen verir; sonuç da buydu. **Gerçek boşluk buydu ve Exp13'te test edildi** |
| A3'te 31 hücre, N_min = 10 ile 31–69 nm | sayılar doğru. Yorumu Exp11 raporunda: ortamdan doğum, Coulomb'dan değil; boyut tek başına kuantum için yetmez |

## 1. Model

Exp12'nin 4 K indirgenmiş modeline yerel bir inhibitör eklendi. Her hücrenin drenaj gerilimi artık yerel; drenaj hattının akım
dengesinden (KCL) hesaplanıyor. Darbe iyonlaşması V_ds'ye çok duyarlı (M: 1,2 V'ta 6·10⁻⁶, 2,4 V'ta 7·10⁻³), bu yüzden iletken bir
filaman komşularının V_ds'ini düşürüp onları bastırıyor. İki geometri kullanıldı:

- **Uçlardan beslenen hat:** kontaklar iki uçta, hat direnci r_D (Ω/µm).
- **Dağıtık balast:** her hücre ortak düğüme kendi balast direnciyle (r_v, Ω·µm) bağlı. İnhibisyon menzili λ = √(r_v / r_D) ve konumdan
  bağımsız.

Doğrulama (39/39 test):
- iki uçtan beslenen hatta merkez düşümü analitik r_D·I·W/8 ile birebir;
- dağıtık balastta sönüm boyu geniş alanda 20,000 µm (beklenen 20; dar alanda cosh etkisiyle 23,6 µm çıkmıştı);
- r_D = 0'da Exp12 aynen korunuyor.

## 2. Sonuçlar

| Koşu | Sonuç |
|---|---|
| A: global bağlaşım, 5 farklı bistabil parametre seti (W = 100 µm) | 5/5'inde tek filaman (doğum 1, bölünme 0). Exp12 sonucu tek bir parametre setine bağlı değil |
| D: R_□ = 10¹⁰, dy = 0,1 ve 0,05 µm | ikisi aynı: tek filaman. Exp12'deki yakınsamama kaba ızgaradan geliyordu |
| B: uçlardan beslenen hat (r_D 30–1000 Ω/µm; W 100 ve 300 µm) | filamanlar kontaklarda **doğuyor** (N_max = 2), konumları tasarımdan geliyor. Yukarı taramada bölünme yok. W = 100'de aşağı taramada tam kaplamış domen **ortadan ikiye** ayrılıyor (merkezde hat düşümü en büyük) |
| C: dağıtık balast (W = 300 µm, λ = 5 / 10 / 20 / 40 µm; yük 0 ve 5 kΩ) | çoğu durumda tek domen. λ = 5, yük 0'daki "3 filaman" ince sürüş adımında kayboldu (adım artefaktı). **λ = 40 µm, 5 kΩ: büyüme → ikiye bölünme** |
| C λ = 40 µm sağlamlık | 3 tohumun 3'ünde de bölünme (V_DD = 2,79 / 3,02 / 3,04 V). Filaman toplam ≈ 150–170 µm'ye (≈ 4λ) ulaşınca ortası sönüyor; kızlar 72–98 µm. Mekanizma: yerel drenaj gerilimi filamanın ortasında en düşük (≈ 0,785 V), darbe iyonlaşması orada bastırılıyor. Gray–Scott'taki "merkezi tükenen leke" ile aynı |
| C λ = 40 µm, W = 600 µm (10 kΩ) | V_DD ≤ 3,4 V'ta bölünme yok; sürüş aralığı kritik genişliğe yetmiyor olabilir (belirsiz) |

## 3. Okuma

- **Bölünme, otokatalize ek olarak yerel bir inhibitör gerektiriyor.** Yalnız global bağlaşımla tek filaman genişliyor; yerel
  inhibitörle (dağıtık balast, λ ≈ 40 µm) filaman büyüyüp ikiye bölünüyor. Dış incelemenin sezgisi bu ek koşulla doğrulandı.
- **Tekrarlayan çoğalma görülmedi.** Taranan aralıkta en fazla bir bölünme oldu; kızlar yeniden bölünmedi.
- **Kızlar büyük ve klasik.** Bölünme ≈ 4λ'da oluyor ve kızlar 70–100 µm genişliğinde. "Küçülerek çoğalma" yok; kuantum ilgisi yok.
  Hücre boyunu λ belirliyor. λ küçültülürse (5–20 µm) bu modelde tek domen kaldı, bölünme görülmedi.

## 4. Tasarım önerisi (SKY130, 4 K)

- **Test yapısı:** bölümlenmiş drenajlı geniş nFET. W ≈ 300 µm, L = 0,15 µm.
- **Drenaj hattı:** silisitli şerit, r_D ≈ 10 Ω/µm.
- **Balast:** her bölüme poly direnç; λ = 40 µm için yaklaşık 16 kΩ·µm, örneğin her 10 µm'lik bölüme 1,6 kΩ.
- **Sürüş:** V_gs ≈ V_T − 0,1 V; drenaj ortak bir yük direnci (≈ 5 kΩ) üzerinden. Cihaz V_ds'i ≈ 0,8 V'ta kalıyor, yani nominal
  1,8 V'un altında.
- **Beklenen imza:** bölümlerden ayrı akım ölçülürse filamanın konumu doğrudan görülür. Bölünmede iki ayrı bölüm grubunun
  iletkenleşmesi ve aradakilerin sönmesi beklenir.
- **Kontrol yapısı:** balastsız ve uçtan beslenen bir eşi. Bunda bölünme beklenmez, filamanlar kontaklarda doğar.

## 5. Geçerlilik sınırları

- Exp12'nin bütün varsayımları geçerli: BSIM4 300 K'de çıkarıldı; 4 K düğmeleri varsayım; ısınma, tuzak kinetiği, parazitik bipolar ve
  kuyu kırılımı modelde yok.
- Yerel inhibitör yalnız drenaj hattı. Kaynak hattı direnci ve 2B kanal yapısı dahil değil.
- Tek bir temel parametre seti var (Exp12'nin temsilci seti). Bölünme yalnız λ = 40 µm'de görüldü; λ ile yük direncinin tam haritası
  çıkarılmadı.
- Soy takibi V_DD adımları arasındaki örtüşmeye dayanıyor. Bölünme olayı 4 kat ince adımla doğrulandı.

## Kaynaklar

- Gray–Scott leke çoğalması: [Pearson 1993](https://arxiv.org/abs/patt-sol/9304003)
- S-NDC akım filamanları ve yerel ve global bağlaşım: [Hüpper vd. PRB 47, 15515](https://link.aps.org/doi/10.1103/PhysRevB.47.15515), [Schöll & Drasdo 1990](https://link.springer.com/article/10.1007/BF01309347)
- ESD balast ve filaman (oda sıcaklığı, karşılaştırma): [IISc yayınları](https://faculty.dese.iisc.ac.in/mayank/wp-content/uploads/sites/5/2022/12/55.pdf)
- Kriyojenik CMOS kink: [arXiv 2003.00207](https://arxiv.org/pdf/2003.00207), [arXiv 2104.08467](https://arxiv.org/pdf/2104.08467)
