"""Exp15 - d-boyutlu elektrostatik PIC cekirdegi (1D1V, 2D2V/2D3V, 3D3V), ayni kod yolu.

Normalize birimler: omega_p = 1 (toplam elektron yogunlugu n0 = 1, e = m_e = eps0 = 1).
Elektronlar q/m = -1; iyonlar hareketsiz, duzgun +1 notrlestirici arka plan.
  Poisson  : laplacian(phi) = n_e - 1   (rho = 1 - n_e),  E = -grad(phi)   [FFT, periyodik kutu]
  Yukleme  : bulut-hucre (CIC, dogrusal); interpolasyon ayni agirliklarla (momentum korunumlu sema)
  Zaman    : leapfrog (B = 0) veya Boris (duzgun B, +x yonunde; 3 hiz bileseni gerekir)
Parcacik dizileri bilesen-once tutulur: x (d, Np), v (nv, Np), nv >= d.
"""
import numpy as np
from scipy.special import wofz


class PIC:
    def __init__(self, L, N, x, v, dt, B=0.0, smooth=0, solver="spectral"):
        self.L = np.asarray(L, dtype=float).reshape(-1)
        self.N = tuple(int(n) for n in N)
        self.d = len(self.N)
        if self.L.size != self.d:
            raise ValueError("L ve N boyutlari uyusmuyor")
        self.h = self.L / np.asarray(self.N, dtype=float)
        self.ncell = int(np.prod(self.N))
        self.x = np.array(x, dtype=float, copy=True, order="C")   # cagirana ait diziler degismesin
        self.v = np.array(v, dtype=float, copy=True, order="C")
        if self.x.shape[0] != self.d or self.v.shape[0] < self.d:
            raise ValueError("x (d,Np), v (nv>=d,Np) olmali")
        self.nv = self.v.shape[0]
        self.np = self.x.shape[1]
        self.dt = float(dt)
        self.B = float(B)
        if self.B != 0.0 and self.nv < 3:
            raise ValueError("B != 0 icin 3 hiz bileseni gerekir")
        self.V = float(np.prod(self.L))
        self.w = self.V / self.np          # makro-parcacik agirligi (n0 = 1)
        self.dV = self.V / self.ncell
        self.t = 0.0
        self.solver = solver
        self.smooth = int(smooth)
        self._wrap(full=True)
        # Fourier izgarasi (rfftn: son eksen yarim)
        kk = []
        for a in range(self.d):
            if a < self.d - 1:
                k = 2 * np.pi * np.fft.fftfreq(self.N[a], d=self.h[a])
            else:
                k = 2 * np.pi * np.fft.rfftfreq(self.N[a], d=self.h[a])
            shape = [1] * self.d
            shape[a] = k.size
            kk.append(k.reshape(shape))
        self.k = kk
        if solver == "spectral":
            K2 = sum(k ** 2 for k in kk)
            grad = []
            for a, k in enumerate(kk):
                g = k.copy()
                if self.N[a] % 2 == 0:  # Nyquist bileseninin turevi sifir (reel alan)
                    g[np.isclose(np.abs(g), np.pi / self.h[a])] = 0.0
                grad.append(g)
        elif solver == "fd":
            K2 = sum((2 * np.sin(k * h / 2) / h) ** 2 for k, h in zip(kk, self.h))
            grad = [np.sin(k * h) / h for k, h in zip(kk, self.h)]
        else:
            raise ValueError(solver)
        K2 = np.asarray(K2, dtype=float)
        inv = np.zeros_like(K2)
        nz = K2 > 0
        inv[nz] = 1.0 / K2[nz]
        filt = np.ones_like(K2)
        if self.smooth:
            for k, h in zip(kk, self.h):
                filt = filt * np.cos(k * h / 2) ** (2 * self.smooth)
        self.G = inv * filt          # phi_k = G * rho_k
        self.grad = grad
        self._idx = None
        self._wts = None

    # ------------------------------------------------------------------ yardimcilar
    def _wrap(self, full=False):
        """Periyodik sarma. Adim basina |v dt| < L varsayimiyla maske ile (np.remainder'dan hizli)."""
        for a in range(self.d):
            xa = self.x[a]
            if full:
                np.remainder(xa, self.L[a], out=xa)
                continue
            m = xa >= self.L[a]
            xa[m] -= self.L[a]
            m = xa < 0.0
            xa[m] += self.L[a]

    def _prep(self):
        """Mevcut x icin 2^d kose duz-indeksi ve CIC agirliklari (C sirasi)."""
        idx = wts = None
        for a in range(self.d):
            s = self.x[a] * (1.0 / self.h[a])
            j = s.astype(np.int64)            # s >= 0 -> kesme = taban
            f1 = s
            f1 -= j                           # s yerinde f1'e donusur
            n = self.N[a]
            j[j >= n] -= n                    # x == L kenar durumu
            jp = j + 1
            jp[jp == n] = 0
            f0 = 1.0 - f1
            if idx is None:
                idx, wts = [j, jp], [f0, f1]
            else:
                ni, nw = [], []
                for I, W in zip(idx, wts):
                    base = I * n
                    ni.append(base + j)
                    ni.append(base + jp)
                    nw.append(W * f0)
                    nw.append(W * f1)
                idx, wts = ni, nw
        self._idx, self._wts = idx, wts

    def deposit(self):
        cnt = np.zeros(self.ncell)
        for I, W in zip(self._idx, self._wts):
            cnt += np.bincount(I, weights=W, minlength=self.ncell)
        return (cnt * (self.w / self.dV)).reshape(self.N)

    def solve(self, ne):
        rk = np.fft.rfftn(1.0 - ne)
        phik = rk * self.G
        ax = tuple(range(self.d))
        phi = np.fft.irfftn(phik, s=self.N, axes=ax)
        E = [np.fft.irfftn(-1j * g * phik, s=self.N, axes=ax) for g in self.grad]
        return phi, E

    def gather(self, E):
        if self.d == 1:
            e = E[0].ravel()
            Ep = np.zeros(self.np)
            for I, W in zip(self._idx, self._wts):
                Ep += W * e[I]
            return Ep[None, :]
        fl = [e.ravel() for e in E]
        Ep = np.zeros((self.d, self.np))
        for I, W in zip(self._idx, self._wts):
            for a in range(self.d):
                Ep[a] += W * fl[a][I]
        return Ep

    def push(self, Ep):
        dt = self.dt
        if self.B == 0.0:
            for a in range(self.d):
                self.v[a] -= dt * Ep[a]          # q/m = -1
        else:
            hdt = 0.5 * dt
            for a in range(self.d):
                self.v[a] -= hdt * Ep[a]
            t = -self.B * hdt                    # (q/m) B dt/2, B || +x
            s = 2.0 * t / (1.0 + t * t)
            vy, vz = self.v[1], self.v[2]
            vyp = vy + vz * t
            vzp = vz - vy * t
            vy += vzp * s
            vz -= vyp * s
            for a in range(self.d):
                self.v[a] -= hdt * Ep[a]
        for a in range(self.d):
            self.x[a] += dt * self.v[a]
        self._wrap()

    # ------------------------------------------------------------------ zaman adimi
    def fields(self):
        self._prep()
        ne = self.deposit()
        phi, E = self.solve(ne)
        return ne, phi, E

    def init_leapfrog(self):
        """v(0) -> v(-dt/2) (yalniz E ile yarim adim geri)."""
        ne, phi, E = self.fields()
        Ep = self.gather(E)
        for a in range(self.d):
            self.v[a] += 0.5 * self.dt * Ep[a]

    def step(self, diag=False):
        """Bir leapfrog adimi. diag=True ise t^n zamanindaki alanlar ve enerjiler dondurulur."""
        ne, phi, E = self.fields()
        Ep = self.gather(E)
        out = None
        if diag:
            K0 = 0.5 * self.w * np.einsum("ij,ij->i", self.v, self.v)
        self.push(Ep)
        if diag:
            K1 = 0.5 * self.w * np.einsum("ij,ij->i", self.v, self.v)
            WE = 0.5 * self.dV * sum(float(np.sum(e * e)) for e in E)
            out = dict(t=self.t, WE=WE, KEc=0.5 * (K0 + K1), ne=ne, phi=phi, E=E)
        self.t += self.dt
        return out


# ---------------------------------------------------------------------- yukleme
def bitrev(n_bits, n):
    """0..n-1 icin taban-2 van der Corput dizisi (bit ters cevirme), (0,1) icinde."""
    i = np.arange(n, dtype=np.uint64)
    r = np.zeros(n, dtype=np.uint64)
    for b in range(n_bits):
        r |= ((i >> np.uint64(b)) & np.uint64(1)) << np.uint64(n_bits - 1 - b)
    return (r.astype(float) + 0.5) / float(1 << n_bits)


def beams(L, npart, v0, vth, nv, rng):
    """Iki simetrik Maxwell demeti (+-v0, x yonunde), her biri yogunluk 1/2; rastgele konumlar."""
    L = np.asarray(L, float)
    d = L.size
    x = rng.random((d, npart)) * L[:, None]
    v = rng.normal(0.0, vth, (nv, npart))
    half = npart // 2
    v[0, :half] += v0
    v[0, half:] -= v0
    return x, v


def displace_x(x, L0, modes, amps, phases=None):
    """x-yonunde duzlemsel yerdegistirme: x += sum (a_m/k_m) sin(k_m x + p_m) -> dn/n ~ -a_m cos(...)."""
    out = x[0].copy()
    if phases is None:
        phases = np.zeros(len(modes))
    for m, a, p in zip(modes, amps, phases):
        k = 2 * np.pi * m / L0
        out += (a / k) * np.sin(k * x[0] + p)
    x[0] = np.remainder(out, L0)
    return x


# ---------------------------------------------------------------------- dispersiyon
def Zfun(z):
    """Plazma dispersiyon fonksiyonu Z(z) = i sqrt(pi) w(z) (tum karmasik duzlemde analitik devam)."""
    return 1j * np.sqrt(np.pi) * wofz(z)


def eps_beams(omega, k, comps):
    """Maxwell bilesenleri [(n_b, v_b, vth_b), ...] icin 1B elektrostatik dielektrik fonksiyonu."""
    e = 1.0 + 0j
    for nb, vb, vt in comps:
        z = (omega - k * vb) / (np.sqrt(2.0) * k * vt)
        e += nb / (k * vt) ** 2 * (1.0 + z * Zfun(z))
    return e


def newton_root(f, w0, tol=1e-12, maxit=100):
    w = complex(w0)
    for _ in range(maxit):
        fw = f(w)
        hstep = 1e-7 * max(1.0, abs(w))
        df = (f(w + hstep) - f(w - hstep)) / (2 * hstep)
        dw = fw / df
        w -= dw
        if abs(dw) < tol:
            return w
    raise RuntimeError("Newton yakinsamadi")


def landau_root(k, w0=1.4 - 0.15j):
    return newton_root(lambda w: eps_beams(w, k, [(1.0, 0.0, 1.0)]), w0)


def cold_two_stream_growth(k, v0=1.0):
    """1 = (1/2)[1/(w-kv0)^2 + 1/(w+kv0)^2] -> w^4 - (2K^2+1) w^2 + (K^4-K^2) = 0 (K = k v0); max Im(w)."""
    K = k * v0
    r = np.roots([1.0, 0.0, -(2 * K * K + 1.0), 0.0, K ** 4 - K * K])
    return float(np.max(r.imag))


def warm_two_stream_growth(k, v0, vth):
    """Simetrik Maxwell demetleri: saf buyuyen kok w = i*gamma (varsa)."""
    g0 = cold_two_stream_growth(k, v0)
    if g0 <= 0:
        return 0.0
    comps = [(0.5, v0, vth), (0.5, -v0, vth)]
    try:
        w = newton_root(lambda w: eps_beams(w, k, comps), 1j * g0)
    except RuntimeError:
        return float("nan")
    return float(w.imag) if w.imag > 0 else 0.0


def warm_two_stream_curve(ks, v0, vth):
    """k uzerinde devamlilik (onceki kok baslangic tahmini) ile simetrik sicak demetlerin buyume orani."""
    comps = [(0.5, v0, vth), (0.5, -v0, vth)]
    out = []
    w = 1j * max(cold_two_stream_growth(ks[0], v0), 1e-3)
    for k in ks:
        try:
            w = newton_root(lambda om: eps_beams(om, k, comps), w)
        except (RuntimeError, ZeroDivisionError, FloatingPointError):
            out.append(np.nan)
            continue
        out.append(w.imag if w.imag > 0 else 0.0)
        if w.imag <= 1e-6:
            w = 1j * 1e-3
    return np.array(out)


def cic_shape(k, h):
    return np.sinc(k * h / (2 * np.pi)) ** 2       # np.sinc(x) = sin(pi x)/(pi x)


def langmuir_freq_pred(k, h, dt, solver="spectral", pmax=5):
    """Soguk plazma, CIC + leapfrog: w_g^2 = (1/k) sum_p q_p S(q_p)^2 (spektral), sonra leapfrog."""
    kg = 2 * np.pi / h
    s = 0.0
    for p in range(-pmax, pmax + 1):
        q = k + p * kg
        s += q * cic_shape(q, h) ** 2
    if solver == "spectral":
        wg2 = s / k
    else:
        K2 = (2 * np.sin(k * h / 2) / h) ** 2
        wg2 = s * (np.sin(k * h) / h) / K2
    wg = np.sqrt(wg2)
    return float(wg), float(2.0 / dt * np.arcsin(0.5 * dt * wg))
