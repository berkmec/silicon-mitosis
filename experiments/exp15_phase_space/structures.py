"""Exp15 yapi analizi: periyodik bilesen etiketleme, ortusme ile birlesme/bolunme takibi, enine uyum.

Her boyutta (1D/2D/3D) ayni yontem: phi'yi (Gauss, sigma = 1 hucre, periyodik) yumusat, phi > esik
bolgelerini (yuz-komsulugu) etiketle, periyodik sinirlarda birlestir, kucuk parcalari at.
Ardisik karelerde ortusme grafigi: bir ebeveynin >= 2 cocugu -> bolunme; bir cocugun >= 2 ebeveyni -> birlesme.
"""
import numpy as np
from scipy import ndimage


def label_periodic(mask, min_size=1):
    lab, n = ndimage.label(mask)
    if n == 0:
        return lab, 0
    parent = np.arange(n + 1)

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for ax in range(mask.ndim):
        if mask.shape[ax] < 2:
            continue
        a = np.take(lab, 0, axis=ax).ravel()
        b = np.take(lab, -1, axis=ax).ravel()
        m = (a > 0) & (b > 0)
        if m.any():
            for p, q in np.unique(np.stack([a[m], b[m]], 1), axis=0):
                rp, rq = find(p), find(q)
                if rp != rq:
                    parent[max(rp, rq)] = min(rp, rq)
    roots = np.array([find(i) for i in range(n + 1)])
    _, inv = np.unique(roots, return_inverse=True)
    lab = inv[lab]
    nlab = int(inv.max())
    if min_size > 1 and nlab > 0:
        sizes = np.bincount(lab.ravel(), minlength=nlab + 1)
        small = sizes < min_size
        small[0] = False
        if small.any():
            keep = ~small
            newid = np.cumsum(keep) - 1
            newid[small] = 0
            lab = newid[lab]
            nlab = int(newid.max())
    return lab, nlab


def overlap_events(labA, nA, labB, nB):
    m = (labA > 0) & (labB > 0)
    if m.any():
        pairs = np.unique(np.stack([labA[m], labB[m]], 1), axis=0)
    else:
        pairs = np.zeros((0, 2), dtype=int)
    childs = np.bincount(pairs[:, 0], minlength=nA + 1)[1:]
    parents = np.bincount(pairs[:, 1], minlength=nB + 1)[1:]
    return dict(splits=int(np.sum(np.maximum(childs - 1, 0))),
                merges=int(np.sum(np.maximum(parents - 1, 0))),
                births=int(np.sum(parents == 0)), deaths=int(np.sum(childs == 0)))


def smooth(phi, sigma=1.0):
    return ndimage.gaussian_filter(phi, sigma=sigma, mode="wrap")


def _pairs(labA, labB):
    m = (labA > 0) & (labB > 0)
    if not m.any():
        return np.zeros((0, 2), dtype=int)
    return np.unique(np.stack([labA[m], labB[m]], 1), axis=0)


def persistent_splits(pairs, P):
    """Kalici bolunme sayisi: ebeveyn -> >= 2 cocuk; en az iki cocugun torun kumeleri P kare boyunca
    bos olmayan ve ayrik kalmali (yeniden birlesme ya da kaybolma = titresim, sayilmaz). Seri doner."""
    nfr = len(pairs)
    out = np.zeros(nfr, dtype=int)
    for n, pr in enumerate(pairs):
        if pr.size == 0:
            continue
        parents, cnt = np.unique(pr[:, 0], return_counts=True)
        for a in parents[cnt >= 2]:
            kids = pr[pr[:, 0] == a, 1]
            desc = [{int(b)} for b in kids]
            ok = n + P < nfr
            for j in range(1, P + 1):
                if not ok:
                    break
                pj = pairs[n + j]
                new = []
                for s in desc:
                    if not s:
                        new.append(set())
                        continue
                    sel = np.isin(pj[:, 0], list(s))
                    new.append(set(int(c) for c in pj[sel, 1]))
                desc = new
            if not ok:
                continue
            alive = [s for s in desc if s]
            # ayriklik: hicbir torun iki farkli cocuga ait olmasin
            seen, disjoint = set(), True
            for s in alive:
                if seen & s:
                    disjoint = False
                    break
                seen |= s
            if len(alive) >= 2 and disjoint:
                out[n] += len(alive) - 1
    return out


def track(phis, theta, min_size, sigma=1.0, P=5):
    """phis: (nt, *N). Donus: N(t), olay serileri (her adimda bolunme/birlesme/dogum/olum/kalici bolunme)."""
    labs = []
    counts = []
    for ph in phis:
        lab, n = label_periodic(smooth(ph.astype(float), sigma) > theta, min_size)
        labs.append((lab, n))
        counts.append(n)
    ev = dict(splits=[], merges=[], births=[], deaths=[])
    pairs = []
    for (la, na), (lb, nb) in zip(labs[:-1], labs[1:]):
        e = overlap_events(la, na, lb, nb)
        pairs.append(_pairs(la, lb))
        for kk in ev:
            ev[kk].append(e[kk])
    ev = {k: np.array(v) for k, v in ev.items()}
    ev["splits_persistent"] = persistent_splits(pairs, P)
    return np.array(counts), ev


def coherence(phi, h):
    """Alan enerjisinin k_perp = 0 (x'e bagli, dik yonde duzgun) kesri ve dik enerji.
    phi: (*N) gercek; h: hucre boyutlari. W ~ sum k^2 |phi_k|^2."""
    d = phi.ndim
    if d == 1:
        return 1.0, 0.0, 0.0
    pk = np.fft.fftn(phi)
    k2 = 0.0
    for a in range(d):
        k = 2 * np.pi * np.fft.fftfreq(phi.shape[a], d=h[a])
        shape = [1] * d
        shape[a] = k.size
        k2 = k2 + k.reshape(shape) ** 2
    w = k2 * np.abs(pk) ** 2
    par = w[(slice(None),) + (0,) * (d - 1)].sum()
    tot = w.sum()
    return float(par / tot) if tot > 0 else 1.0, float(par), float(tot - par)


def hole_centers_1d(prof, h, theta):
    """1D profil (periyodik) icin esigi asan yerel maksimumlarin konumlari (parabolik alt-hucre)."""
    n = prof.size
    pl, pr = np.roll(prof, 1), np.roll(prof, -1)
    idx = np.where((prof > pl) & (prof >= pr) & (prof > theta))[0]
    out = []
    for i in idx:
        a, b, c = pl[i], prof[i], pr[i]
        den = a - 2 * b + c
        off = 0.5 * (a - c) / den if den != 0 else 0.0
        out.append(((i + off) % n) * h)
    return np.array(out)


def kink_amplitude(phi2d, h, X, half_window):
    """2D (Nx, Ny): her duzlemsel delik merkezi X_j cevresinde her y-satiri icin x_c(y) = argmax phi;
    donus: deliklerin x_c(y) standart sapmalarinin ortalamasi (enine kayma/kink genligi)."""
    Nx = phi2d.shape[0]
    xs = np.arange(Nx) * h
    amps = []
    for Xj in X:
        dxp = (xs - Xj + 0.5 * Nx * h) % (Nx * h) - 0.5 * Nx * h
        sel = np.where(np.abs(dxp) <= half_window)[0]
        sub = phi2d[sel, :]
        im = np.argmax(sub, axis=0)
        # parabolik alt-hucre
        xc = []
        for j, i in enumerate(im):
            if 0 < i < len(sel) - 1:
                a, b, c = sub[i - 1, j], sub[i, j], sub[i + 1, j]
                den = a - 2 * b + c
                off = 0.5 * (a - c) / den if den != 0 else 0.0
            else:
                off = 0.0
            xc.append(dxp[sel[i]] + off * h)
        amps.append(np.std(xc))
    return float(np.mean(amps)) if amps else float("nan")
