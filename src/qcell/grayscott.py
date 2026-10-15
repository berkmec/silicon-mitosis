"""Gray-Scott reaction-diffusion model exactly as in Pearson, Science 261, 189 (1993) (arXiv patt-sol/9304003):
    U_t = Du lap U - U V^2 + F (1 - U),   V_t = Dv lap V + U V^2 - (F + k) V
Du = 2e-5, Dv = 1e-5, periodic box, 5-point Laplacian, forward Euler, dt = 1, mesh 2.5/256."""
import numpy as np
from scipy import ndimage


def uniform_states(F, k):
    """Homogeneous steady states: (1, 0) and, if F >= 4 (F + k)^2, v = [F +- sqrt(F^2 - 4F(F+k)^2)] / (2(F+k)), u = (F+k)/v."""
    out = [(1.0, 0.0)]
    disc = F**2 - 4 * F * (F + k) ** 2
    if disc >= 0:
        for sgn in (1.0, -1.0):
            v = (F + sgn * np.sqrt(disc)) / (2 * (F + k))
            out.append(((F + k) / v, v))
    return out


def rhs_kinetics(u, v, F, k):
    uvv = u * v * v
    return -uvv + F * (1 - u), uvv - (F + k) * v


class GrayScott:
    def __init__(self, n=128, h=2.5 / 256, Du=2e-5, Dv=1e-5, F=0.02, k=0.059):
        self.n, self.h, self.Du, self.Dv, self.F, self.k = n, h, Du, Dv, F, k

    def seed(self, rng, half=10, noise=0.01):
        n = self.n
        U, V = np.ones((n, n)), np.zeros((n, n))
        s = slice(n // 2 - half, n // 2 + half)
        U[s, s], V[s, s] = 0.5, 0.25
        U *= 1 + noise * (2 * rng.random((n, n)) - 1)
        V *= 1 + noise * (2 * rng.random((n, n)) - 1)
        return U, V

    def lap(self, Z):
        return (np.roll(Z, 1, 0) + np.roll(Z, -1, 0) + np.roll(Z, 1, 1) + np.roll(Z, -1, 1) - 4 * Z) / self.h**2

    def run(self, U, V, steps, dt=1.0, every=None, cb=None):
        for t in range(1, steps + 1):
            du, dv = rhs_kinetics(U, V, self.F, self.k)
            U, V = U + dt * (self.Du * self.lap(U) + du), V + dt * (self.Dv * self.lap(V) + dv)
            if cb is not None and t % every == 0 and cb(t, U, V) is False:
                break
        return U, V

    def spots(self, V, thr=0.1):
        """Count V > thr components on the periodic grid; perc = some component wraps the box (stripes/labyrinth)."""
        n = self.n
        big = np.tile(V > thr, (3, 3))
        lab, k = ndimage.label(big)
        if k == 0:
            return dict(n=0, perc=False, area_frac=0.0)
        com = ndimage.center_of_mass(big, lab, range(1, k + 1))
        slices = ndimage.find_objects(lab)
        cnt, perc = 0, False
        for i in range(k):
            if not (n <= com[i][0] < 2 * n and n <= com[i][1] < 2 * n):
                continue
            if any(sl.stop - sl.start >= 2 * n for sl in slices[i]):
                perc = True
            else:
                cnt += 1
        return dict(n=cnt, perc=perc, area_frac=float((V > thr).mean()))
