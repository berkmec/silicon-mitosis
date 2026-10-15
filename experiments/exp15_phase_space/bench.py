"""Hiz olcumu (kisa): 1D, 2D, 2D+B, 3D adim suresi."""
import os
os.environ["OMP_NUM_THREADS"] = "1"
import sys
if sys.platform == "win32":  # below-normal priority, Windows only
    import ctypes
    from ctypes import wintypes

    k32 = ctypes.windll.kernel32
    k32.GetCurrentProcess.restype = wintypes.HANDLE
    k32.SetPriorityClass.argtypes = [wintypes.HANDLE, wintypes.DWORD]
    k32.SetPriorityClass(k32.GetCurrentProcess(), 0x4000)
import sys
import time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pic

rng = np.random.default_rng(0)
for (L, N, npart, nv, B) in [((102.6,), (512,), 204800, 1, 0.0),
                             ((41.05, 41.05), (128, 128), 1048576, 2, 0.0),
                             ((41.05, 41.05), (128, 128), 1048576, 3, 2.0),
                             ((20.52, 20.52, 20.52), (32, 32, 32), 1048576, 3, 0.0)]:
    x, v = pic.beams(L, npart, 1.0, 0.2, nv, rng)
    sim = pic.PIC(L, N, x, v, 0.1, B=B)
    sim.init_leapfrog()
    t0 = time.perf_counter()
    nst = 4
    for _ in range(nst):
        sim.step(diag=True)
    dt = (time.perf_counter() - t0) / nst
    print(f"d={len(N)} N={N} Np={npart} B={B}: {dt*1e3:.0f} ms/step", flush=True)
    del sim, x, v
