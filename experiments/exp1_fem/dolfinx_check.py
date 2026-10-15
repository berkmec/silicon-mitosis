"""Independent cross-check (dolfinx, runs inside WSL): same mesh as scikit-fem, concentric spheres."""
import sys

import basix.ufl
import numpy as np
import ufl
from dolfinx import fem, mesh
from dolfinx.fem import petsc
from mpi4py import MPI

d = np.load(sys.argv[1])
p, t, c_skfem = d["p"], d["t"].astype(np.int64), float(d["C_fem"])
domain = mesh.create_mesh(MPI.COMM_WORLD, t, p, ufl.Mesh(basix.ufl.element("Lagrange", "tetrahedron", 1, shape=(3,))))
V = fem.functionspace(domain, ("Lagrange", 1))
xd = V.tabulate_dof_coordinates()
r = np.linalg.norm(xd, axis=1)
dofs_in = np.where(r < 0.5 + 1e-7)[0].astype(np.int32)
dofs_out = np.where(r > 1.0 - 1e-7)[0].astype(np.int32)
bcs = [fem.dirichletbc(fem.Constant(domain, 1.0), dofs_in, V),
       fem.dirichletbc(fem.Constant(domain, 0.0), dofs_out, V)]
u, v = ufl.TrialFunction(V), ufl.TestFunction(V)
a = ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx
L = fem.Constant(domain, 0.0) * v * ufl.dx
uh = petsc.LinearProblem(a, L, bcs=bcs, petsc_options={"ksp_type": "preonly", "pc_type": "lu"}).solve()
C = domain.comm.allreduce(fem.assemble_scalar(fem.form(ufl.inner(ufl.grad(uh), ufl.grad(uh)) * ufl.dx)), op=MPI.SUM)
exact = 4 * np.pi
print(f"dolfinx {__import__('dolfinx').__version__}: nodes={p.shape[0]} dofs_in={len(dofs_in)} dofs_out={len(dofs_out)}")
print(f"C_dolfinx = {C:.10f}")
print(f"C_skfem   = {c_skfem:.10f}")
print(f"|dolfinx-skfem|/skfem = {abs(C - c_skfem) / c_skfem:.3e}")
print(f"C_analytic = {exact:.10f}  (dolfinx hatası {abs(C - exact) / exact:.3e})")
