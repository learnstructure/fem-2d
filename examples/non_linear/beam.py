import numpy as np
from fem2d import Structure, Node, ElasticMaterial, NewtonRaphsonSolver, BeamElementNL

# -------------------------------------------------------------------
# Example 7.1.1 – Cantilever beam with end vertical load
# Data from paper: L = 10 m, b = 0.5 m, h = 0.25 m, E = 200 GPa
# -------------------------------------------------------------------

# Geometry
L = 10.0
b = 0.5
h = 0.25
A = b * h
I = b * h**3 / 12
E = 200e9  # Pa

material = ElasticMaterial(E)

# Create 5 elements (6 nodes)
n_elem = 5
x = np.linspace(0, L, n_elem + 1)
nodes = []
for i, xi in enumerate(x):
    nodes.append(Node(i + 1, xi, 0.0))

structure = Structure()
for node in nodes:
    structure.add_node(node)

# Create beam elements
for i in range(n_elem):
    elem = BeamElementNL(i + 1, nodes[i], nodes[i + 1], material, A, I)
    structure.add_element(elem)

# Support: node 1 fully fixed
nodes[0].set_support(ux_fixed=True, uy_fixed=True, rz_fixed=True)

# Load: vertical down at tip, corresponding to PL²/EI = 10
EI = E * I
P_target = 10.0 * EI / L**2
nodes[-1].set_load(fx=0.0, fy=-P_target, mz=0.0)

# Solve (full Newton‑Raphson without load stepping)
structure.solve_nonlinear(tolerance=1e-8, max_iter=30)

# Results
print("Cantilever beam with end vertical load")
print(f"PL²/EI = 10")
print("Tip displacements:")
disp = structure.disp[nodes[-1].dofs]
print(f"  ux = {disp[0]:.6f} m   (U/L = {disp[0]/L:.6f})")
print(f"  uy = {disp[1]:.6f} m   (W/L = {disp[1]/L:.6f})")
print(f"  theta = {disp[2]:.6f} rad")

# Optional: compare with paper values (Table 1)
paper_U_L = 0.55131
paper_W_L = 0.81471
print(f"\nPaper (Present) U/L = {paper_U_L:.6f}, W/L = {paper_W_L:.6f}")
