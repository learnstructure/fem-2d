"""
Example: Cantilever Beam with Shear Deformation (Timoshenko Beam)
Reference: Kassimali, A. (2022). Matrix Analysis of Structures, 3rd ed.,
           Chapter 9: Special Topics, Section 9.7: Shear Deformations,
           Example 9.6, pp. 551-552.

Units: kN, m

Problem Description:
A cantilever beam of rectangular cross section (depth d = 0.5 m, width b = 0.2 m)
and span L = 2.5 m (depth-to-span ratio d/L = 0.2) is subjected to a vertical
concentrated load P = 100 kN downward at its free tip.

Material properties:
  E = 200 GPa = 200e6 kN/m^2
  G = 77.2 GPa = 77.2e6 kN/m^2  (Poisson's ratio nu = 0.2953)

Section properties:
  A = b * d = 0.10 m^2
  I = b * d^3 / 12 = 0.00208333 m^4
  For a rectangular section, shear shape factor fs = 1.2 => As = A / 1.2 = 5/6 * A

Analytical solutions from Kassimali (Eq. 9.60 & Example 9.6):
  beta_s = 12 * E * I / (G * As * L^2) = 0.124352
  Bending tip deflection:  Delta_M = P * L^3 / (3 * E * I) = 1.25000 mm
  Shear tip deflection:    Delta_s = (beta_s / 4) * Delta_M = 0.03886 mm (+3.11%)
  Total tip deflection:    Delta = Delta_M + Delta_s = 1.28886 mm
  Tip rotation:            theta = -P * L^2 / (2 * E * I) = -0.000750 rad
"""

from fem2d import (
    Structure,
    Node,
    ElasticMaterial,
    Section,
    BeamElement,
    SimpleFrame,
    Results,
    DrawStructure,
)
from fem2d.loads import PointLoad

# ==============================================================================
# 1. Parameters
# ==============================================================================
L = 2.5       # m
b = 0.2       # m
d = 0.5       # m
P = 100.0     # kN
E = 200.0e6   # kN/m^2
G = 77.2e6    # kN/m^2

# Material and section
material = ElasticMaterial(E=E, G=G)
section = Section.from_rectangle(width=b, depth=d)

beta_s = (12.0 * E * section.I) / (G * section.As * L**2)
delta_M_theory = (P * L**3) / (3.0 * E * section.I)
delta_s_theory = (beta_s / 4.0) * delta_M_theory
delta_total_theory = delta_M_theory + delta_s_theory
theta_theory = -(P * L**2) / (2.0 * E * section.I)

print("=" * 70)
print("CANTILEVER BEAM WITH SHEAR DEFORMATION (KASSIMALI EXAMPLE 9.6)")
print("=" * 70)
print(f"Span L = {L:.2f} m, Section = {b:.2f} m x {d:.2f} m (d/L = {d/L:.2f})")
print(f"beta_s parameter = {beta_s:.6f}")
print(f"Analytical Bending Deflection (Delta_M) : {delta_M_theory * 1e3:.5f} mm")
print(f"Analytical Shear Deflection   (Delta_s) : {delta_s_theory * 1e3:.5f} mm (+{(delta_s_theory / delta_M_theory) * 100:.2f}%)")
print(f"Analytical Total Deflection   (Delta)   : {delta_total_theory * 1e3:.5f} mm")
print(f"Analytical Tip Rotation       (theta)   : {theta_theory:.6f} rad")
print("-" * 70)

# ==============================================================================
# 2. Case A: Classical Euler-Bernoulli (Shear deformation neglected)
# ==============================================================================
n1_eb = Node(1, 0.0, 0.0)
n2_eb = Node(2, L, 0.0)
n1_eb.set_support(ux_fixed=True, uy_fixed=True, rz_fixed=True)

elem_eb = BeamElement(1, n1_eb, n2_eb, material, section, include_shear_deformation=False)

struct_eb = Structure()
struct_eb.add_node(n1_eb)
struct_eb.add_node(n2_eb)
struct_eb.add_element(elem_eb)
struct_eb.add_load(PointLoad(n2_eb, fx=0.0, fy=-P))
struct_eb.solve()

res_eb = Results(struct_eb)
disp_eb = res_eb.node_displacements()
delta_eb_fem = abs(disp_eb["uy"][1])
theta_eb_fem = disp_eb["theta"][1]

print(f"Case A: Euler-Bernoulli FEA Tip Deflection : {delta_eb_fem * 1e3:.5f} mm")
print(f"        Euler-Bernoulli FEA Tip Rotation   : {theta_eb_fem:.6f} rad")

# ==============================================================================
# 3. Case B: Timoshenko Beam (Shear deformation included via BeamElement)
# ==============================================================================
n1_timo = Node(1, 0.0, 0.0)
n2_timo = Node(2, L, 0.0)
n1_timo.set_support(ux_fixed=True, uy_fixed=True, rz_fixed=True)

elem_timo = BeamElement(1, n1_timo, n2_timo, material, section, include_shear_deformation=True)

struct_timo = Structure()
struct_timo.add_node(n1_timo)
struct_timo.add_node(n2_timo)
struct_timo.add_element(elem_timo)
struct_timo.add_load(PointLoad(n2_timo, fx=0.0, fy=-P))
struct_timo.solve()

res_timo = Results(struct_timo)
disp_timo = res_timo.node_displacements()
delta_timo_fem = abs(disp_timo["uy"][1])
theta_timo_fem = disp_timo["theta"][1]

print(f"Case B: Timoshenko FEA Tip Deflection      : {delta_timo_fem * 1e3:.5f} mm")
print(f"        Timoshenko FEA Tip Rotation        : {theta_timo_fem:.6f} rad")
print(f"        Error in Deflection vs Theory      : {abs(delta_timo_fem - delta_total_theory):.2e} m")
print(f"        Error in Rotation vs Theory        : {abs(theta_timo_fem - theta_theory):.2e} rad")

# ==============================================================================
# 4. Case C: Timoshenko Beam using SimpleFrame API
# ==============================================================================
frame = SimpleFrame()
frame.add_node(1, 0.0, 0.0)
frame.add_node(2, L, 0.0)
frame.add_support(1, [1, 1, 1])
frame.add_frame(
    1,
    1,
    2,
    E=E,
    A=section.A,
    I=section.I,
    As=section.As,
    G=G,
    include_shear_deformation=True,
)
frame.add_node_load(2, [0.0, -P, 0.0])
frame.solve()

res_frame = Results(frame)
disp_frame = res_frame.node_displacements()
delta_frame_fem = abs(disp_frame["uy"][1])

print(f"Case C: SimpleFrame API Tip Deflection     : {delta_frame_fem * 1e3:.5f} mm")
print("=" * 70)

# Optional visualization
# drawer = DrawStructure(struct_timo, scale=1000)
# drawer.draw()
