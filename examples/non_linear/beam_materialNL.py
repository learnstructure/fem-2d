from fem2d import Structure, Node
from fem2d.materials import BilinearMaterial
from fem2d.sections import FiberSection
from fem2d import BeamElementMaterialNL

#units in kN, m

# Create structure
node1 = Node(1, 0.0, 0.0)
node2 = Node(2, 2.0, 0.0)

# Support Node 1 (fixed cantilever)
node1.set_support(ux_fixed=True, uy_fixed=True, rz_fixed=True)

# Material & Section
E = 200e6
sigma_y = 250e3
E_t = 2e6 
material = BilinearMaterial(E, sigma_y, E_t)
section = FiberSection.from_rectangular_fibers(width=0.1, depth=0.2, material=material, n_fibers=20)

# Create element
elem = BeamElementMaterialNL(1, node1, node2, section, n_int_pts=3)

structure = Structure()
structure.add_node(node1)
structure.add_node(node2)
structure.add_element(elem)

# Linear load: -10 kN (stress is below yield)
node2.set_load(fx=0.0, fy=-10, mz=0.0)
structure.solve_nonlinear(n_steps=1)
disp = structure.disp[node2.dofs]
print(disp)

#Large load
node2.set_load(fx=0.0, fy=-500, mz=0.0)
# Solve with multiple steps
structure.solve_nonlinear(n_steps=50)
disp_plastic = structure.disp[node2.dofs]
print(disp_plastic)