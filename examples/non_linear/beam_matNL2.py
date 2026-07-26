from fem2d import Structure, Node
from fem2d.materials import BilinearMaterial
from fem2d.sections import FiberSection
from fem2d import BeamElementMaterialNL

#units in kN, m

# Create structure
node1 = Node(1, 0.0, 0.0)
node2 = Node(2, 5.0, 0.0)

# Support Node 1 (fixed cantilever)
node1.set_support(1, 1, 1)

# Material & Section
E = 200e6
sigma_y = 400000000e3
E_t = 0.05 * E
material = BilinearMaterial(E, sigma_y, E_t, beta=0.0)
section = FiberSection.from_rectangular_fibers(width=0.05, depth=0.05, material=material, n_fibers=20)
# Create element
elem = BeamElementMaterialNL(1, node1, node2, section, n_int_pts=3)

structure = Structure()
structure.add_node(node1)
structure.add_node(node2)
structure.add_element(elem)

node2.set_load(fx=0.0, fy=-2.5, mz=0.0)
structure.solve_nonlinear(n_steps=1, tolerance=1e-8, max_iter=100)
disp = structure.disp[node2.dofs]
print(disp)      
