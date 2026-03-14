import numpy as np
from fem2d import Structure, Node, ElasticMaterial, Section, Results, BeamElement

# 1. Create nodes
node1 = Node(1, 0, 0)
node2 = Node(2, 0, 400)
node3 = Node(3, 200, 400)

# 2. Define material and section properties
E = 12000
A = 1000
I = 25
material = ElasticMaterial(E)
section = Section(A, I)

# 3. Create elements
elem1 = BeamElement(1, node1, node2, material, section.A, section.I)
elem2 = BeamElement(2, node2, node3, material, section.A, section.I)

# 4. Create structure and add nodes/elements
structure = Structure()
structure.add_node(node1)
structure.add_node(node2)
structure.add_node(node3)
structure.add_element(elem1)
structure.add_element(elem2)

# 5. Apply supports (all DOFs fixed)
node1.set_support(ux_fixed=True, uy_fixed=True, rz_fixed=True)
node3.set_support(ux_fixed=True, uy_fixed=True, rz_fixed=True)

# 6. Apply nodal load at node 2 (moment -300 kip-in)
node2.set_load(fx=0, fy=0, mz=-300)

# 7. Solve
structure.solve()

results = Results(structure)
print("Node Displacements:\n", results.node_displacements())
print("Reactions:\n", results.reactions())
print("Element Forces:\n", results.element_forces())
