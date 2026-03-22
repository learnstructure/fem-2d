# Example 5.1: Logan D. L. (2017), A First Course in the Finite Element Method
# units in kips, inches

from fem2d import (
    Structure,
    Node,
    ElasticMaterial,
    Section,
    BeamElement,
    Results,
    DrawStructure,
)

import matplotlib.pyplot as plt

# 1. Create nodes
node1 = Node(1, 0, 0)
node2 = Node(2, 0, 10 * 12)  # 10 ft = 120 in
node3 = Node(3, 10 * 12, 10 * 12)
node4 = Node(4, 10 * 12, 0)

# 2. Material and section properties
E = 30e3
A = 10
I_full = 200
I_half = I_full / 2

material = ElasticMaterial(E)  # same material for all elements
section_full = Section(A, I_full)
section_half = Section(A, I_half)

# 3. Create elements
elem1 = BeamElement(1, node1, node2, material, section_full.A, section_full.I)
elem2 = BeamElement(2, node2, node3, material, section_half.A, section_half.I)
elem3 = BeamElement(3, node3, node4, material, section_full.A, section_full.I)

# 4. Create structure and add nodes/elements
structure = Structure()
structure.add_node(node1)
structure.add_node(node2)
structure.add_node(node3)
structure.add_node(node4)
structure.add_element(elem1)
structure.add_element(elem2)
structure.add_element(elem3)

# 5. Apply supports (all DOFs fixed)
node1.set_support(ux_fixed=True, uy_fixed=True, rz_fixed=True)
node4.set_support(ux_fixed=True, uy_fixed=True, rz_fixed=True)

# 6. Apply nodal loads
node2.set_load(fx=10, fy=0, mz=0)  # horizontal point load at node 2
node3.set_load(fx=0, fy=0, mz=5)  # moment at node 3

# 7. Solve
structure.solve()

# 8. Print results
results = Results(structure)
disp = results.node_displacements()
reactions = results.reactions()
el_forces = results.element_forces()
print(disp["theta"][2])  # result = -0.0014859999862147104
print(el_forces["m_i"][2])  # result = 226.19833955969983

drawer = DrawStructure(structure, scale=100)
drawer.draw()
