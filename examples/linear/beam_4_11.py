# Example 4.11: Logan D. L. (2022), A First Course in the Finite Element Method
# units in kN, m

from fem2d.loads import PointLoad, DistributedLoad, ElementPointLoad
from fem2d import (
    Structure,
    Node,
    ElasticMaterial,
    Section,
    BeamElement, BeamWithHingesElement,
    Results,
    DrawStructure,
)

import matplotlib.pyplot as plt

# 1. Create nodes
node1 = Node(1, 0, 0)
node2 = Node(2, 2, 0)
node3 = Node(3, 3, 0)
node4 = Node(4, 4, 0)

# 2. Material and section properties
E = 210e6  # kN/m^2
A = 6500e-6  # m^2
I = 2e-4  # m^4
material = ElasticMaterial(E)  # same material for all elements
section = Section(A, I)

# 3. Create elements
elem1 = BeamElement(1, node1, node2, material, section.A, section.I)
elem2 = BeamWithHingesElement(2, node2, node3, material, section.A, section.I, hinge_i=False, hinge_j=True)
elem3 = BeamElement(3, node3, node4, material, section.A, section.I)

# elem2 = BeamElement(2, node2, node3, material, section.A, section.I)
# elem3 = BeamWithHingesElement(3, node3, node4, material, section.A, section.I, hinge_i=True, hinge_j=False)

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
node2.set_support(ux_fixed=True, uy_fixed=True, rz_fixed=False)
node4.set_support(ux_fixed=True, uy_fixed=True, rz_fixed=True)

# 6. Apply loads
dl = DistributedLoad(elem3, wx=0, wy=10)  # uniform distributed load on element 3
structure.add_load(dl)

# 7. Solve
structure.solve()

# 8. Print results
results = Results(structure)
disp = results.node_displacements()
reactions = results.reactions()
el_forces = results.element_forces()

# print(disp)
# print(reactions)  
print(disp["theta"][1])  #result = 1.2755102040816325e-05
print(el_forces["m_j"][2])  #result = 3.928571428571429
# print(el_forces) 

# drawer = DrawStructure(structure, scale=100)
# drawer.draw()
