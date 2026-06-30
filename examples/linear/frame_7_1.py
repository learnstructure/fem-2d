# Example 7.1: Kassimali A. (2022), Matrix analysis of structures
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
node2 = Node(2, 0, 5)
node3 = Node(3, 5, 5)
node4 = Node(4, 5, 0)

# 2. Material and section properties
E = 200e6  # kN/m^2
A = 6500e-6  # m^2
I = 150e-6  # m^4
material = ElasticMaterial(E)  # same material for all elements
section = Section(A, I)

# 3. Create elements
# elem1 = BeamElement(1, node1, node2, material, section.A, section.I)
# elem2 = BeamWithHingesElement(2, node2, node3, material, section.A, section.I, hinge_i=True, hinge_j=False)
elem1 = BeamWithHingesElement(1, node1, node2, material, section.A, section.I, hinge_i=False, hinge_j=True)
elem2 = BeamElement(2, node2, node3, material, section.A, section.I)
elem3 = BeamElement(3, node3, node4, material, section.A, section.I)

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
node4.set_support(ux_fixed=True, uy_fixed=True, rz_fixed=False)

# 6. Apply loads
# node2.set_load(fx=90.75, fy=0, mz=0)  # horizontal point load at node 2
pl1 = PointLoad(node2, fx=90.75, fy=0, mz=0)
structure.add_load(pl1)

pl2 = ElementPointLoad(elem2, x=2.5, py = -300)  # point load at midspan of element 2
structure.add_load(pl2)

dl = DistributedLoad(elem1, wx=0, wy=-19.2)  # uniform distributed load on element 2
structure.add_load(dl)

# 7. Solve
structure.solve()

# 8. Print results
results = Results(structure)
disp = results.node_displacements()
reactions = results.reactions()
el_forces = results.element_forces()

print(disp)  
print(el_forces)  
print(disp["theta"][2])     # result = -0.0013662379393815254
print(el_forces["m_i"][2])      # result = 304.15771709113193

drawer = DrawStructure(structure, scale=10)
drawer.draw(show_deformed=False)
