# Example 10.2: Kassimali A. (2022), Matrix analysis of structures
# units in kN, m
import numpy as np
from fem2d import Structure, Node, ElasticMaterial, TrussElement
from fem2d.results import Results
from fem2d.sections import Section

# 1. Create nodes
node1 = Node(1, 0.0, 0.0)
node2 = Node(2, 4.0, 3.0)
node3 = Node(3, 8.0, 0.0)

# 2. Material properties
E = 200e6  # 200 GPa = 200e6 kN/m^2
EA = 45155.0  # axial stiffness [kN]
A = EA / E  # back-computed area [m^2]

material = ElasticMaterial(E)
section = Section(A)

# 3. Create truss elements.
elem1 = TrussElement(1, node1, node2, material, section)
elem2 = TrussElement(2, node2, node3, material, section)
elem3 = TrussElement(3, node3, node1, material, section)

# 4. Build structure
structure = Structure()
structure.add_node(node1)
structure.add_node(node2)
structure.add_node(node3)
structure.add_element(elem1)
structure.add_element(elem2)
structure.add_element(elem3)

# 5. Apply supports
node1.set_support(1, 1)  # pinned (ux, uy fixed)
node3.set_support(0, 1)  # roller (only uy fixed)

# 6. Apply external load (vertical -2000 kN at node 2)
node2.set_load(fx=0.0, fy=-2000.0, mz=0.0)

# 7. Run geometrically-nonlinear analysis.
structure.solve(geometric_nonlinear=True, tolerance=1e-8, max_iter=30)

# 8. Post-process
results = Results(structure)
disp = results.node_displacements()
reactions = results.reactions()
el_forces = results.element_forces()

print("Node Displacements:\n", disp)
print("Reactions:\n", reactions)
print("Element Forces:\n", el_forces)

print(disp["uy"][1])  # result = -0.6499578914694384
print(el_forces["fx_i"][2])  # result = -1768.7701456115713
