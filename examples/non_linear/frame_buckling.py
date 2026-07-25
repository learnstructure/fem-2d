import numpy as np
from fem2d import Structure, Node, ElasticMaterial, BeamElement, TrussElement, buckling_analysis, Results

# 1. Create nodes
node1 = Node(1, 0.0, 0.0)
node2 = Node(2, 4.0, 3.0)
node3 = Node(3, 8.0, 0.0)

# 2. Material properties
E = 70e6  
A = 645.6e-6
I = 0.379389e-6  # m^4

material = ElasticMaterial(E)
from fem2d.sections import Section
section = Section(A, I)

# 3. Create beam elements
elem1 = BeamElement(1, node1, node2, material, section)
elem2 = BeamElement(2, node2, node3, material, section)
elem3 = BeamElement(3, node3, node1, material, section)

# 4. Build structure
structure = Structure()
structure.add_node(node1)
structure.add_node(node2)
structure.add_node(node3)
structure.add_element(elem1)
structure.add_element(elem2)
structure.add_element(elem3)

# 5. Apply supports 
node1.set_support(ux_fixed=True, uy_fixed=True)  # pinned
node3.set_support(ux_fixed=False, uy_fixed=True)  # roller

# 6. Apply external load
node2.set_load(fx=0.0, fy=-10.0, mz=0.0)

# 7. Run anaysis
structure.solve()

# results = Results(structure)
# disp = results.node_displacements()

# 8. Perform buckling analysis
factors, modes = buckling_analysis(structure, num_modes=2)

print("Buckling Factors:\n", factors)       # result = [ 2.87630061 -8.3502087 ]
# print("Buckling Modes:\n", modes)
