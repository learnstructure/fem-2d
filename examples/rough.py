import numpy as np
from fem2d import Structure, Node, ElasticMaterial, Section, BeamElement
from fem2d.loads import DistributedLoad
from fem2d.results import Results

# -------------------------------------------------------------------
# Units: kN, m
# -------------------------------------------------------------------

# Create structure
structure = Structure()

# Add nodes
node1 = Node(1, 0.0, 0.0)
node2 = Node(2, 0.0, 6.0)
node3 = Node(3, 0.0, 12.0)
node4 = Node(4, 9.0, 6.0)
node5 = Node(5, 9.0, 0.0)

structure.add_node(node1)
structure.add_node(node2)
structure.add_node(node3)
structure.add_node(node4)
structure.add_node(node5)

# Material and section properties
E = 30e6
A = 75e-3
I = 4.8e-4

material = ElasticMaterial(E)
section = Section(A, I)

# Create beam elements
elem1 = BeamElement(1, node1, node2, material, section.A, section.I)
elem2 = BeamElement(2, node2, node3, material, section.A, section.I)
elem3 = BeamElement(3, node3, node4, material, section.A, section.I)
elem4 = BeamElement(4, node2, node4, material, section.A, section.I)
elem5 = BeamElement(5, node4, node5, material, section.A, section.I)

structure.add_element(elem1)
structure.add_element(elem2)
structure.add_element(elem3)
structure.add_element(elem4)
structure.add_element(elem5)

# Supports: nodes 1 and 5 fully fixed
node1.set_support(ux_fixed=True, uy_fixed=True, rz_fixed=True)
node5.set_support(ux_fixed=True, uy_fixed=True, rz_fixed=True)

# Nodal loads: horizontal forces at nodes 2 and 3
node2.set_load(fx=80.0, fy=0.0, mz=0.0)
node3.set_load(fx=40.0, fy=0.0, mz=0.0)

# Distributed load on element 3 (from node 3 to node 4)
udl = DistributedLoad(elem3, wy=12)
structure.add_load(udl)

# Solve
structure.solve()

# Post‑processing
results = Results(structure)
disp = results.node_displacements()
reactions = results.reactions()
el_forces = results.element_forces()

print("Node Displacements:\n", disp)
print("\nReactions:\n", reactions)
print("\nElement Forces:\n", el_forces)

# Check specific values
print("\nNode 3 rotation:", disp["theta"][2])
print("Element 3 moment at i‑end:", el_forces["m_i"][2])
