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


# 1. Create nodes
node1 = Node(1, 0, 0)
node2 = Node(2, 2, 0)


# 2. Material and section properties
E = 200e6  # kN/m^2
material = ElasticMaterial(E)  # same material for all elements
section = Section.from_rectangle(width=0.1, depth=0.2)

# 3. Create elements
elem1 = BeamElement(1, node1, node2, material, section.A, section.I)

# 4. Create structure and add nodes/elements
structure = Structure()
structure.add_node(node1)
structure.add_node(node2)
structure.add_element(elem1)


# 5. Apply supports (all DOFs fixed)
node1.set_support(ux_fixed=True, uy_fixed=True, rz_fixed=True)


# 6. Apply loads
node2.set_load(fx=0.0, fy=-800, mz=0.0)

# 7. Solve
structure.solve()

# 8. Print results
results = Results(structure)
disp = results.node_displacements()

print(disp)

