import numpy as np
from fem2d import Structure, Node, ElasticMaterial, NewtonRaphsonSolver, TrussElementNL
from fem2d.results import Results

# -------------------------------------------------------------------
# Example 10.2 – Geometrically nonlinear analysis of a plane truss
# Units: kN, m
# -------------------------------------------------------------------

# 1. Create nodes
node1 = Node(1, 0.0, 0.0)
node2 = Node(2, 4.0, 3.0)
node3 = Node(3, 8.0, 0.0)

# 2. Material properties
E = 200e6  # 200 GPa = 200×10⁶ kN/m²
EA = 45155.0
A = EA / E

material = ElasticMaterial(E)

# 3. Create nonlinear truss elements
elem1 = TrussElementNL(1, node1, node2, material, A)
elem2 = TrussElementNL(2, node2, node3, material, A)
elem3 = TrussElementNL(3, node3, node1, material, A)

# 4. Build structure
structure = Structure()
structure.add_node(node1)
structure.add_node(node2)
structure.add_node(node3)
structure.add_element(elem1)
structure.add_element(elem2)
structure.add_element(elem3)

# 5. Apply supports – FIX ROTATIONS FOR ALL NODES!
node1.set_support(ux_fixed=True, uy_fixed=True)  # pinned
node3.set_support(ux_fixed=False, uy_fixed=True)  # roller

# 6. Apply external load (vertical -2000 kN at node 2)
node2.set_load(fx=0.0, fy=-2000.0, mz=0.0)

# 7. Run nonlinear analysis (one line!)
structure.solve_nonlinear(tolerance=1e-8, max_iter=30)

results = Results(structure)
print("Node Displacements:\n", results.node_displacements())
print("Reactions:\n", results.reactions())
print("Element Forces:\n", results.element_forces())


# # 8. Print results
# print("\nConverged displacements (free DOFs):")
# free_dofs = structure.free_dofs
# for i, dof in enumerate(free_dofs):
#     print(f"  DOF {dof+1}: {structure.disp[dof]:.6f} m")

# print("\nMember forces (axial, compression positive):")
# for el in structure.elements.values():
#     # el.update_state(structure.disp)   # already updated during solver
#     print(f"  Member {el.id}: Q = {el.Q:.1f} kN   (L_current = {el.L:.4f} m)")

# print("\nReactions:")
# for node in [node1, node3]:
#     print(
#         f"  Node {node.id}: Fx = {structure.reactions[node.dofs[0]]:8.2f} kN, Fy = {structure.reactions[node.dofs[1]]:8.2f} kN"
#     )
