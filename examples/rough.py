# HW4 - Earthquake resistant design of structures

import numpy as np
import scipy.linalg as la
from fem2d import Structure, Node, ElasticMaterial, BeamElement

# units in kips-in
# -------------------------------------------------------------------
# Model parameters
# -------------------------------------------------------------------
width = 16 * 12  # bay width
story_height = 10 * 12  # story height
n_stories = 4

# Material properties (steel)
E = 6000

# Section properties (example: HEA 300)
A_beam = 20 * 35
I_beam = 20 * 35**3 / 12

A_column = 20 * 20
I_column = 20 * 20**3 / 12

# Extra mass on beams
beam_extra_mass = 450 / (width * 386.4)  # 450 kips is total weight lumped at each floor

# -------------------------------------------------------------------
# Create nodes
# -------------------------------------------------------------------
nodes = []
structure = Structure()

# Bottom nodes (ground level)
node_left_bottom = Node(1, 0.0, 0.0)
node_right_bottom = Node(2, width, 0.0)
structure.add_node(node_left_bottom)
structure.add_node(node_right_bottom)
nodes.append(node_left_bottom)
nodes.append(node_right_bottom)

# Create nodes for each story
for story in range(1, n_stories + 1):
    y = story * story_height
    node_left = Node(2 * story + 1, 0.0, y)
    node_right = Node(2 * story + 2, width, y)
    structure.add_node(node_left)
    structure.add_node(node_right)
    nodes.append(node_left)
    nodes.append(node_right)

# -------------------------------------------------------------------
# Material
# -------------------------------------------------------------------
material = ElasticMaterial(E)

# -------------------------------------------------------------------
# Elements: beams and columns
# -------------------------------------------------------------------
# Beams (horizontal) – one per story, connecting left and right nodes at same height
for story in range(1, n_stories + 1):
    left_node = nodes[2 * story]  # zero‑based: story1 -> nodes[2] and [3]
    right_node = nodes[2 * story + 1]
    beam = BeamElement(
        story,
        left_node,
        right_node,
        material,
        A_beam,
        I_beam,
        extra_mass=beam_extra_mass,
    )
    structure.add_element(beam)

# Columns (vertical) – left and right columns
# Left column: connect nodes at consecutive stories
for story in range(n_stories):
    bottom_node = nodes[2 * story]  # story 0: node 0 (index 0) – bottom left
    top_node = nodes[2 * (story + 1)]  # story 1: node 2
    col = BeamElement(
        n_stories + story + 1, bottom_node, top_node, material, A_column, I_column
    )
    structure.add_element(col)

# Right column
for story in range(n_stories):
    bottom_node = nodes[2 * story + 1]  # bottom right
    top_node = nodes[2 * (story + 1) + 1]  # top right
    col = BeamElement(
        2 * n_stories + story + 1, bottom_node, top_node, material, A_column, I_column
    )
    structure.add_element(col)

# -------------------------------------------------------------------
# Supports: pin at bottom nodes (ux, uy fixed, rz free)
# For a pinned base, we fix translations but allow rotation.
# If you prefer fixed base, set rz_fixed=True as well.
# -------------------------------------------------------------------
node_left_bottom.set_support(ux_fixed=True, uy_fixed=True, rz_fixed=True)
node_right_bottom.set_support(ux_fixed=True, uy_fixed=True, rz_fixed=True)

# -------------------------------------------------------------------
# Prepare matrices for eigenvalue analysis
# -------------------------------------------------------------------
structure.number_dofs()
structure.apply_boundary_conditions()
structure.assemble_stiffness()
structure.assemble_mass_matrix()

# Reduce to free degrees of freedom
K_ff, M_ff = structure.get_reduced_matrices()

print("Free DOFs:", structure.free_dofs)
print("K_ff shape:", K_ff.shape)
print("M_ff shape:", M_ff.shape)

# -------------------------------------------------------------------
# Solve eigenvalue problem: (K_ff - ω² M_ff) φ = 0
# -------------------------------------------------------------------
# Using scipy.linalg.eig (generalized eigenvalue problem)
eigvals, eigvecs = la.eig(K_ff, M_ff)

# Extract real parts (should be positive real)
omega = np.sqrt(np.real(eigvals))
freqs = omega / (2 * np.pi)  # Hz

# Sort in increasing order
idx = np.argsort(omega)
omega = omega[idx]
freqs = freqs[idx]
eigvecs = eigvecs[:, idx]
T = 1 / freqs

# Print first 5 frequencies
# print("\nNatural frequencies (Hz):")
# for i in range(min(5, len(freqs))):
#     print(f"{i+1}: {freqs[i]:.3f} Hz")

# Print first 5 time periods
print("\n Time periods:")
for i in range(min(5, len(T))):
    print(f"{i+1}: {T[i]:.3f} s")

# Optional: print mode shapes (for verification)
# print("\nFirst mode shape (free DOFs):\n", eigvecs[:, 0])
