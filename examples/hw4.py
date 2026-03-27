# HW4 - Earthquake resistant design of structures
import numpy as np
from fem2d import Structure, Node, ElasticMaterial, BeamElement, Section, DrawStructure
from structdyn.mdf.mdf import MDF
import math

# -------------------------------------------------------------------
# Model parameters (kips, inches, seconds)
# -------------------------------------------------------------------
g = 386.4  # in/s²
width = 16 * 12  # in
story_height = 10 * 12  # in
n_stories = 4

# Concrete properties
E = 6000  # ksi
material = ElasticMaterial(E)

# Section properties
# Beam: 20 in x 35 in
beam_section = Section.from_rectangle(width=20, depth=35)

# Column: 20 in x 20 in
column_section = Section.from_rectangle(width=20, depth=20)

# -------------------------------------------------------------------
# Build nodes
# -------------------------------------------------------------------
nodes = []
structure = Structure()

# Bottom nodes (fixed base)
node_left_bottom = Node(1, 0.0, 0.0)
node_right_bottom = Node(2, width, 0.0)
structure.add_node(node_left_bottom)
structure.add_node(node_right_bottom)
nodes.append(node_left_bottom)
nodes.append(node_right_bottom)

# Story nodes
for story in range(1, n_stories + 1):
    y = story * story_height
    left_node = Node(2 * story + 1, 0.0, y)
    right_node = Node(2 * story + 2, width, y)
    structure.add_node(left_node)
    structure.add_node(right_node)
    nodes.append(left_node)
    nodes.append(right_node)


# -------------------------------------------------------------------
# Beams (horizontal)
# -------------------------------------------------------------------
for story in range(1, n_stories + 1):
    left_node = nodes[2 * story]  # zero‑based: story1 -> nodes[2] and [3]
    right_node = nodes[2 * story + 1]
    beam = BeamElement(
        story,
        left_node,
        right_node,
        material,
        beam_section.A,
        beam_section.I,
        extra_mass=0.0,  # no distributed mass
    )
    structure.add_element(beam)

# -------------------------------------------------------------------
# Columns (vertical)
# -------------------------------------------------------------------
# Left columns
for story in range(n_stories):
    bottom_node = nodes[2 * story]
    top_node = nodes[2 * (story + 1)]
    col = BeamElement(
        n_stories + story + 1,
        bottom_node,
        top_node,
        material,
        column_section.A,
        column_section.I,
    )
    structure.add_element(col)

# Right columns
for story in range(n_stories):
    bottom_node = nodes[2 * story + 1]
    top_node = nodes[2 * (story + 1) + 1]
    col = BeamElement(
        2 * n_stories + story + 1,
        bottom_node,
        top_node,
        material,
        column_section.A,
        column_section.I,
    )
    structure.add_element(col)

# -------------------------------------------------------------------
# Supports: fixed at base (all DOFs fixed)
# -------------------------------------------------------------------
node_left_bottom.set_support(ux_fixed=True, uy_fixed=True, rz_fixed=True)
node_right_bottom.set_support(ux_fixed=True, uy_fixed=True, rz_fixed=True)

# -------------------------------------------------------------------
# Lumped masses at each story node
# -------------------------------------------------------------------
# Each floor carries 450 kips, distributed equally to the two nodes
mass_per_node = 450.0 / 2.0 / g  # slugs
small_inertia = 1e-6

for story in range(1, n_stories + 1):
    left_node = nodes[2 * story]  # left node of this story
    right_node = nodes[2 * story + 1]  # right node of this story
    left_node.set_mass(mass=mass_per_node, inertia=small_inertia)
    right_node.set_mass(mass=mass_per_node, inertia=small_inertia)

# -------------------------------------------------------------------
# Assemble matrices and reduce
# -------------------------------------------------------------------
structure.number_dofs()
structure.apply_boundary_conditions()
structure.assemble_stiffness()
# structure.solve()
structure.assemble_mass_matrix()

K_ff, M_ff = structure.get_reduced_matrices()

print("K_ff shape:", K_ff.shape)
print("M_ff shape:", M_ff.shape)

# -------------------------------------------------------------------
# Eigenvalue analysis
# -------------------------------------------------------------------
mdf = MDF(M_ff, K_ff)
omega, phi = mdf.modal.modal_analysis(dof_normalize=-1)
T = 2 * math.pi / omega

print("\nFirst 5 periods (seconds):")
for i in range(min(5, len(T))):
    print(f"Mode {i+1}: {T[i]:.3f} s")

# Choose which mode to plot (first mode: index 0)
mode_index = 0
mode_shape_reduced = phi[:, mode_index]

# Build full displacement vector for the structure
disp_full = np.zeros(structure.neq)  # all DOFs initially zero
disp_full[structure.free_dofs] = mode_shape_reduced  # fill free DOFs with mode shape

structure.disp = disp_full
drawer = DrawStructure(structure, scale=0.05)
drawer.draw()
