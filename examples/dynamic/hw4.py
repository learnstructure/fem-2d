# HW4 - Earthquake resistant design of structures
# first mode period by sap2000 is 0.740.
# To match with sap, shear stiffness modifier in sap needs to be very high to ignore shear deformation. also self weight needs to be zeroed out.
import numpy as np
from fem2d import Structure, Node, ElasticMaterial, BeamElement, Section, DrawStructure, plot_mode_shape
from structdyn.mdf.mdf import MDF
import math

# -------------------------------------------------------------------
# Model parameters (kips, inches, seconds)
# -------------------------------------------------------------------
g = 386.09  # in/s²
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
        beam_section,
        mass=0.0,  # no distributed mass
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
        column_section
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
        column_section
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
# Create MDF 
# -------------------------------------------------------------------
mdf = MDF.from_fem2d(structure)

# -------------------------------------------------------------------
# Eigenvalue analysis
# -------------------------------------------------------------------

omega, phi = mdf.modal.modal_analysis(dof_normalize=-1)
T = 2 * math.pi / omega

print("\nFirst 5 periods (seconds):")
for i in range(min(5, len(T))):
    print(f"Mode {i+1}: {T[i]:.6f} s")

# Register modal results on the structure so plot_mode_shape can use them directly.
structure.set_modal_results(omega, phi)

# Option 1 — convenience function (recommended for quick use):
plot_mode_shape(structure, mode=1)

# Option 2 — object-oriented style with DrawStructure (useful when you need the axes object
#             or want to customise the plot further, e.g. save to file, change figsize, etc.):
# drawer = DrawStructure(structure)
# drawer.draw_mode_shape(mode=1)
