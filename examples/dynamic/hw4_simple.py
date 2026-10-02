# HW4 - Earthquake resistant design of structures
# first mode period by sap2000 is 0.740.
# To match with sap, shear stiffness modifier in sap needs to be very high to ignore shear deformation. also self weight needs to be zeroed out.
import numpy as np
from fem2d.utils.simple_frame import SimpleFrame
from fem2d import Section, ElasticMaterial, DrawStructure, plot_mode_shape
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
# Build nodes and structure using SimpleFrame
# -------------------------------------------------------------------
frame = SimpleFrame()
nodes = []

# Bottom nodes (fixed base)
node_left_bottom = frame.add_node(1, 0.0, 0.0)
nodes.append(node_left_bottom)
node_right_bottom = frame.add_node(2, width, 0.0)
nodes.append(node_right_bottom)

# Story nodes
for story in range(1, n_stories + 1):
    y = story * story_height
    left_node = frame.add_node(2 * story + 1, 0.0, y)
    right_node = frame.add_node(2 * story + 2, width, y)
    nodes.append(left_node)
    nodes.append(right_node)

# -------------------------------------------------------------------
# Beams (horizontal)
# -------------------------------------------------------------------
for story in range(1, n_stories + 1):
    left_node = nodes[2 * story]     # zero‑based: story1 -> nodes[2] and [3]
    right_node = nodes[2 * story + 1]
    frame.add_frame(
        story,                       # element ID
        left_node.id,                # start node ID
        right_node.id,               # end node ID
        E,                           # Young's modulus
        beam_section.A,              # cross-sectional area
        beam_section.I               # moment of inertia
    )

# -------------------------------------------------------------------
# Columns (vertical)
# -------------------------------------------------------------------
# Left columns
for story in range(n_stories):
    bottom_node = nodes[2 * story]
    top_node = nodes[2 * (story + 1)]
    frame.add_frame(
        n_stories + story + 1,       # element ID
        bottom_node.id,              # start node ID
        top_node.id,                 # end node ID
        E,                           # Young's modulus
        column_section.A,            # cross-sectional area
        column_section.I             # moment of inertia
    )

# Right columns
for story in range(n_stories):
    bottom_node = nodes[2 * story + 1]
    top_node = nodes[2 * (story + 1) + 1]
    frame.add_frame(
        2 * n_stories + story + 1,   # element ID
        bottom_node.id,              # start node ID
        top_node.id,                 # end node ID
        E,                           # Young's modulus
        column_section.A,            # cross-sectional area
        column_section.I             # moment of inertia
    )

# -------------------------------------------------------------------
# Supports: fixed at base (all DOFs fixed)
# -------------------------------------------------------------------
frame.add_support(1, [1, 1, 1])  # left bottom node
frame.add_support(2, [1, 1, 1])  # right bottom node

# -------------------------------------------------------------------
# Lumped masses at each story node
# -------------------------------------------------------------------
# Each floor carries 450 kips, distributed equally to the two nodes
mass_per_node = 450.0 / 2.0 / g  # slugs
small_inertia = 1e-6

# Add masses to the underlying structure
for story in range(1, n_stories + 1):
    left_node_id = 2 * story + 1
    right_node_id = 2 * story + 2
    frame.structure.nodes[left_node_id].set_mass(mass=mass_per_node, inertia=small_inertia)
    frame.structure.nodes[right_node_id].set_mass(mass=mass_per_node, inertia=small_inertia)

# -------------------------------------------------------------------
# Dynamic analysis using MDF (same as original)
# -------------------------------------------------------------------
mdf = MDF.from_fem2d(frame.structure)

# -------------------------------------------------------------------
# Eigenvalue analysis
# -------------------------------------------------------------------
import math

omega, phi = mdf.modal.modal_analysis(dof_normalize=-1)
T = 2 * math.pi / omega

print("\nFirst 5 periods (seconds):")
for i in range(min(5, len(T))):
    print(f"Mode {i+1}: {T[i]:.6f} s")

# Register modal results on the structure so plot_mode_shape can use them directly.
frame.structure.set_modal_results(omega, phi)

# Option 1 — convenience function (recommended for quick use):
plot_mode_shape(frame, mode=1)

# Option 2 — object-oriented style with DrawStructure (useful when you need the axes object
#             or want to customise the plot further, e.g. save to file, change figsize, etc.):
# drawer = DrawStructure(frame.structure)
# drawer.draw_mode_shape(mode=1)