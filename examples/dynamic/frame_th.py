#Value from SAP2000 is at top right node: max = 7.703153 and min = -7.677126

import numpy as np
from fem2d import Structure, Node, ElasticMaterial, BeamElement, Section, DrawStructure
from structdyn import MDF, GroundMotion, elcentro_chopra
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
        column_section.A,
        column_section.I,
    )
    structure.add_element(col)

# -------------------------------------------------------------------
# Supports: fixed at base (all DOFs fixed)
# -------------------------------------------------------------------
node_left_bottom.set_support(1, 1, 1)
node_right_bottom.set_support(1, 1, 1)

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

#elcentro ground motion

elc = elcentro_chopra()
gm = GroundMotion.from_arrays(elc["acc (g)"], dt= 0.02, scale_factor=g)

# influence vector
# in_vec = np.ones(mdf.n_dof)  # set to 1 for all stories for now

# run analysis
response = mdf.find_response_ground_motion(gm, method="newmark_beta")

print(response['u22'].max(), response['u22'].min())    #result from fem2d is max=7.703123133306655, min=-7.67719008289298

# drawer = DrawStructure(structure, scale=0.05)
# drawer.draw()
