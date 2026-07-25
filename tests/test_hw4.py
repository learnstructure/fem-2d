"""Test the dynamic example ``hw4.py`` by reproducing its calculations.

The original example builds a 4‑story frame, assembles stiffness and mass
matrices, performs a modal analysis, and prints the first five natural
periods.  We recreate the same steps here and assert that the first period
matches the reference value from SAP2000 (≈ 0.740 s).
"""

import math

import numpy as np
import pytest

from fem2d import (
    Structure,
    Node,
    ElasticMaterial,
    BeamElement,
    Section,
)
from structdyn.mdf.mdf import MDF


def build_hw4_structure():
    """Construct the structure defined in ``examples/dynamic/hw4.py``.

    Returns:
        tuple: ``(structure, K_ff, M_ff)`` where ``K_ff`` and ``M_ff`` are the
        reduced stiffness and mass matrices ready for modal analysis.
    """
    # Model parameters (kips, inches, seconds)
    g = 386.09  # in/s²
    width = 16 * 12  # in
    story_height = 10 * 12  # in
    n_stories = 4

    # Material and sections
    material = ElasticMaterial(6000)  # ksi
    beam_section = Section.from_rectangle(width=20, depth=35)
    column_section = Section.from_rectangle(width=20, depth=20)

    # Build nodes
    structure = Structure()
    nodes = []
    node_left_bottom = Node(1, 0.0, 0.0)
    node_right_bottom = Node(2, width, 0.0)
    structure.add_node(node_left_bottom)
    structure.add_node(node_right_bottom)
    nodes.extend([node_left_bottom, node_right_bottom])

    for story in range(1, n_stories + 1):
        y = story * story_height
        left_node = Node(2 * story + 1, 0.0, y)
        right_node = Node(2 * story + 2, width, y)
        structure.add_node(left_node)
        structure.add_node(right_node)
        nodes.extend([left_node, right_node])

    # Beams (horizontal)
    for story in range(1, n_stories + 1):
        left_node = nodes[2 * story]
        right_node = nodes[2 * story + 1]
        beam = BeamElement(
            story,
            left_node,
            right_node,
            material,
            beam_section,
            extra_mass=0.0,
        )
        structure.add_element(beam)

    # Columns (vertical)
    for story in range(n_stories):
        bottom_node = nodes[2 * story]
        top_node = nodes[2 * (story + 1)]
        col = BeamElement(
            n_stories + story + 1,
            bottom_node,
            top_node,
            material,
            column_section,
        )
        structure.add_element(col)

    for story in range(n_stories):
        bottom_node = nodes[2 * story + 1]
        top_node = nodes[2 * (story + 1) + 1]
        # The original script passed ``column_section.A`` and ``column_section.I``
        # separately; BeamElement accepts a ``Section`` instance, so we reuse the
        # ``column_section`` directly for consistency.
        col = BeamElement(
            2 * n_stories + story + 1,
            bottom_node,
            top_node,
            material,
            column_section,
        )
        structure.add_element(col)

    # Supports at the base
    node_left_bottom.set_support(ux_fixed=True, uy_fixed=True, rz_fixed=True)
    node_right_bottom.set_support(ux_fixed=True, uy_fixed=True, rz_fixed=True)

    # Lumped masses at each story node (450 kips total per floor)
    mass_per_node = 450.0 / 2.0 / g  # slugs per node
    small_inertia = 1e-6
    for story in range(1, n_stories + 1):
        left_node = nodes[2 * story]
        right_node = nodes[2 * story + 1]
        left_node.set_mass(mass=mass_per_node, inertia=small_inertia)
        right_node.set_mass(mass=mass_per_node, inertia=small_inertia)

    # Assemble matrices and obtain reduced forms
    structure.number_dofs()
    structure.apply_boundary_conditions()
    structure.assemble_stiffness()
    structure.assemble_mass_matrix()
    K_ff, M_ff = structure.get_reduced_matrices()
    return structure, K_ff, M_ff


def test_hw4_modal_analysis():
    """Validate the first natural period of the 4‑story frame.

    The reference period from SAP2000 is 0.740 s.  We assert that the computed
    period is within 5 % of that value.
    """
    _, K_ff, M_ff = build_hw4_structure()
    mdf = MDF(M_ff, K_ff)
    omega, _ = mdf.modal.modal_analysis(dof_normalize=-1)
    periods = 2 * math.pi / omega
    first_period = periods[0]
    # Allow a modest tolerance (1 %) because of modeling differences
    assert math.isclose(first_period, 0.740, rel_tol=0.01), (
        f"First period {first_period:.6f}s deviates from expected 0.740s"
    )
