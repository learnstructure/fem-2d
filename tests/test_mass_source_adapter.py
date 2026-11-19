"""Regression tests for MassSource nodal inertia preservation and include_element_mass in fem2d."""

from __future__ import annotations

import math

import numpy as np
import pytest
from struct_core import (
    BeamElement,
    DistributedLoad,
    ElasticMaterial,
    LoadCase,
    MassSource,
    MassSourceFactor,
    Node,
    PointLoad,
    Project,
    RectangularSection,
    StructuralModel,
    Support,
)

from fem2d.adapters import model_from_core


def _cantilever_2d(
    *,
    node_mass: float = 0.0,
    node_inertia: float = 0.0,
    mass_source: MassSource,
    stiffness_modifiers: dict | None = None,
) -> Project:
    model = StructuralModel(
        dimension="2D",
        vertical_axis="Y",
        units="kN, m",
        nodes=[
            Node(id="1", x=0.0, y=0.0),
            Node(id="2", x=0.0, y=4.0, mass=node_mass, inertia=node_inertia),
        ],
        supports=[Support.fixed("1")],
        materials=[ElasticMaterial(id="M", E=30.0e6, nu=0.2, unit_weight=24.0)],
        sections=[RectangularSection(id="C", material_id="M", b=0.4, h=0.4)],
        elements=[
            BeamElement(
                id="E1",
                start_node="1",
                end_node="2",
                section_id="C",
                stiffness_modifiers=stiffness_modifiers or {},
            )
        ],
        load_cases=[
            LoadCase(
                id="SDL",
                point_loads=[PointLoad(node_id="2", fy=-98.0665)],  # exactly 10.0 tonnes (kN*s^2/m)
                element_loads=[DistributedLoad(element_id="E1", wy=-9.80665)],  # 1.0 tonne/m * 4 m = 4.0 tonnes
            )
        ],
        mass_source=mass_source,
    )
    return Project(model=model, units="kN, m")


def test_mass_source_preserves_nodal_inertia_and_counts_added_mass_once():
    # Node 2: added mass = 5000 kg (5.0 tonnes), added inertia = 1200 kg*m^2 (1.2 tonne*m^2)
    # Plus point load SDL = 98.0665 kN -> 10.0 tonnes -> total nodal mass = 15.0 tonnes
    ms = MassSource(
        id="MS1",
        include_element_mass=False,
        include_added_mass=True,
        factors=[MassSourceFactor(load_case_id="SDL", factor=1.0)],
    )
    proj = _cantilever_2d(node_mass=5000.0, node_inertia=1200.0, mass_source=ms)
    st = model_from_core(proj)

    n2 = st.nodes["2"]
    assert n2.mass == pytest.approx(15.0, rel=1e-12)
    assert n2.inertia == pytest.approx(1.2, rel=1e-12)

    st.number_dofs()
    st.apply_boundary_conditions()
    _, M_ff = st.get_reduced_matrices(lumped=True, include_rotational_inertia=False)
    # Free DOFs at node 2: [ux, uy, rz]. Element has 4.0 tonnes distributed SDL -> 2.0 tonnes at node 2.
    assert np.diag(M_ff) == pytest.approx([17.0, 17.0, 1.2], rel=1e-12)


def test_mass_source_exclude_added_mass_excludes_nodal_mass_and_inertia():
    ms = MassSource(
        id="MS2",
        include_element_mass=False,
        include_added_mass=False,
        factors=[MassSourceFactor(load_case_id="SDL", factor=1.0)],
    )
    proj = _cantilever_2d(node_mass=5000.0, node_inertia=1200.0, mass_source=ms)
    st = model_from_core(proj)

    n2 = st.nodes["2"]
    assert n2.mass == pytest.approx(10.0, rel=1e-12)
    assert n2.inertia == pytest.approx(0.0, abs=1e-15)

    st.number_dofs()
    st.apply_boundary_conditions()
    _, M_ff = st.get_reduced_matrices(lumped=True, include_rotational_inertia=False)
    # 10.0 tonnes point load + 2.0 tonnes half of distributed SDL; zero inertia
    assert np.diag(M_ff) == pytest.approx([12.0, 12.0, 0.0], rel=1e-12)


def test_include_element_mass_false_removes_frame_self_mass_even_with_modifiers():
    # Massless column (include_element_mass=False) with tip mass 10000 kg (10 tonnes), no load factors
    ms = MassSource(id="MS3", include_element_mass=False, include_added_mass=True, factors=[])
    proj = _cantilever_2d(
        node_mass=10000.0,
        node_inertia=0.0,
        mass_source=ms,
        stiffness_modifiers={"A": 0.6, "mass": 2.0},
    )
    st = model_from_core(proj)
    st.number_dofs()
    st.apply_boundary_conditions()
    K_ff, M_ff = st.get_reduced_matrices(lumped=True, include_rotational_inertia=False)

    # Zero frame self-mass: only the 10.0 tonne nodal mass enters M_ff
    assert np.diag(M_ff) == pytest.approx([10.0, 10.0, 0.0], rel=1e-12)

    # Verify first lateral frequency matches closed-form cantilever omega = sqrt(3*E*I / (m*L^3))
    E, I, L, m = 30.0e6, 0.4**4 / 12.0, 4.0, 10.0
    omega_exact = math.sqrt(3.0 * E * I / (m * L**3))
    from structdyn import MDF

    omega, _ = MDF.from_fem(st, lumped=True).modal.modal_analysis(n_modes=1)
    assert float(omega[0]) == pytest.approx(omega_exact, rel=1e-12)


def test_include_element_mass_true_includes_self_mass_once():
    ms = MassSource(id="MS4", include_element_mass=True, include_added_mass=True, factors=[])
    proj = _cantilever_2d(node_mass=5000.0, node_inertia=1200.0, mass_source=ms)
    st = model_from_core(proj)
    st.number_dofs()
    st.apply_boundary_conditions()
    _, M_ff = st.get_reduced_matrices(lumped=True, include_rotational_inertia=False)

    # Half of element self-mass at node 2: 0.5 * (24.0 / 9.80665) * 0.16 * 4.0
    m_self_half = 0.5 * (24.0 / 9.80665) * 0.16 * 4.0
    assert np.diag(M_ff) == pytest.approx([5.0 + m_self_half, 5.0 + m_self_half, 1.2], rel=1e-12)
