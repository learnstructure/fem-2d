"""
Tests for 2D in-plane BeamElement roll_angle (0 and 180 degrees) in fem2d.

For roll_angle = 180 deg, local axes 2 and 3 are negated (keeping a right-handed
local frame). Global stiffness and mass matrices are unchanged; for global loads,
nodal displacements and support reactions are unchanged while internal V2 and M3
are negated; for local-y member loads, the global load direction reverses.
"""

from __future__ import annotations

import numpy as np
import pytest

from struct_core import (
    BeamElement,
    ConcreteMaterial,
    DistributedLoad,
    ElementPointLoad,
    LoadCase,
    Node,
    PointLoad,
    Project,
    RectangularSection,
    StructuralModel,
    Support,
)

from fem2d.adapters import model_from_core, model_to_core, result_to_core
from fem2d.elements.beam import BeamElement as FemBeamElement
from fem2d.materials.elastic import ElasticMaterial as FemElasticMaterial
from fem2d.nodes import Node as FemNode
from fem2d.sections.section import Section as FemSection


def _portal_frame(col_roll: float, *, shear: bool = False) -> Project:
    model = StructuralModel(
        dimension="2D",
        vertical_axis="Y",
        nodes=[
            Node(id="1", x=0.0, y=0.0),
            Node(id="2", x=0.0, y=4.0),
            Node(id="3", x=6.0, y=4.0),
            Node(id="4", x=6.0, y=0.0),
        ],
        materials=[ConcreteMaterial(id="C28", E=24.9e6, nu=0.2, fc=28.0e3, unit_weight=23.6)],
        sections=[
            RectangularSection(id="COL", material_id="C28", b=0.4, h=0.4),
            RectangularSection(id="BM", material_id="C28", b=0.3, h=0.6),
        ],
        elements=[
            BeamElement(
                id="C1", start_node="1", end_node="2", section_id="COL",
                include_shear_deformation=shear, roll_angle=col_roll,
            ),
            BeamElement(
                id="B1", start_node="2", end_node="3", section_id="BM",
                include_shear_deformation=shear, release_j=[False, False, True],
            ),
            BeamElement(
                id="C2", start_node="4", end_node="3", section_id="COL",
                include_shear_deformation=shear, roll_angle=col_roll,
            ),
        ],
        supports=[Support.fixed("1"), Support.pinned("4")],
        load_cases=[
            LoadCase(
                id="GLOBAL",
                include_self_weight=True,
                point_loads=[PointLoad(node_id="2", fx=15.0), PointLoad(node_id="3", mz=5.0)],
                element_loads=[
                    DistributedLoad(element_id="B1", wy=-20.0),
                    DistributedLoad(element_id="C1", wx=4.0, wx_end=2.0, start=0.1, end=0.9),
                ],
                element_point_loads=[
                    ElementPointLoad(element_id="C1", position=0.6, fx=6.0, mz=8.0),
                ],
            ),
            LoadCase(
                id="LOCAL_COL",
                element_loads=[DistributedLoad(element_id="C1", w2=3.0, is_local=True)],
            ),
        ],
    )
    return Project(model=model, units="kN, m")


@pytest.mark.parametrize("shear", [False, True])
def test_portal_frame_columns_0_vs_180_same_displacements_reactions_negated_vm(shear):
    s0 = model_from_core(_portal_frame(0.0, shear=shear), load_case_id="GLOBAL")
    s180 = model_from_core(_portal_frame(180.0, shear=shear), load_case_id="GLOBAL")

    # Global element stiffness matrices are identical
    for eid in ("C1", "B1", "C2"):
        np.testing.assert_allclose(
            s0.elements[eid].global_stiffness(),
            s180.elements[eid].global_stiffness(),
            atol=1e-9,
        )

    s0.solve()
    s180.solve()
    r0 = result_to_core(s0, analysis_case_id="GLOBAL")
    r180 = result_to_core(s180, analysis_case_id="GLOBAL")

    # Displacements and reactions are identical
    for nid in ("1", "2", "3", "4"):
        d0, d180 = r0.node_results[nid].displacement, r180.node_results[nid].displacement
        np.testing.assert_allclose([d0.ux, d0.uy, d0.rz], [d180.ux, d180.uy, d180.rz], atol=1e-12)
        rc0, rc180 = r0.node_results[nid].reaction, r180.node_results[nid].reaction
        np.testing.assert_allclose([rc0.fx, rc0.fy, rc0.mz], [rc180.fx, rc180.fy, rc180.mz], atol=1e-10)

    # Beam B1 internal forces are unchanged
    for sf0, sf180 in zip(r0.element_results["B1"].forces, r180.element_results["B1"].forces):
        np.testing.assert_allclose([sf0.P, sf0.V2, sf0.M3], [sf180.P, sf180.V2, sf180.M3], atol=1e-10)

    # Columns C1 and C2: P unchanged, V2 and M3 negated at every station
    for cid in ("C1", "C2"):
        for sf0, sf180 in zip(r0.element_results[cid].forces, r180.element_results[cid].forces):
            assert sf180.P == pytest.approx(sf0.P, abs=1e-10)
            assert sf180.V2 == pytest.approx(-sf0.V2, abs=1e-10)
            assert sf180.M3 == pytest.approx(-sf0.M3, abs=1e-10)


def test_cantilever_local_y_member_load_reverses_global_direction_for_180():
    """On a vertical column (0,0)->(0,4), local +2 is -X at roll=0 and +X at roll=180."""
    def _run(roll: float):
        model = StructuralModel(
            dimension="2D",
            vertical_axis="Y",
            nodes=[Node(id="1", x=0.0, y=0.0), Node(id="2", x=0.0, y=4.0)],
            materials=[ConcreteMaterial(id="M", E=25.0e6, nu=0.2, fc=30.0e3)],
            sections=[RectangularSection(id="S", material_id="M", b=0.4, h=0.4)],
            elements=[BeamElement(id="C1", start_node="1", end_node="2", section_id="S", roll_angle=roll)],
            supports=[Support.fixed("1")],
            load_cases=[
                LoadCase(
                    id="LOC",
                    element_loads=[DistributedLoad(element_id="C1", w2=3.0, is_local=True)],
                    element_point_loads=[ElementPointLoad(element_id="C1", position=0.5, fy=4.0, mz=2.0, is_local=True)],
                )
            ],
        )
        s = model_from_core(Project(model=model, units="kN, m"), load_case_id="LOC")
        s.solve()
        return result_to_core(s, analysis_case_id="LOC")

    r0 = _run(0.0)
    r180 = _run(180.0)

    # Global deflection and rotation reverse sign (load acts in -X for roll=0, +X for roll=180)
    assert r0.node_results["2"].displacement.ux < 0.0
    assert r180.node_results["2"].displacement.ux == pytest.approx(-r0.node_results["2"].displacement.ux, rel=1e-12)
    assert r180.node_results["2"].displacement.rz == pytest.approx(-r0.node_results["2"].displacement.rz, rel=1e-12)
    assert r180.node_results["1"].reaction.fx == pytest.approx(-r0.node_results["1"].reaction.fx, rel=1e-12)
    assert r180.node_results["1"].reaction.mz == pytest.approx(-r0.node_results["1"].reaction.mz, rel=1e-12)

    # Local V2 and M3 are identical in the element's own right-handed local frame
    for sf0, sf180 in zip(r0.element_results["C1"].forces, r180.element_results["C1"].forces):
        assert sf180.V2 == pytest.approx(sf0.V2, abs=1e-10)
        assert sf180.M3 == pytest.approx(sf0.M3, abs=1e-10)


def test_roll_angle_normalization_and_invalid_angle_refused():
    n1, n2 = FemNode(1, 0.0, 0.0), FemNode(2, 0.0, 4.0)
    sec = FemSection.from_rectangle(0.4, 0.4, FemElasticMaterial(E=25e6))

    assert FemBeamElement(1, n1, n2, section=sec, roll_angle=0.0).roll_angle == 0.0
    assert FemBeamElement(1, n1, n2, section=sec, roll_angle=360.0).roll_angle == 0.0
    assert FemBeamElement(1, n1, n2, section=sec, roll_angle=-180.0).roll_angle == 180.0
    assert FemBeamElement(1, n1, n2, section=sec, roll_angle=540.0).roll_angle == 180.0

    for bad in (45.0, 90.0, -90.0, 120.0):
        with pytest.raises(ValueError, match="roll_angle of 0 or 180"):
            FemBeamElement(1, n1, n2, section=sec, roll_angle=bad)


def test_model_to_core_preserves_roll_angle_and_mass_matrix_invariant():
    s0 = model_from_core(_portal_frame(0.0), load_case_id="GLOBAL")
    s180 = model_from_core(_portal_frame(180.0), load_case_id="GLOBAL")

    for eid in ("C1", "B1", "C2"):
        np.testing.assert_allclose(
            s0.elements[eid].mass_matrix(lumped=False),
            s180.elements[eid].mass_matrix(lumped=False),
            atol=1e-12,
        )
        np.testing.assert_allclose(
            s0.elements[eid].mass_matrix(lumped=True),
            s180.elements[eid].mass_matrix(lumped=True),
            atol=1e-12,
        )

    back = model_to_core(s180).model
    by_id = {e.id: e for e in back.elements}
    assert by_id["C1"].roll_angle == 180.0
    assert by_id["C2"].roll_angle == 180.0
    assert by_id["B1"].roll_angle == 0.0
