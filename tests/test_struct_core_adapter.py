"""
Unit tests for fem2d.adapters.struct_core_adapter.

Covers:

- :func:`model_from_core` -- struct_core -> fem2d.Structure, accepting
  either a bare ``StructuralModel`` or a ``Project`` wrapper.
- :func:`model_to_core` -- fem2d.Structure -> struct_core.StructuralModel,
  accepting ``Structure`` / ``SimpleFrame`` / ``Results`` and rejecting
  foreign types. Includes a full round-trip (fem2d -> core -> fem2d) check.
- :func:`result_to_core` -- end-to-end on a cantilever beam, with reaction
  suppression, envelope checks, and ``AnalysisResult`` schema validation.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from struct_core import (
    AnalysisResult,
    BeamElement,
    ElasticMaterial,
    LinearStaticAnalysis,
    LoadCase,
    Node,
    PointLoad,
    Project,
    RectangularSection,
    StructuralModel,
    Support,
    validate_structural_model,
)

from fem2d import SimpleFrame
from fem2d.adapters import model_from_core, model_to_core, result_to_core
from fem2d.results import Results
from fem2d.structure import Structure


# ---------------------------------------------------------------------------
# Shared helpers / fixtures
# ---------------------------------------------------------------------------


def _build_solved_beam() -> Structure:
    """Build a 2 m cantilever beam with a tip point load and solve it."""
    frame = SimpleFrame()
    frame.add_node(1, 0.0, 0.0)
    frame.add_node(2, 2.0, 0.0)
    frame.add_frame(1, 1, 2, 30000.0, 10.0, 1e-4)
    frame.add_support(1, [True, True, True])
    frame.add_node_load(2, [0.0, -5000.0, 0.0])
    frame.solve()
    return frame.structure


@pytest.fixture()
def cantilever_results():
    E_val = 200e9
    b_val, h_val = 0.2, 0.4
    I_val = (b_val * (h_val**3)) / 12.0
    L_val = 4.0
    P_val = -10000.0

    model = StructuralModel(
        nodes=[Node(id=1, x=0.0, y=0.0), Node(id=2, x=L_val, y=0.0)],
        materials=[ElasticMaterial(id="Steel", E=E_val, rho=7850.0)],
        sections=[RectangularSection(id="BeamSection", b=b_val, h=h_val)],
        elements=[
            BeamElement(
                id=1,
                start_node=1,
                end_node=2,
                material_id="Steel",
                section_id="BeamSection",
            )
        ],
        supports=[Support.fixed(node_id=1)],
        load_cases=[
            LoadCase(id="TipLoad", point_loads=[PointLoad(node_id=2, fy=P_val)])
        ],
        analysis_cases=[
            LinearStaticAnalysis(id="StaticSolve", load_case_id="TipLoad")
        ],
    )
    project = Project()
    project.model = model
    project.metadata.title = "Cantilever Adapter Test"
    structure = model_from_core(project)
    structure.solve()
    return Results(structure), {"E": E_val, "I": I_val, "L": L_val, "P": P_val}


# ===========================================================================
# model_from_core
# ===========================================================================


def test_model_from_core_builds_and_solves_cantilever():
    E_val = 200e9
    b_val, h_val = 0.2, 0.4
    I_val = (b_val * (h_val**3)) / 12.0
    L_val = 4.0
    P_val = -10000.0

    model = StructuralModel(
        nodes=[Node(id=1, x=0.0, y=0.0), Node(id=2, x=L_val, y=0.0)],
        materials=[ElasticMaterial(id="Steel", E=E_val, rho=7850.0)],
        sections=[RectangularSection(id="BeamSection", b=b_val, h=h_val)],
        elements=[
            BeamElement(
                id=1,
                start_node=1,
                end_node=2,
                material_id="Steel",
                section_id="BeamSection",
            )
        ],
        supports=[Support.fixed(node_id=1)],
        load_cases=[
            LoadCase(id="TipLoad", point_loads=[PointLoad(node_id=2, fy=P_val)])
        ],
        analysis_cases=[
            LinearStaticAnalysis(id="StaticSolve", load_case_id="TipLoad")
        ],
    )
    project = Project()
    project.model = model
    project.metadata.title = "Cantilever Beam Test"

    structure = model_from_core(project)
    assert len(structure.nodes) == 2
    assert len(structure.elements) == 1

    structure.solve()
    analytical_delta = (P_val * (L_val**3)) / (3.0 * E_val * I_val)
    np.testing.assert_allclose(structure.disp[4], analytical_delta, rtol=1e-4)


def test_model_from_core_accepts_structural_model_directly():
    model = StructuralModel(
        nodes=[Node(id=1, x=0.0, y=0.0), Node(id=2, x=1.0, y=0.0)],
        materials=[ElasticMaterial(id="M", E=200e9, rho=0.0)],
        sections=[RectangularSection(id="S", b=0.1, h=0.2)],
        elements=[
            BeamElement(id=1, start_node=1, end_node=2, material_id="M", section_id="S")
        ],
        supports=[Support.fixed(node_id=1)],
        load_cases=[LoadCase(id="L", point_loads=[PointLoad(node_id=2, fy=-1000.0)])],
        analysis_cases=[LinearStaticAnalysis(id="A", load_case_id="L")],
    )
    structure = model_from_core(model)
    structure.solve()
    tip_dofs = structure.nodes[2].dofs
    assert structure.disp[tip_dofs[1]] != 0.0


# ===========================================================================
# model_to_core
# ===========================================================================


def test_model_to_core_returns_validated_structural_model():
    frame = SimpleFrame()
    frame.add_node(1, 0.0, 0.0)
    frame.add_node(2, 120.0, 0.0)
    frame.add_frame(1, 1, 2, 30000.0, 10.0, 200.0)
    frame.add_support(1, [True, True, True])
    frame.add_node_load(2, [10.0, 0.0, 0.0])

    sc_model = model_to_core(frame)

    assert isinstance(sc_model, StructuralModel)
    assert len(sc_model.nodes) == 2
    assert len(sc_model.elements) == 1
    assert len(sc_model.materials) == 1
    assert len(sc_model.sections) == 1
    assert len(sc_model.supports) == 1
    assert len(sc_model.load_cases) == 1

    validate_structural_model(sc_model, raise_error=True)

    elem = sc_model.elements[0]
    assert any(m.id == elem.material_id for m in sc_model.materials)
    assert any(s.id == elem.section_id for s in sc_model.sections)


def test_model_to_core_accepts_structure_and_results():
    structure = _build_solved_beam()
    from_struct = model_to_core(structure)
    from_results = model_to_core(Results(structure))
    assert isinstance(from_struct, StructuralModel)
    assert isinstance(from_results, StructuralModel)
    assert len(from_struct.nodes) == len(from_results.nodes) == 2


def test_model_to_core_round_trip_through_struct_core():
    frame = SimpleFrame()
    frame.add_node(1, 0.0, 0.0)
    frame.add_node(2, 100.0, 0.0)
    frame.add_node(3, 200.0, 0.0)
    frame.add_frame(1, 1, 2, 29000.0, 8.5, 150.0)
    frame.add_frame(2, 2, 3, 29000.0, 8.5, 150.0)
    frame.add_support(1, [True, True, True])
    frame.add_support(3, [True, True, True])
    frame.add_node_load(2, [0.0, -1000.0, 0.0])
    frame.solve()

    sc_model = model_to_core(frame)
    structure2 = model_from_core(sc_model)
    structure2.solve()

    assert set(structure2.nodes.keys()) == set(frame.structure.nodes.keys())
    assert len(structure2.elements) == len(frame.structure.elements)
    assert any(structure2.nodes[1].support)
    assert any(structure2.nodes[3].support)

    np.testing.assert_allclose(
        structure2.disp, frame.structure.disp, atol=1e-9, rtol=1e-9
    )


def test_model_to_core_collects_node_loads_into_load_case():
    frame = SimpleFrame()
    frame.add_node(1, 0.0, 0.0)
    frame.add_node(2, 100.0, 0.0)
    frame.add_node(3, 200.0, 0.0)
    frame.add_frame(1, 1, 2, 29000.0, 8.5, 150.0)
    frame.add_frame(2, 2, 3, 29000.0, 8.5, 150.0)
    frame.add_support(1, [True, True, True])
    frame.add_support(3, [True, True, True])
    frame.add_node_load(2, [0.0, -1000.0, 500.0])

    sc_model = model_to_core(frame)
    assert len(sc_model.load_cases) == 1
    lc = sc_model.load_cases[0]
    assert len(lc.point_loads) == 1
    pl = lc.point_loads[0]
    assert pl.node_id == 2
    np.testing.assert_allclose(pl.fy, -1000.0)
    np.testing.assert_allclose(pl.mz, 500.0)


def test_model_to_core_rejects_non_fem2d_object():
    with pytest.raises(TypeError):
        model_to_core({"not": "a structure"})
def test_model_to_core_rejects_non_fem2d_object():
    with pytest.raises(TypeError):
        model_to_core({"not": "a structure"})


# ===========================================================================
# result_to_core
# ===========================================================================


def test_result_to_core_returns_analysis_result(cantilever_results):
    results, _ = cantilever_results
    ar = result_to_core(results, analysis_case_id="StaticSolve")
    assert isinstance(ar, AnalysisResult)
    assert ar.analysis_case_id == "StaticSolve"


def test_result_to_core_node_displacements_match_results(cantilever_results):
    results, _ = cantilever_results
    df = results.node_displacements().set_index("node")
    ar = result_to_core(results)

    for nid, row in df.iterrows():
        nr = ar.node_results[nid]
        np.testing.assert_allclose(nr.displacement.ux, row["ux"], atol=1e-12)
        np.testing.assert_allclose(nr.displacement.uy, row["uy"], atol=1e-12)
        np.testing.assert_allclose(nr.displacement.rz, row["theta"], atol=1e-12)


def test_result_to_core_analytical_cantilever_deflection(cantilever_results):
    results, params = cantilever_results
    analytical_delta = (params["P"] * params["L"] ** 3) / (
        3.0 * params["E"] * params["I"]
    )
    ar = result_to_core(results)
    tip = ar.node_results[2]
    assert math.isclose(tip.displacement.uy, analytical_delta, rel_tol=1e-4)


def test_result_to_core_reactions_are_populated(cantilever_results):
    results, _ = cantilever_results
    ar = result_to_core(results)
    assert ar.node_results[1].reaction is not None
    np.testing.assert_allclose(ar.node_results[1].reaction.fy, 10000.0, rtol=1e-6)


def test_result_to_core_reactions_can_be_suppressed(cantilever_results):
    results, _ = cantilever_results
    ar = result_to_core(results, include_reactions=False)
    for nr in ar.node_results.values():
        assert nr.reaction is None


def test_result_to_core_accepts_structure_and_simple_frame():
    structure = _build_solved_beam()
    ar_struct = result_to_core(structure)
    ar_results = result_to_core(Results(structure))
    assert isinstance(ar_struct, AnalysisResult)
    assert isinstance(ar_results, AnalysisResult)

    frame = SimpleFrame()
    frame.add_node(1, 0.0, 0.0)
    frame.add_node(2, 2.0, 0.0)
    frame.add_frame(1, 1, 2, 30000.0, 10.0, 1e-4)
    frame.add_support(1, [True, True, True])
    frame.add_node_load(2, [0.0, -5000.0, 0.0])
    frame.solve()

    ar_frame = result_to_core(frame)
    assert isinstance(ar_frame, AnalysisResult)
    assert 2 in ar_frame.node_results


def test_result_to_core_element_forces_and_envelopes(cantilever_results):
    results, _ = cantilever_results
    ar = result_to_core(results)

    assert 1 in ar.element_results
    er = ar.element_results[1]
    assert len(er.forces) == 2
    assert math.isclose(er.forces[0].station, 0.0)
    assert math.isclose(er.forces[1].station, 1.0)

    df = results.element_forces().set_index("element").loc[1]
    np.testing.assert_allclose(er.forces[0].P, df["fx_i"], atol=1e-12)
    np.testing.assert_allclose(er.forces[0].Vy, df["fy_i"], atol=1e-12)
    np.testing.assert_allclose(er.forces[0].Mz, df["m_i"], atol=1e-12)
    np.testing.assert_allclose(er.forces[1].P, df["fx_j"], atol=1e-12)
    np.testing.assert_allclose(er.forces[1].Vy, df["fy_j"], atol=1e-12)
    np.testing.assert_allclose(er.forces[1].Mz, df["m_j"], atol=1e-12)

    stations_mag_P = [abs(er.forces[0].P), abs(er.forces[1].P)]
    expected_max_axial = er.forces[stations_mag_P.index(max(stations_mag_P))].P
    np.testing.assert_allclose(er.max_axial, expected_max_axial, atol=1e-12)


def test_result_to_core_has_no_units_field():
    structure = _build_solved_beam()
    ar = result_to_core(structure)
    ar.model_dump_json()
    assert "units" not in ar.model_dump()


def test_result_to_core_rejects_non_fem2d_object():
    with pytest.raises(TypeError):
        result_to_core({"not": "a structure"})