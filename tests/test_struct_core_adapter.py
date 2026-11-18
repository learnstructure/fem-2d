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
        sections=[RectangularSection(id="BeamSection", material_id="Steel", b=b_val, h=h_val)],
        elements=[
            BeamElement(
                id=1,
                start_node=1,
                end_node=2,
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
        sections=[RectangularSection(id="BeamSection", material_id="Steel", b=b_val, h=h_val)],
        elements=[
            BeamElement(
                id=1,
                start_node=1,
                end_node=2,
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
        sections=[RectangularSection(id="S", material_id="M", b=0.1, h=0.2)],
        elements=[
            BeamElement(id=1, start_node=1, end_node=2, section_id="S")
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

    assert isinstance(sc_model, (StructuralModel, Project))
    assert len(sc_model.nodes) == 2
    assert len(sc_model.elements) == 1
    assert len(sc_model.materials) == 1
    assert len(sc_model.sections) == 1
    assert len(sc_model.supports) == 1
    assert len(sc_model.load_cases) == 1

    validate_structural_model(sc_model.model if isinstance(sc_model, Project) else sc_model, raise_error=True)

    elem = sc_model.elements[0]
    assert sc_model.material_of(elem) in sc_model.materials
    assert any(s.id == elem.section_id for s in sc_model.sections)


def test_model_to_core_accepts_structure_and_results():
    structure = _build_solved_beam()
    from_struct = model_to_core(structure)
    from_results = model_to_core(Results(structure))
    assert isinstance(from_struct, (StructuralModel, Project))
    assert isinstance(from_results, (StructuralModel, Project))
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
    from struct_core.analysis_result import STANDARD_STATIONS

    assert [f.station for f in er.forces] == list(STANDARD_STATIONS)

    # Internal-force convention: the i-end action is the negated node force, the j-end equals the node force
    df = results.element_forces().set_index("element").loc[1]
    np.testing.assert_allclose(er.forces[0].P, -df["fx_i"], atol=1e-9)
    np.testing.assert_allclose(er.forces[0].Vy, -df["fy_i"], atol=1e-9)
    np.testing.assert_allclose(er.forces[0].Mz, -df["m_i"], atol=1e-9)
    np.testing.assert_allclose(er.forces[-1].P, df["fx_j"], atol=1e-9)
    np.testing.assert_allclose(er.forces[-1].Vy, df["fy_j"], atol=1e-9)
    np.testing.assert_allclose(er.forces[-1].Mz, df["m_j"], atol=1e-9)

    np.testing.assert_allclose(er.max_axial, max((f.P for f in er.forces), key=abs), atol=1e-12)
    np.testing.assert_allclose(er.max_moment, max((f.Mz for f in er.forces), key=abs), atol=1e-12)


def test_result_to_core_has_no_units_field():
    structure = _build_solved_beam()
    ar = result_to_core(structure)
    ar.model_dump_json()
    assert "units" not in ar.model_dump()


def test_result_to_core_rejects_non_fem2d_object():
    with pytest.raises(TypeError):
        result_to_core({"not": "a structure"})


def test_load_combination_and_envelope_solving():
    from struct_core import LoadCombination, LoadCaseFactor, DistributedLoad
    from fem2d.adapters import solve_combinations_and_envelope

    nodes = [
        Node(id=1, x=0.0, y=0.0),
        Node(id=2, x=0.0, y=3.0),
        Node(id=3, x=5.0, y=3.0),
        Node(id=4, x=5.0, y=0.0),
    ]
    mat = ElasticMaterial(id="Conc", E=2.5e7, rho=2400.0)
    sec_beam = RectangularSection(id="B1", material_id="Conc", b=0.30, h=0.50)
    sec_col = RectangularSection(id="C1", material_id="Conc", b=0.40, h=0.40)
    elements = [
        BeamElement(id=1, start_node=1, end_node=2, section_id="C1"),
        BeamElement(id=2, start_node=2, end_node=3, section_id="B1"),
        BeamElement(id=3, start_node=4, end_node=3, section_id="C1"),
    ]
    supports = [Support.fixed(node_id=1), Support.fixed(node_id=4)]

    # Dead load: 10 kN/m on beam 2
    lc_dead = LoadCase(id="Dead", element_loads=[DistributedLoad(element_id=2, wy=-10.0)])
    # Live load: 15 kN/m on beam 2
    lc_live = LoadCase(id="Live", element_loads=[DistributedLoad(element_id=2, wy=-15.0)])

    comb1 = LoadCombination(
        id="1.4D",
        factors=[LoadCaseFactor(load_case_id="Dead", factor=1.4)],
    )
    comb2 = LoadCombination(
        id="1.2D+1.6L",
        factors=[
            LoadCaseFactor(load_case_id="Dead", factor=1.2),
            LoadCaseFactor(load_case_id="Live", factor=1.6),
        ],
    )

    model = StructuralModel(
        nodes=nodes,
        materials=[mat],
        sections=[sec_beam, sec_col],
        elements=elements,
        supports=supports,
        load_cases=[lc_dead, lc_live],
        load_combinations=[comb1, comb2],
    )
    project = Project(model=model)

    comb_results, env_result = solve_combinations_and_envelope(project)

    assert "1.4D" in comb_results
    assert "1.2D+1.6L" in comb_results
    assert env_result.analysis_case_id == "Envelope"

    # In 1.4D: w = 1.4 * 10 = 14 kN/m
    # In 1.2D+1.6L: w = 1.2 * 10 + 1.6 * 15 = 12 + 24 = 36 kN/m
    # Envelope moment on beam 2 should be governed by 1.2D+1.6L
    b2_comb1 = comb_results["1.4D"].element_results[2]
    b2_comb2 = comb_results["1.2D+1.6L"].element_results[2]
    b2_env = env_result.element_results[2]

    assert abs(b2_comb2.max_moment) > abs(b2_comb1.max_moment)
    assert abs(b2_env.max_moment) == pytest.approx(abs(b2_comb2.max_moment), rel=1e-5)


def test_model_from_core_with_self_weight_2d():
    from fem2d.adapters import model_from_core
    from fem2d.results import Results
    from struct_core import (
        StructuralModel,
        Node,
        Support,
        ElasticMaterial,
        RectangularSection,
        BeamElement,
        LoadCase,
        LoadCombination,
        LoadCaseFactor,
    )

    nodes = [
        Node(id=1, x=0.0, y=0.0),
        Node(id=2, x=10.0, y=0.0),
    ]
    mat = ElasticMaterial(id="Conc", E=3e7, rho=2400.0)
    sec = RectangularSection(id="B1", material_id="Conc", b=0.5, h=1.0)
    elements = [
        BeamElement(id=1, start_node=1, end_node=2, section_id="B1"),
    ]
    supports = [Support.fixed(node_id=1)]
    lc_dead = LoadCase(id="DEAD", name="Dead Load", include_self_weight=True)
    comb_14d = LoadCombination(
        id="1.4D",
        factors=[LoadCaseFactor(load_case_id="DEAD", factor=1.4)],
    )

    model = StructuralModel(
        nodes=nodes,
        materials=[mat],
        sections=[sec],
        elements=elements,
        supports=supports,
        load_cases=[lc_dead],
        load_combinations=[comb_14d],
    )

    # 1. Automatic self-weight included
    struct = model_from_core(model, load_case_id="DEAD", include_self_weight=True)
    assert len(struct.loads) == 1
    struct.solve()
    reac_df = Results(struct).reactions()
    node1_fy = reac_df.loc[reac_df["node"] == 1, "Fy"].values[0]
    expected_w = (2400.0 * 0.5 * 9.80665 / 1000.0) * 10.0
    assert pytest.approx(node1_fy, rel=1e-3) == expected_w

    # 2. Self-weight disabled
    struct_no_sw = model_from_core(model, load_case_id="DEAD", include_self_weight=False)
    assert len(struct_no_sw.loads) == 0

    # 3. In load combination with factor 1.4
    struct_comb = model_from_core(model, load_combination_id="1.4D", include_self_weight=True)
    assert len(struct_comb.loads) == 1
    struct_comb.solve()
    comb_reac_df = Results(struct_comb).reactions()
    comb_node1_fy = comb_reac_df.loc[comb_reac_df["node"] == 1, "Fy"].values[0]
    assert pytest.approx(comb_node1_fy, rel=1e-3) == expected_w * 1.4


def _vertical_column_model(load: "DistributedLoad") -> StructuralModel:
    return StructuralModel(
        nodes=[Node(id=1, x=0.0, y=0.0), Node(id=2, x=0.0, y=4.0)],
        materials=[ElasticMaterial(id="Conc", E=3e7, rho=2400.0)],
        sections=[RectangularSection(id="C1", material_id="Conc", b=0.4, h=0.4)],
        elements=[BeamElement(id=1, start_node=1, end_node=2, section_id="C1")],
        supports=[Support.fixed(node_id=1)],
        load_cases=[LoadCase(id="DEAD", element_loads=[load])],
    )


def test_model_from_core_global_distributed_load_on_vertical_member_is_axial():
    """A global -Y line load on a column must act axially (no lateral reaction, no base moment)."""
    from struct_core import DistributedLoad

    model = _vertical_column_model(DistributedLoad(element_id=1, wy=-5.0))
    struct = model_from_core(model, load_case_id="DEAD")
    struct.solve()
    reac = Results(struct).reactions().set_index("node").loc[1]
    assert reac["Fy"] == pytest.approx(5.0 * 4.0, rel=1e-9)
    assert reac["Fx"] == pytest.approx(0.0, abs=1e-9)
    assert reac["M"] == pytest.approx(0.0, abs=1e-9)


def test_model_from_core_local_distributed_load_on_vertical_member_is_transverse():
    """A local transverse (w2) load on a column produces a lateral reaction of w*L."""
    from struct_core import DistributedLoad

    model = _vertical_column_model(DistributedLoad(element_id=1, w2=-5.0, is_local=True))
    struct = model_from_core(model, load_case_id="DEAD")
    struct.solve()
    reac = Results(struct).reactions().set_index("node").loc[1]
    assert abs(reac["Fx"]) == pytest.approx(5.0 * 4.0, rel=1e-9)
    assert reac["Fy"] == pytest.approx(0.0, abs=1e-9)


def test_model_to_core_exports_distributed_loads_in_local_axes():
    """fem2d loads are local, so model_to_core must export them with is_local=True and survive a round trip."""
    from struct_core import DistributedLoad

    struct = model_from_core(_vertical_column_model(DistributedLoad(element_id=1, wy=-5.0)), load_case_id="DEAD")
    core = model_to_core(struct)
    dl = core.load_cases[0].element_loads[0]
    assert dl.is_local
    assert dl.w1 == pytest.approx(-5.0)
    assert dl.w2 == pytest.approx(0.0, abs=1e-12)

    struct_rt = model_from_core(core, load_case_id=core.load_cases[0].id)
    struct_rt.solve()
    reac = Results(struct_rt).reactions().set_index("node").loc[1]
    assert reac["Fy"] == pytest.approx(20.0, rel=1e-9)
    assert reac["Fx"] == pytest.approx(0.0, abs=1e-9)


def test_model_from_core_shear_deformation():
    from struct_core import (
        StructuralModel,
        Node,
        Support,
        ElasticMaterial,
        GeneralSection,
        BeamElement,
        PointLoad,
        LoadCase,
    )

    nodes = [
        Node(id=1, x=0.0, y=0.0),
        Node(id=2, x=5.0, y=0.0),
    ]
    mat = ElasticMaterial(id="Steel", E=200e6, nu=0.3)
    sec = GeneralSection(id="Sec1", material_id="Steel", A=0.01, I3=1e-4, As2=0.008)
    elements = [
        BeamElement(id=1, start_node=1, end_node=2, section_id="Sec1"),
    ]
    supports = [Support.fixed(node_id=1)]
    lc = LoadCase(id="LC1", point_loads=[PointLoad(node_id=2, fy=-100.0)])
    model = StructuralModel(nodes=nodes, materials=[mat], sections=[sec], elements=elements, supports=supports, load_cases=[lc])

    # 1. With shear deformation explicitly disabled
    struct_eb = model_from_core(model, include_shear_deformation=False)
    assert struct_eb.elements[1].include_shear_deformation is False
    struct_eb.solve()
    disp_eb = Results(struct_eb).node_displacements()
    tip_dy_eb = abs(disp_eb.loc[disp_eb["node"] == 2, "uy"].values[0])

    # 2. With shear deformation explicitly enabled
    struct_timo = model_from_core(model, include_shear_deformation=True)
    assert struct_timo.elements[1].include_shear_deformation is True
    assert struct_timo.elements[1].section.As == 0.008
    struct_timo.solve()
    disp_timo = Results(struct_timo).node_displacements()
    tip_dy_timo = abs(disp_timo.loc[disp_timo["node"] == 2, "uy"].values[0])

    # Timoshenko beam includes additional shear flexibility: tip deflection must be strictly greater
    assert tip_dy_timo > tip_dy_eb


def test_model_from_core_mass_source():
    from struct_core import (
        BeamElement,
        ElasticMaterial,
        GeneralSection,
        LoadCase,
        MassSource,
        MassSourceFactor,
        Node,
        PointLoad,
        StructuralModel,
        Support,
    )

    nodes = [
        Node(id=1, x=0.0, y=0.0),
        Node(id=2, x=0.0, y=4.0),
    ]
    mat = ElasticMaterial(id="Steel", E=200e6, rho=7850.0)
    sec = GeneralSection(id="ColSec", material_id="Steel", A=0.02, I3=2e-4)
    elements = [
        BeamElement(id=1, start_node=1, end_node=2, section_id="ColSec"),
    ]
    supports = [Support.fixed(node_id=1)]

    # Dead load: 100 kN tip vertical load
    # Live load: 50 kN tip vertical load
    lc_dead = LoadCase(id="DEAD", point_loads=[PointLoad(node_id=2, fy=-100.0)])
    lc_live = LoadCase(id="LIVE", point_loads=[PointLoad(node_id=2, fy=-50.0)])

    ms = MassSource(
        id="MS1",
        include_element_mass=True,
        include_added_mass=True,
        lateral_only=False,
        factors=[
            MassSourceFactor(load_case_id="DEAD", factor=1.0),
            MassSourceFactor(load_case_id="LIVE", factor=0.25),
        ],
    )

    model = StructuralModel(
        nodes=nodes,
        materials=[mat],
        sections=[sec],
        elements=elements,
        supports=supports,
        load_cases=[lc_dead, lc_live],
        mass_source=ms,
    )

    struct = model_from_core(model, mass_source=ms)
    assert hasattr(struct, "mass_source_result")
    # Total tip load = 100 + 0.25 * 50 = 112.5 kN
    # Lumped mass at node 2 should be approx 112.5 / 9.80665 = 11.4717 tonne
    m_node2 = struct.nodes[2].mass
    assert pytest.approx(m_node2, rel=1e-3) == 112.5 / 9.80665

    # Run modal analysis
    modal = struct.modal_analysis(num_modes=2)
    assert len(modal["frequencies"]) == 2
    assert modal["periods"][0] > 0.0


def test_model_from_core_consistent_mass_conversion():
    """Verify materials with rho in kg/m^3 and nodes with mass in kg convert to consistent metric tonnes."""
    mat = ElasticMaterial(id="Conc", E=2.5e7, rho=2400.0)
    sec = RectangularSection(id="ColSec", material_id="Conc", b=0.4, h=0.4)
    nodes = [
        Node(id=1, x=0.0, y=0.0),
        Node(id=2, x=0.0, y=3.5, mass=5000.0),
    ]
    elements = [
        BeamElement(id=1, start_node=1, end_node=2, section_id="ColSec"),
    ]
    supports = [Support.fixed(node_id=1)]
    model = StructuralModel(nodes=nodes, materials=[mat], sections=[sec], elements=elements, supports=supports)

    struct = model_from_core(model)
    # Density should be converted from 2400.0 kg/m^3 to 2.4 tonne/m^3
    assert pytest.approx(struct.elements[1].section.material.rho, rel=1e-5) == 2.4
    # Node mass should be converted from 5000.0 kg to 5.0 tonnes
    assert pytest.approx(struct.nodes[2].mass, rel=1e-5) == 5.0

    # Modal period of ~5-tonne tip on 3.5m column (K = 12*E*I/L^3 = 14840 kN/m):
    # omega ~ sqrt(14840 / 5.0) ~ 54.5 rad/s -> T ~ 2*pi / 54.5 ~ 0.115 s (NOT 3.65 s!)
    modal = struct.modal_analysis(num_modes=1)
    T1 = modal["periods"][0]
    assert 0.05 < T1 < 0.3, f"Expected realistic period ~0.1s, got {T1}s"


def test_fem2d_model_from_core_load_cases_filtering():
    """Verify 2D load case transfer filtering and default summation."""
    from struct_core import (
        StructuralModel,
        Node,
        Support,
        ElasticMaterial,
        RectangularSection,
        BeamElement,
        PointLoad,
        LoadCase,
    )

    nodes = [
        Node(id=1, x=0.0, y=0.0),
        Node(id=2, x=0.0, y=3.0),
    ]
    mat = ElasticMaterial(id="M1", E=3e7)
    sec = RectangularSection(id="S1", material_id="M1", b=0.3, h=0.3)
    elements = [
        BeamElement(id=1, start_node=1, end_node=2, section_id="S1"),
    ]
    supports = [Support.fixed(node_id=1)]

    lc_dead = LoadCase(id="DEAD", name="Dead Load", point_loads=[PointLoad(node_id=2, fy=-100.0)])
    lc_live = LoadCase(id="LIVE", name="Live Load", point_loads=[PointLoad(node_id=2, fy=-50.0)])

    model = StructuralModel(
        nodes=nodes,
        materials=[mat],
        sections=[sec],
        elements=elements,
        supports=supports,
        load_cases=[lc_dead, lc_live],
    )

    # 1. Default: sum all load cases (-100 + -50 = -150)
    struct_all = model_from_core(model)
    struct_all.solve()
    assert pytest.approx(struct_all.reactions[struct_all.nodes[1].dofs[1]], rel=1e-5) == 150.0

    # 2. Explicit 'DEAD' only
    struct_dead = model_from_core(model, load_case_id="DEAD")
    struct_dead.solve()
    assert pytest.approx(struct_dead.reactions[struct_dead.nodes[1].dofs[1]], rel=1e-5) == 100.0

    # 3. Explicit 'LIVE' only
    struct_live = model_from_core(model, load_case_id="LIVE")
    struct_live.solve()
    assert pytest.approx(struct_live.reactions[struct_live.nodes[1].dofs[1]], rel=1e-5) == 50.0

    # 4. Explicit list of cases ['DEAD', 'LIVE']
    struct_list = model_from_core(model, load_case_id=["DEAD", "LIVE"])
    struct_list.solve()
    assert pytest.approx(struct_list.reactions[struct_list.nodes[1].dofs[1]], rel=1e-5) == 150.0

    # 5. Invalid case raises ValueError
    with pytest.raises(ValueError):
        model_from_core(model, load_case_id="WIND")


def test_fem2d_ecosystem_roundtrip_units_and_unit_weight():
    from fem2d import Structure, Node, ElasticMaterial, Section, BeamElement
    from fem2d.loads import DistributedLoad
    from fem2d.adapters import model_to_core, model_from_core
    from struct_load import apply_self_weight

    struct = Structure(units="kN, m")
    n1 = Node(1, 0.0, 0.0)
    n1.set_support([True, True, True])
    n2 = Node(2, 5.0, 0.0)
    struct.add_node(n1)
    struct.add_node(n2)

    mat = ElasticMaterial(E=30e6, nu=0.2, unit_weight=24.0)
    sec = Section.from_rectangle(width=0.3, depth=0.5, material=mat)  # A = 0.15 m^2
    elem = BeamElement(1, n1, n2, section=sec)
    struct.add_element(elem)

    # Superimposed live load (-10 kN/m along global -Y)
    struct.add_load(DistributedLoad(elem, wy=-10.0))

    # Convert to struct_core
    core_model = model_to_core(struct)
    assert core_model.units.force == "kN"
    assert core_model.units.length == "m"
    assert core_model.materials[0].unit_weight == 24.0

    # Apply self weight: w = 24.0 * 0.15 = 3.6 kN/m along -y
    core_model = apply_self_weight(core_model)
    dead_case = next(lc for lc in core_model.load_cases if lc.id == "DEAD")
    assert len(dead_case.element_loads) == 1
    assert pytest.approx(dead_case.element_loads[0].wy, rel=1e-5) == -3.6

    # Convert back to fem2d Structure
    struct_rebuilt = model_from_core(core_model)
    assert struct_rebuilt.units.force == "kN"
    assert struct_rebuilt.elements[1].section.material.unit_weight == 24.0
    # Both DEAD and LIVE transferred: total w = -10.0 + -3.6 = -13.6 kN/m, total load = 5.0 * 13.6 = 68.0 kN
    struct_rebuilt.solve()
    assert pytest.approx(struct_rebuilt.reactions[struct_rebuilt.nodes[1].dofs[1]], rel=1e-5) == 68.0

def test_result_to_core_reports_internal_forces_along_members():
    """Fixed-fixed beam under UDL: |M| = wL^2/12 at the ends, wL^2/24 at midspan; V = wL/2 -> 0."""
    from struct_core import DistributedLoad
    from struct_core.analysis_result import STANDARD_STATIONS

    w, L = 10.0, 6.0
    model = StructuralModel(
        nodes=[Node(id=1, x=0.0, y=0.0), Node(id=2, x=L, y=0.0)],
        materials=[ElasticMaterial(id="Conc", E=3e7, rho=2400.0)],
        sections=[RectangularSection(id="B", material_id="Conc", b=0.3, h=0.5)],
        elements=[BeamElement(id=1, start_node=1, end_node=2, section_id="B")],
        supports=[Support.fixed(node_id=1), Support.fixed(node_id=2)],
        load_cases=[LoadCase(id="DEAD", element_loads=[DistributedLoad(element_id=1, wy=-w)])],
    )
    struct = model_from_core(model, load_case_id="DEAD")
    struct.solve()
    stations = result_to_core(struct).element_results[1].forces
    assert [s.station for s in stations] == list(STANDARD_STATIONS)
    assert abs(stations[0].M3) == pytest.approx(w * L**2 / 12.0, rel=1e-6)
    assert abs(stations[5].M3) == pytest.approx(w * L**2 / 24.0, rel=1e-6)
    assert abs(stations[-1].M3) == pytest.approx(w * L**2 / 12.0, rel=1e-6)
    assert abs(stations[0].V2) == pytest.approx(w * L / 2.0, rel=1e-6)
    assert stations[5].V2 == pytest.approx(0.0, abs=1e-6)
    # sagging and hogging moments have opposite signs
    assert stations[0].M3 * stations[5].M3 < 0.0


def test_flexural_stiffness_modifiers_scale_bending_deflection():
    """Cantilever tip deflection P L^3 / (3 E I) doubles when I is halved (ACI cracked-section modifiers)."""
    from struct_core import PointLoad

    model = StructuralModel(
        nodes=[Node(id=1, x=0.0, y=0.0), Node(id=2, x=0.0, y=3.0)],
        materials=[ElasticMaterial(id="Conc", E=3e7, rho=0.0)],
        sections=[RectangularSection(id="C", material_id="Conc", b=0.4, h=0.4)],
        elements=[BeamElement(id=1, start_node=1, end_node=2, section_id="C")],
        supports=[Support.fixed(node_id=1)],
        load_cases=[LoadCase(id="H", point_loads=[PointLoad(node_id=2, fx=10.0)])],
    )
    gross = model_from_core(model, load_case_id="H")
    gross.solve()
    cracked = model_from_core(model, load_case_id="H", flexural_stiffness_modifiers={1: 0.5})
    cracked.solve()
    u_gross = gross.disp[gross.nodes[2].dofs[0]]
    u_cracked = cracked.disp[cracked.nodes[2].dofs[0]]
    assert u_cracked == pytest.approx(2.0 * u_gross, rel=1e-9)
    # area (axial stiffness) is not modified
    assert cracked.elements[1].section.A == gross.elements[1].section.A


def test_element_stiffness_modifier_i3_scales_in_plane_bending():
    """The struct_core key I3 (in-plane bending of a 2D model) is applied; the non-canonical keys I and As that fem2d
    also accepted are refused by struct_core."""
    from struct_core import PointLoad

    def tip(mods):
        model = StructuralModel(
            nodes=[Node(id=1, x=0.0, y=0.0), Node(id=2, x=0.0, y=3.0)],
            materials=[ElasticMaterial(id="Conc", E=3e7, rho=0.0)],
            sections=[RectangularSection(id="C", material_id="Conc", b=0.4, h=0.4)],
            elements=[BeamElement(id=1, start_node=1, end_node=2, section_id="C", stiffness_modifiers=mods)],
            supports=[Support.fixed(node_id=1)],
            load_cases=[LoadCase(id="H", point_loads=[PointLoad(node_id=2, fx=10.0)])],
        )
        structure = model_from_core(model, load_case_id="H")
        structure.solve()
        return structure.disp[structure.nodes[2].dofs[0]]

    assert tip({"I3": 0.5}) == pytest.approx(2.0 * tip(None), rel=1e-9)
    assert tip({"I2": 0.5, "J": 0.1}) == pytest.approx(tip(None), rel=1e-12)  # out-of-plane keys do not apply
    with pytest.raises(ValueError, match="Unknown stiffness modifiers"):
        BeamElement(id=1, start_node=1, end_node=2, section_id="C", stiffness_modifiers={"I": 0.5})


def _cantilever_modes(stiffness_modifiers):
    """(bending, axial) frequencies of a 4-element horizontal cantilever, the axial modes found by their shape."""
    model = StructuralModel(
        dimension="2D",
        nodes=[Node(id=i, x=float(i), y=0.0) for i in range(5)],
        materials=[ElasticMaterial(id="C", E=25e6, nu=0.2, unit_weight=24.0)],
        sections=[RectangularSection(id="R", material_id="C", b=0.3, h=0.5)],
        elements=[BeamElement(id=i, start_node=i, end_node=i + 1, section_id="R",
                              stiffness_modifiers=stiffness_modifiers) for i in range(4)],
        supports=[Support.fixed(node_id=0)],
    )
    structure = model_from_core(model, units="kN, m")
    result = structure.modal_analysis(num_modes=12)
    axial = np.abs(np.asarray(result["modes"])[structure.nodes[4].dofs[0], :]) > 0.99
    frequencies = np.asarray(result["frequencies"])
    return frequencies[~axial], frequencies[axial]


def test_mass_and_area_modifiers_act_as_in_sap2000():
    """SAP2000 27.1 (measured, struct_link tests): the mass modifier scales the mass only, the A modifier the axial
    stiffness only. Regression: the mass modifier was ignored and the A modifier also scaled the element mass."""
    bending, axial = _cantilever_modes(None)
    bending_m, axial_m = _cantilever_modes({"mass": 0.5})
    np.testing.assert_allclose(bending_m, bending * math.sqrt(2.0), rtol=1e-9)
    np.testing.assert_allclose(axial_m, axial * math.sqrt(2.0), rtol=1e-9)
    bending_a, axial_a = _cantilever_modes({"A": 0.5})
    np.testing.assert_allclose(bending_a, bending, rtol=1e-9)
    np.testing.assert_allclose(axial_a, axial / math.sqrt(2.0), rtol=1e-9)
    assert len(axial) >= 1 and len(bending) >= 4


def test_round_trip_keeps_masses_densities_and_bilinear_materials():
    """model_to_core wrote fem2d's consistent masses (t, t/m3) as model masses (kg, kg/m3), and read a bilinear
    material's E_t / sigma_y under names it does not have (Et = 0.01 E, fy = nan); model_from_core passed them to the
    constructor under those names too."""
    from struct_core import BilinearMaterial as ScBilinear, Metadata, Project, TrussElement, parse_units

    model = StructuralModel(
        nodes=[Node(id=1, x=0.0, y=0.0), Node(id=2, x=3.0, y=0.0, mass=500.0)],
        supports=[Support.fixed(node_id=1)],
        materials=[ElasticMaterial(id="C", E=2.5e7, rho=2400.0),
                   ScBilinear(id="S", E=2e8, Et=2e6, fy=2.5e5, rho=7850.0)],
        sections=[RectangularSection(id="R", b=0.3, h=0.5, material_id="C"),
                  RectangularSection(id="T", b=0.1, h=0.1, material_id="S")],
        elements=[BeamElement(id=1, start_node=1, end_node=2, section_id="R"),
                  TrussElement(id=2, start_node=1, end_node=2, section_id="T")],
    )
    struct = model_from_core(Project(model=model, metadata=Metadata(units=parse_units("kN, m"))))
    assert struct.elements[2].section.material.E_t == pytest.approx(2e6)
    back = model_to_core(struct)
    mats = sorted(back.materials, key=lambda m: m.type)
    assert mats[0].type == "bilinear" and (mats[0].Et, mats[0].fy, mats[0].rho) == pytest.approx((2e6, 2.5e5, 7850.0))
    assert mats[1].rho == pytest.approx(2400.0)
    assert next(n for n in back.nodes if n.id == 2).mass == pytest.approx(500.0)
