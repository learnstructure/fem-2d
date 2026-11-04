"""
Tests for support springs and shear area requirements in fem2d.
"""

import pytest
import numpy as np

from fem2d import Structure, Node, BeamElement, Section, ElasticMaterial
from struct_core import (
    StructuralModel, Node as ScNode, BeamElement as ScBeamElement,
    RectangularSection, GeneralSection, ConcreteMaterial, Support as ScSupport,
    PointLoad as ScPointLoad, LoadCase as ScLoadCase, Units
)
from fem2d.adapters.struct_core_adapter import model_from_core


def test_fem2d_cantilever_tip_spring():
    """
    Cantilever: L = 4.0 m, EI = 20,000 kN*m^2.
    Tip spring: ky = 1,000 kN/m.
    Tip load: Fy = -10 kN.
    Cantilever tip flexibility: f_cant = L^3 / (3*EI) = 64 / 60,000 = 0.00106667 m/kN.
    Cantilever tip stiffness: k_cant = 1 / f_cant = 937.5 kN/m.
    Total stiffness with spring: K_tot = k_cant + k_spring = 937.5 + 1000 = 1937.5 kN/m.
    Tip deflection: delta = 10 / 1937.5 = 0.00516129 m.
    Spring reaction: R_spring = 1000 * 0.00516129 = 5.16129 kN (upward).
    Base reaction: R_base = 937.5 * 0.00516129 = 4.83871 kN (upward).
    Total reaction = 10.0 kN.
    """
    struct = Structure()
    n1 = Node(1, 0.0, 0.0)
    n2 = Node(2, 4.0, 0.0)
    struct.add_node(n1)
    struct.add_node(n2)
    n1.set_support(ux_fixed=True, uy_fixed=True, rz_fixed=True)
    n2.set_springs(ky=1000.0)

    # EI = 20,000 -> E = 2e7, I = 0.001, A = 0.1
    mat = ElasticMaterial(E=2.0e7, nu=0.2)
    sec = Section(area=0.1, moi=0.001, material=mat)
    struct.add_element(BeamElement(1, n1, n2, sec))

    n2.set_load(fy=-10.0)
    struct.solve()

    disp_y = struct.disp[n2.dofs[1]]
    assert abs(disp_y) == pytest.approx(0.00516129, rel=1e-4)

    # Check reactions
    rxn_base_y = struct.reactions[n1.dofs[1]]
    rxn_tip_y = struct.reactions[n2.dofs[1]]
    assert rxn_tip_y == pytest.approx(5.16129, rel=1e-4)
    assert rxn_base_y == pytest.approx(4.83871, rel=1e-4)
    assert rxn_base_y + rxn_tip_y == pytest.approx(10.0, rel=1e-5)


def test_fem2d_adapter_support_springs():
    """Verify that struct_core.Support with spring stiffnesses transfers to fem2d correctly."""
    # EI: E = 2.0e7 kPa, I3 for GeneralSection = 0.001 m^4 -> EI = 20,000 kN*m^2
    model = StructuralModel(
        dimension="2D",
        units=Units(force="kN", length="m"),
        nodes=[
            ScNode(id=1, x=0.0, y=0.0),
            ScNode(id=2, x=4.0, y=0.0),
        ],
        materials=[
            ConcreteMaterial(id="C1", fc=30000.0, E=2.0e7, unit_weight=24.0),
        ],
        sections=[
            GeneralSection(id="S1", A=0.1, I3=0.001, material_id="C1"),
        ],
        elements=[
            ScBeamElement(id="B1", start_node=1, end_node=2, section_id="S1"),
        ],
        supports=[
            ScSupport(node_id=1, ux=True, uy=True, rz=True),
            ScSupport(node_id=2, ky=1000.0),
        ],
        load_cases=[
            ScLoadCase(
                id="L1",
                point_loads=[ScPointLoad(node_id=2, fy=-10.0)],
            )
        ]
    )

    struct = model_from_core(model, load_case_id="L1")
    struct.solve()

    disp_y = struct.disp[struct.nodes[2].dofs[1]]
    assert abs(disp_y) == pytest.approx(0.00516129, rel=1e-4)
    assert struct.reactions[struct.nodes[2].dofs[1]] == pytest.approx(5.16129, rel=1e-4)
    assert struct.reactions[struct.nodes[1].dofs[1]] == pytest.approx(4.83871, rel=1e-4)


def test_fem2d_support_cannot_be_both_fixed_and_spring():
    """Verify that a DOF cannot be both fixed and have a spring."""
    struct = Structure()
    n1 = Node(1, 0.0, 0.0)
    struct.add_node(n1)
    n1.set_support(uy_fixed=True)
    n1.set_springs(ky=1000.0)

    struct.number_dofs()
    with pytest.raises(ValueError, match="cannot be both restrained"):
        struct.apply_boundary_conditions()


def test_fem2d_general_section_missing_shear_area_raises_when_shear_deformation_requested():
    """Verify F12: GeneralSection without As2 raises ValueError when shear deformation is requested."""
    model = StructuralModel(
        dimension="2D",
        units=Units(force="kN", length="m"),
        nodes=[ScNode(id=1, x=0.0, y=0.0), ScNode(id=2, x=4.0, y=0.0)],
        materials=[ConcreteMaterial(id="C1", fc=30000.0, E=2.0e7, unit_weight=24.0)],
        sections=[GeneralSection(id="G1", A=0.1, I3=0.001, As2=None, material_id="C1")],
        elements=[ScBeamElement(id="B1", start_node=1, end_node=2, section_id="G1", include_shear_deformation=True)],
        supports=[ScSupport(node_id=1, ux=True, uy=True, rz=True)],
    )
    with pytest.raises(ValueError, match="without shear area 'As2'"):
        model_from_core(model)
