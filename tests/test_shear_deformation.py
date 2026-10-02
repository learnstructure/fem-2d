"""
Tests for Timoshenko shear deformation in 2D beam elements.

References:
- Kassimali, A. (2022). Matrix Analysis of Structures, 3rd ed.,
  Chapter 9: Special Topics, Section 9.7: Shear Deformations, pp. 546-552.
- Przemieniecki, J. S. (1968). Theory of Matrix Structural Analysis.
"""

import math
import numpy as np
import pytest

from fem2d import (
    Structure,
    Node,
    ElasticMaterial,
    Section,
    BeamElement,
    SimpleFrame,
    Results,
)
from fem2d.loads import PointLoad


def test_cantilever_rectangular_exact_kassimali_example_9_6():
    """
    Test cantilever beam against Kassimali (2022) Example 9.6:
    Rectangular section: b = 0.2 m, d = 0.5 m (d/L = 0.2)
    L = 2.5 m, P = 100 kN tip load downward
    E = 200 GPa = 200e6 kPa, G = 77.2 GPa = 77.2e6 kPa
    """
    L = 2.5
    b = 0.2
    d = 0.5
    P = 100.0
    E = 200.0e6
    G = 77.2e6

    mat = ElasticMaterial(E=E, G=G)
    sec = Section.from_rectangle(width=b, depth=d)

    # 1. Euler-Bernoulli (shear neglected)
    n1_eb = Node(1, 0.0, 0.0)
    n2_eb = Node(2, L, 0.0)
    n1_eb.set_support(True, True, True)
    elem_eb = BeamElement(1, n1_eb, n2_eb, mat, sec, include_shear_deformation=False)
    s_eb = Structure()
    s_eb.add_node(n1_eb)
    s_eb.add_node(n2_eb)
    s_eb.add_element(elem_eb)
    s_eb.add_load(PointLoad(n2_eb, fx=0.0, fy=-P))
    s_eb.solve()

    res_eb = Results(s_eb)
    disp_eb = res_eb.node_displacements()
    delta_eb_expected = (P * L**3) / (3.0 * E * sec.I)
    theta_expected = -(P * L**2) / (2.0 * E * sec.I)

    assert abs(abs(disp_eb["uy"][1]) - delta_eb_expected) < 1e-12
    assert abs(disp_eb["theta"][1] - theta_expected) < 1e-12

    # 2. Timoshenko (shear included)
    n1_s = Node(1, 0.0, 0.0)
    n2_s = Node(2, L, 0.0)
    n1_s.set_support(True, True, True)
    elem_s = BeamElement(1, n1_s, n2_s, mat, sec, include_shear_deformation=True)
    s_timo = Structure()
    s_timo.add_node(n1_s)
    s_timo.add_node(n2_s)
    s_timo.add_element(elem_s)
    s_timo.add_load(PointLoad(n2_s, fx=0.0, fy=-P))
    s_timo.solve()

    res_timo = Results(s_timo)
    disp_timo = res_timo.node_displacements()

    # Analytical values from Kassimali Eq. 9.60 & Example 9.6:
    # beta_s = 12 * E * I / (G * As * L^2)
    # Delta = (P * L^3 / (3 * E * I)) * (1 + beta_s / 4)
    beta_s = (12.0 * E * sec.I) / (G * sec.As * L**2)
    assert abs(beta_s - 0.1243523316) < 1e-6

    delta_timo_expected = delta_eb_expected * (1.0 + beta_s / 4.0)

    # Deflection should include the exact +3.11% shear contribution
    assert abs(abs(disp_timo["uy"][1]) - delta_timo_expected) < 1e-12
    # Concentrated tip load rotation is unaffected by shear deformation:
    assert abs(disp_timo["theta"][1] - theta_expected) < 1e-12
    # Shear deflection must strictly exceed Euler-Bernoulli deflection
    assert abs(disp_timo["uy"][1]) > abs(disp_eb["uy"][1])


def test_cantilever_circular_section():
    """
    Test cantilever beam with circular cross section (fs = 10/9 => As = 0.9 * A).
    """
    L = 3.0
    dia = 0.3
    P = 50.0
    E = 210.0e6
    nu = 0.25
    G = E / (2.0 * (1.0 + nu))

    mat = ElasticMaterial(E=E, nu=nu)
    sec = Section.from_circle(diameter=dia)

    assert abs(sec.fs - (10.0 / 9.0)) < 1e-10
    assert abs(sec.As - 0.9 * sec.A) < 1e-10

    n1 = Node(1, 0.0, 0.0)
    n2 = Node(2, L, 0.0)
    n1.set_support(True, True, True)

    elem = BeamElement(1, n1, n2, mat, sec, include_shear_deformation=True)
    struct = Structure()
    struct.add_node(n1)
    struct.add_node(n2)
    struct.add_element(elem)
    struct.add_load(PointLoad(n2, fx=0.0, fy=-P))
    struct.solve()

    res = Results(struct)
    disp = res.node_displacements()

    beta_s = (12.0 * E * sec.I) / (G * sec.As * L**2)
    delta_theory = (P * L**3 / (3.0 * E * sec.I)) * (1.0 + beta_s / 4.0)

    assert abs(abs(disp["uy"][1]) - delta_theory) < 1e-12


def test_simply_supported_beam_shear_deformation():
    """
    Test a simply-supported beam with a central point load P modeled using 2 elements.
    Analytical midspan deflection: Delta = P * L^3 / (48 * E * I) + P * L / (4 * G * As)
                                         = (P * L^3 / (48 * E * I)) * (1 + beta_s / 4)
    """
    L = 4.0
    b = 0.25
    d = 0.6
    P = 120.0
    E = 30.0e6
    G = 12.5e6

    mat = ElasticMaterial(E=E, G=G)
    sec = Section.from_rectangle(width=b, depth=d)

    n1 = Node(1, 0.0, 0.0)
    n2 = Node(2, L / 2.0, 0.0)
    n3 = Node(3, L, 0.0)

    # Pin at node 1, roller at node 3
    n1.set_support(True, True, False)
    n3.set_support(False, True, False)

    e1 = BeamElement(1, n1, n2, mat, sec, include_shear_deformation=True)
    e2 = BeamElement(2, n2, n3, mat, sec, include_shear_deformation=True)

    struct = Structure()
    struct.add_node(n1)
    struct.add_node(n2)
    struct.add_node(n3)
    struct.add_element(e1)
    struct.add_element(e2)
    struct.add_load(PointLoad(n2, fx=0.0, fy=-P))
    struct.solve()

    res = Results(struct)
    disp = res.node_displacements()

    # Midspan deflection:
    delta_mid_fem = abs(disp["uy"][1])

    delta_M_theory = (P * L**3) / (48.0 * E * sec.I)
    delta_shear_theory = (P * L) / (4.0 * G * sec.As)
    delta_mid_theory = delta_M_theory + delta_shear_theory

    assert abs(delta_mid_fem - delta_mid_theory) < 1e-12


def test_simple_frame_shear_deformation():
    """
    Test that SimpleFrame helper works seamlessly with include_shear_deformation=True.
    """
    L = 2.5
    b = 0.2
    d = 0.5
    P = 100.0
    E = 200.0e6
    G = 77.2e6

    sec = Section.from_rectangle(width=b, depth=d)

    frame = SimpleFrame()
    frame.add_node(1, 0.0, 0.0)
    frame.add_node(2, L, 0.0)
    frame.add_support(1, [1, 1, 1])
    frame.add_frame(
        1,
        1,
        2,
        E=E,
        A=sec.A,
        I=sec.I,
        As=sec.As,
        G=G,
        include_shear_deformation=True,
    )
    frame.add_node_load(2, [0.0, -P, 0.0])
    frame.solve()

    res = Results(frame)
    disp = res.node_displacements()

    beta_s = (12.0 * E * sec.I) / (G * sec.As * L**2)
    delta_theory = (P * L**3 / (3.0 * E * sec.I)) * (1.0 + beta_s / 4.0)

    assert abs(abs(disp["uy"][1]) - delta_theory) < 1e-12


def test_shear_deformation_with_end_hinge():
    """
    Test combining Timoshenko shear deformation with static condensation for end releases.
    """
    L = 3.0
    b = 0.2
    d = 0.4
    E = 200.0e6
    nu = 0.3
    G = E / (2.0 * (1.0 + nu))

    mat = ElasticMaterial(E=E, nu=nu)
    sec = Section.from_rectangle(width=b, depth=d)

    n1 = Node(1, 0.0, 0.0)
    n2 = Node(2, L, 0.0)

    elem = BeamElement(
        1,
        n1,
        n2,
        mat,
        sec,
        hinge_j=True,
        include_shear_deformation=True,
    )

    k_local = elem.local_stiffness()

    # Released DOF 5 (theta_j) should have 0 stiffness in the condensed matrix
    assert abs(k_local[5, 5]) < 1e-14
    assert np.all(np.abs(k_local[5, :]) < 1e-14)
    assert np.all(np.abs(k_local[:, 5]) < 1e-14)

    # Matrix must be symmetric
    assert np.allclose(k_local, k_local.T)


def test_shear_deformation_disabled_by_default():
    """
    Verify that include_shear_deformation defaults to False and produces
    classical Euler-Bernoulli stiffness matrix identically.
    """
    n1 = Node(1, 0.0, 0.0)
    n2 = Node(2, 4.0, 0.0)
    mat = ElasticMaterial(E=200e6)
    sec = Section.from_rectangle(width=0.2, depth=0.4)

    elem_default = BeamElement(1, n1, n2, mat, sec)
    assert elem_default.include_shear_deformation is False

    elem_eb = BeamElement(2, n1, n2, mat, sec, include_shear_deformation=False)
    elem_timo = BeamElement(3, n1, n2, mat, sec, include_shear_deformation=True)

    assert np.allclose(elem_default.local_stiffness(), elem_eb.local_stiffness())
    # Timoshenko stiffness must differ from Euler-Bernoulli
    assert not np.allclose(elem_default.local_stiffness(), elem_timo.local_stiffness())
