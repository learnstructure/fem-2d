import pytest
from fem2d import (
    SimpleFrame,
    Structure,
    Node,
    ElasticMaterial,
    Section,
    BeamElement,
    Results,
)
import math


def test_frame_5_1():
    # Example 5.1: Logan D. L. (2017), A First Course in the Finite Element Method
    # units in kips, inches
    frame = SimpleFrame()

    frame.add_node(1, 0, 0)
    frame.add_node(2, 0, 10 * 12)
    frame.add_node(3, 10 * 12, 10 * 12)
    frame.add_node(4, 10 * 12, 0)

    E, A = 30e3, 10
    I = 200

    frame.add_frame(1, 1, 2, E, A, I)
    frame.add_frame(2, 2, 3, E, A, I / 2)
    frame.add_frame(3, 3, 4, E, A, I)

    frame.add_support(1, [1, 1, 1])
    frame.add_support(4, [1, 1, 1])

    frame.add_node_load(2, [10, 0, 0])
    frame.add_node_load(3, [0, 0, 5])

    frame.solve()

    results = Results(frame)
    disp = results.node_displacements()
    el_forces = results.element_forces()

    assert abs(disp["theta"][2] - -0.0014859999862147104) < 1e-9
    assert abs(el_forces["m_i"][2] - 226.19833955969983) < 1e-9


def test_frame_5_1_v2():
    # Example 5.1: Logan D. L. (2017), A First Course in the Finite Element Method
    # units in kips, inches

    # 1. Create nodes
    node1 = Node(1, 0, 0)
    node2 = Node(2, 0, 10 * 12)  # 10 ft = 120 in
    node3 = Node(3, 10 * 12, 10 * 12)
    node4 = Node(4, 10 * 12, 0)

    # 2. Material and section properties
    E = 30e3
    A = 10
    I_full = 200
    I_half = I_full / 2

    material = ElasticMaterial(E)  # same material for all elements
    section_full = Section(A, I_full)
    section_half = Section(A, I_half)

    # 3. Create elements
    elem1 = BeamElement(1, node1, node2, material, section_full.A, section_full.I)
    elem2 = BeamElement(2, node2, node3, material, section_half.A, section_half.I)
    elem3 = BeamElement(3, node3, node4, material, section_full.A, section_full.I)

    # 4. Create structure and add nodes/elements
    structure = Structure()
    structure.add_node(node1)
    structure.add_node(node2)
    structure.add_node(node3)
    structure.add_node(node4)
    structure.add_element(elem1)
    structure.add_element(elem2)
    structure.add_element(elem3)

    # 5. Apply supports (all DOFs fixed)
    node1.set_support(ux_fixed=True, uy_fixed=True, rz_fixed=True)
    node4.set_support(ux_fixed=True, uy_fixed=True, rz_fixed=True)

    # 6. Apply nodal loads
    node2.set_load(fx=10, fy=0, mz=0)  # horizontal point load at node 2
    node3.set_load(fx=0, fy=0, mz=5)  # moment at node 3

    # 7. Solve
    structure.solve()

    # 8. Print results
    results = Results(structure)
    disp = results.node_displacements()
    el_forces = results.element_forces()
    assert abs(disp["theta"][2] - -0.0014859999862147104) < 1e-9
    assert abs(el_forces["m_i"][2] - 226.19833955969983) < 1e-9


def test_frame_5_2():
    # Example 5.2: Logan D. L. (2017), A First Course in the Finite Element Method
    # units in kips, inches
    frame = SimpleFrame()

    # units in kips, in
    frame.add_node(1, 0, 0)
    frame.add_node(2, 30 * 12, 30 * 12)
    frame.add_node(3, 70 * 12, 30 * 12)

    E = 30 * 1000
    I = 1000
    A = 100

    frame.add_frame(1, 1, 2, E=E, A=A, I=I)
    frame.add_frame(2, 2, 3, E=E, A=A, I=I)

    frame.add_support(1, [1, 1, 1])
    frame.add_support(3, [1, 1, 1])

    frame.add_distributed_load(2, wy=-1 / 12)  # udl
    frame.solve()

    results = Results(frame)
    disp = results.node_displacements()
    el_forces = results.element_forces()

    assert abs(disp["theta"][1] - -0.0032917095717463736) < 1e-9
    assert abs(el_forces["fy_i"][1] - 17.39663896899728) < 1e-9


def test_frame_6_7():
    # Example 6.7: Kassimali A. (2022), Matrix analysis of structures
    # units in kN, m
    frame = SimpleFrame()

    frame.add_node(1, 0, 0)
    frame.add_node(2, 0, 6)
    frame.add_node(3, 0, 12)
    frame.add_node(4, 9, 6)
    frame.add_node(5, 9, 0)

    E = 30e6
    I = 4.8e-4
    A = 75e-3

    frame.add_frame(1, 1, 2, E, A, I)
    frame.add_frame(2, 2, 3, E, A, I)
    frame.add_frame(3, 3, 4, E, A, I)
    frame.add_frame(4, 2, 4, E, A, I)
    frame.add_frame(5, 4, 5, E, A, I)

    frame.add_support(1, [1, 1, 1])
    frame.add_support(5, [1, 1, 1])

    frame.add_node_load(2, [80, 0, 0])  # point load
    frame.add_node_load(3, [40, 0, 0])
    frame.add_distributed_load(3, wy=12)  # udl

    frame.solve()

    results = Results(frame)
    disp = results.node_displacements()
    el_forces = results.element_forces()

    assert abs(disp["theta"][2] - 0.01789119564872733) < 1e-9
    assert abs(el_forces["m_i"][2] - -90.05939530594357) < 1e-9


def test_truss1():
    # Problem: Logan D. L. (2017), A First Course in the Finite Element Method
    # units in kips, inches
    frame = SimpleFrame()
    ft = 12
    frame.add_node(1, 0, 0)
    frame.add_node(2, 60 * ft, 0)
    frame.add_node(3, 30 * ft, 40 * ft)
    frame.add_node(4, 30 * ft, 60 * ft)

    E = 30e3
    A = 3

    frame.add_truss(1, 1, 3, E, A)
    frame.add_truss(2, 2, 3, E, A)
    frame.add_truss(3, 3, 4, E, A)

    frame.add_support(1, [1, 1, 1])
    frame.add_support(2, [1, 1, 1])
    frame.add_support(4, [1, 1, 1])

    frame.add_node_load(3, [5, -10, 0])

    frame.solve()

    results = Results(frame.structure)
    assert abs((results.node_displacements())["uy"][2] - -0.01763668430335097) < 1e-9
    assert abs(results.element_forces()["fx_i"][0] - -2.0502645502645502) < 1e-9
