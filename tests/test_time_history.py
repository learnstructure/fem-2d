import numpy as np
import pytest
from fem2d import Structure, Node, ElasticMaterial, TrussElementNL, BeamElementMaterialNL
from fem2d.sections import Section, FiberSection
from fem2d.materials import BilinearMaterial
from structdyn import MDF, GroundMotion, LoadHistory


def test_fem2d_structdyn_linear_bridge():
    # Create 2D portal frame with lumped mass
    node1 = Node(1, 0.0, 0.0)
    node2 = Node(2, 0.0, 3.0)
    node2.set_mass(mass=1000.0)
    node3 = Node(3, 4.0, 3.0)
    node3.set_mass(mass=1000.0)
    node4 = Node(4, 4.0, 0.0)

    node1.set_support(ux_fixed=True, uy_fixed=True, rz_fixed=True)
    node4.set_support(ux_fixed=True, uy_fixed=True, rz_fixed=True)

    E = 200e9
    A = 0.01
    Izz = 0.0001
    mat = ElasticMaterial(E)
    sec = Section(A, Izz)

    from fem2d.elements import BeamElement
    el1 = BeamElement(1, node1, node2, mat, sec)
    el2 = BeamElement(2, node2, node3, mat, sec)
    el3 = BeamElement(3, node3, node4, mat, sec)

    struct = Structure()
    for n in [node1, node2, node3, node4]:
        struct.add_node(n)
    for el in [el1, el2, el3]:
        struct.add_element(el)

    mdf = MDF.from_fem2d(struct, zeta=0.05)
    assert mdf.ndof == 6  # 2 free nodes * 3 DOFs
    assert mdf.M.shape == (6, 6)
    assert mdf.K.shape == (6, 6)

    # Solve response under sinusoidal force on node2 ux (DOF 0)
    t = np.linspace(0, 1.0, 101)
    p = 50e3 * np.sin(2 * np.pi * 5 * t)  # 50 kN
    load_history = LoadHistory(time_steps=t, load_values=p, dof=[0])

    res = mdf.find_response(load_history, method="newmark_beta")
    assert "u1" in res.columns
    assert len(res) == 101
    assert not np.isnan(res["u1"].values).any()


def test_fem2d_structdyn_nonlinear_bridge():
    # Non-linear 2-bar truss under dynamic ground motion / pulse load
    node1 = Node(1, 0.0, 0.0)
    node2 = Node(2, 4.0, 3.0)
    node2.set_mass(mass=500.0)
    node3 = Node(3, 8.0, 0.0)

    node1.set_support(ux_fixed=True, uy_fixed=True)
    node3.set_support(ux_fixed=True, uy_fixed=True)

    E = 200e9
    sigma_y = 250e6
    E_t = 2e9  # 1% post-yield stiffness
    A = 0.002

    mat = BilinearMaterial(E, sigma_y, E_t)
    sec = Section(A)

    el1 = TrussElementNL(1, node1, node2, mat, sec)
    el2 = TrussElementNL(2, node2, node3, mat, sec)

    struct = Structure()
    for n in [node1, node2, node3]:
        struct.add_node(n)
    struct.add_element(el1)
    struct.add_element(el2)

    mdf = MDF.from_fem2d(struct)
    assert mdf.ndof == 2  # Node 2 free DOFs: ux, uy

    # Transient impulse load pushing node2 in ux (DOF 0) — pulse that causes yield
    t = np.linspace(0, 0.5, 51)
    p = np.zeros(len(t))
    p[1:10] = 400e3  # 400 kN in ux for first 0.1s
    load_history = LoadHistory(time_steps=t, load_values=p, dof=[0])

    res = mdf.find_response(load_history, method="newmark_beta")
    assert len(res) == 51
    assert not np.isnan(res["u1"].values).any()
    # Check that residual displacement exists after unloading due to plasticity
    assert abs(res["u1"].iloc[-1]) > 0.0