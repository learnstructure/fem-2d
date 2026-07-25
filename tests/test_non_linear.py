from fem2d import BeamElementMaterialNL
import pytest
from fem2d import Structure, Node, ElasticMaterial, TrussElementNL, Results


def test_truss_nonlinear():
    # Example 10.2: Kassimali A. (2022), Matrix analysis of structures
    # units in kN, m

    # 1. Create nodes
    node1 = Node(1, 0.0, 0.0)
    node2 = Node(2, 4.0, 3.0)
    node3 = Node(3, 8.0, 0.0)

    # 2. Material properties
    E = 200e6  # 200 GPa = 200×10⁶ kN/m²
    EA = 45155.0
    A = EA / E

    material = ElasticMaterial(E)
    from fem2d.sections import Section
    section = Section(A)

    # 3. Create nonlinear truss elements
    elem1 = TrussElementNL(1, node1, node2, material, section)
    elem2 = TrussElementNL(2, node2, node3, material, section)
    elem3 = TrussElementNL(3, node3, node1, material, section)

    # 4. Build structure
    structure = Structure()
    structure.add_node(node1)
    structure.add_node(node2)
    structure.add_node(node3)
    structure.add_element(elem1)
    structure.add_element(elem2)
    structure.add_element(elem3)

    # 5. Apply supports – FIX ROTATIONS FOR ALL NODES!
    node1.set_support(ux_fixed=True, uy_fixed=True)  # pinned
    node3.set_support(ux_fixed=False, uy_fixed=True)  # roller

    # 6. Apply external load (vertical -2000 kN at node 2)
    node2.set_load(fx=0.0, fy=-2000.0, mz=0.0)

    # 7. Run nonlinear analysis (one line!)
    structure.solve_nonlinear(tolerance=1e-8, max_iter=30)

    results = Results(structure)
    disp = results.node_displacements()
    el_forces = results.element_forces()

    assert abs(disp["uy"][1] - -0.6499578914694384) < 1e-9
    assert abs(el_forces["fx_i"][2] - -1768.7701456115713) < 1e-9


def test_bilinear_material():
    from fem2d.materials import BilinearMaterial

    E = 200e9
    sigma_y = 250e6
    E_t = 20e9

    # pure kinematic hardening
    mat = BilinearMaterial(E, sigma_y, E_t, beta=0.0)

    # 1. Elastic strain (0.001 < eps_y = 0.00125)
    s1, t1, state1 = mat.get_response(0.001)
    assert abs(s1 - 200e6) < 1e-7
    assert abs(t1 - E) < 1e-7

    # 2. Plastic strain (0.002 > eps_y)
    # Expected stress: 250e6 + 20e9 * (0.002 - 0.00125) = 265e6
    s2, t2, state2 = mat.get_response(0.002, state=state1)
    assert abs(s2 - 265e6) < 1e-7
    assert abs(t2 - E_t) < 1e-7

    # 3. Unloading from 0.002 back to 0.0015 (should be elastic with slope E)
    # Expected stress: 265e6 - 200e9 * (0.002 - 0.0015) = 165e6
    s3, t3, state3 = mat.get_response(0.0015, state=state2)
    assert abs(s3 - 165e6) < 1e-7
    assert abs(t3 - E) < 1e-7


def test_csv_material(tmp_path):
    import pandas as pd
    from fem2d.materials import CSVMaterial

    csv_file = tmp_path / "curve.csv"
    # Symmetric positive-only quadrant data points
    df = pd.DataFrame({
        "strain": [0.0, 0.00125, 0.002, 0.003],
        "stress": [0.0, 250e6, 265e6, 285e6]
    })
    df.to_csv(csv_file, index=False)

    # Path-independent
    mat = CSVMaterial(str(csv_file), x_col="strain", y_col="stress", path_dependent=False)

    # Interpolation within bounds
    s1, t1, _ = mat.get_response(0.001)
    assert abs(s1 - 200e6) < 1e-7
    assert abs(t1 - 200e9) < 1e-7

    # Negative strain (symmetric reflection)
    s2, t2, _ = mat.get_response(-0.001)
    assert abs(s2 - -200e6) < 1e-7
    assert abs(t2 - 200e9) < 1e-7

    # Path-dependent
    mat_pd = CSVMaterial(str(csv_file), x_col="strain", y_col="stress", path_dependent=True)
    _, _, state1 = mat_pd.get_response(0.002) # yield
    # Unload: initial slope E_init = 200e9
    s_unload, t_unload, _ = mat_pd.get_response(0.0015, state=state1)
    # Expected: 265e6 - 200e9 * (0.002 - 0.0015) = 165e6
    assert abs(s_unload - 165e6) < 1e-7
    assert abs(t_unload - 200e9) < 1e-7


def test_beam_materialNL_cantilever():
    from fem2d import Structure, Node
    from fem2d.materials import BilinearMaterial
    from fem2d.sections import FiberSection

    # Create structure
    node1 = Node(1, 0.0, 0.0)
    node2 = Node(2, 2.0, 0.0)

    # Support Node 1 (fixed cantilever)
    node1.set_support(ux_fixed=True, uy_fixed=True, rz_fixed=True)

    # Material & Section
    E = 200e9
    sigma_y = 250e6
    E_t = 2e9  # 1% hardening
    material = BilinearMaterial(E, sigma_y, E_t)
    section = FiberSection.from_rectangular_fibers(width=0.1, depth=0.2, material=material, n_fibers=20)

    # Create element
    elem = BeamElementMaterialNL(1, node1, node2, section, n_int_pts=3)

    structure = Structure()
    structure.add_node(node1)
    structure.add_node(node2)
    structure.add_element(elem)

    # Linear load: -10 kN (stress is below yield)
    node2.set_load(fx=0.0, fy=-10e3, mz=0.0)
    structure.solve_nonlinear(n_steps=1)
    disp = structure.disp[node2.dofs]

    # Analytical elastic displacement: w = P * L^3 / (3 * E * I) = -0.002
    assert abs(disp[1] - (-0.002)) < 1e-4

    # Large load: -500 kN (causes plastic yielding)
    node2.set_load(fx=0.0, fy=-500e3, mz=0.0)
    # Solve with multiple steps
    structure.solve_nonlinear(n_steps=50)
    disp_plastic = structure.disp[node2.dofs]

    # Displacement should be significantly larger than the linear elastic extrapolation (-0.16 m)
    assert disp_plastic[1] < -0.25

