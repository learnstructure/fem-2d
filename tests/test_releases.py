import pytest
import numpy as np

from fem2d import Node, BeamElement, Structure, ElasticMaterial, Section
from fem2d.loads import DistributedLoad, ElementPointLoad, TriangularLoad

from fem3d import Node as Node3D, FrameElement as Frame3D, Structure as Structure3D
from fem3d.loads import DistributedLoad as DistributedLoad3D
from fem3d.materials.elastic import ElasticMaterial as ElasticMaterial3D
from fem3d.sections.section import Section as Section3D


def portal_reference(rel_i, rel_j, E=200e6, A=0.01, I=1e-4, q=-10.0, L=8.0, h=3.0, load_type="uniform", py=-50.0, a_ratio=0.3):
    """Beam M3 (sagging +) at x = 0, L/2, L for a fixed-base portal with beam end releases
    rel_i / rel_j = [u1, u2, r3] in beam local axes. Euler-Bernoulli."""
    nodes = {1: (0, 0), 2: (0, h), 3: (L, h), 4: (L, 0)}

    def kl(Le):
        a, b, c, d = E * A / Le, 12 * E * I / Le**3, 6 * E * I / Le**2, 4 * E * I / Le
        return np.array(
            [
                [a, 0, 0, -a, 0, 0],
                [0, b, c, 0, -b, c],
                [0, c, d, 0, -c, d / 2],
                [-a, 0, 0, a, 0, 0],
                [0, -b, -c, 0, b, -c],
                [0, c, d / 2, 0, -c, d],
            ]
        )

    rel = [k for k, f in enumerate(list(rel_i) + list(rel_j)) if f]
    K = np.zeros((12, 12))
    F = np.zeros(12)
    keep = {}
    for eid, i, j in (("C1", 1, 2), ("B", 2, 3), ("C2", 4, 3)):
        (xi, yi), (xj, yj) = nodes[i], nodes[j]
        Le = np.hypot(xj - xi, yj - yi)
        cx, cy = (xj - xi) / Le, (yj - yi) / Le
        k = kl(Le)
        if eid == "B":
            if load_type == "uniform":
                f = np.array([0, q * Le / 2, q * Le**2 / 12, 0, q * Le / 2, -q * Le**2 / 12])
            elif load_type == "point":
                a = a_ratio * Le
                b = Le - a
                f = np.array([
                    0,
                    py * b**2 * (3 * a + b) / Le**3,
                    py * a * b**2 / Le**2,
                    0,
                    py * a**2 * (a + 3 * b) / Le**3,
                    -py * a**2 * b / Le**2,
                ])
            elif load_type == "triangular":
                f = np.array([
                    0,
                    3 * q * Le / 20,
                    q * Le**2 / 30,
                    0,
                    7 * q * Le / 20,
                    -q * Le**2 / 20,
                ])
        else:
            f = np.zeros(6)

        if eid == "B" and rel:  # condense stiffness AND fixed-end forces
            r = [n for n in range(6) if n not in rel]
            inv = np.linalg.inv(k[np.ix_(rel, rel)])
            kc = np.zeros((6, 6))
            kc[np.ix_(r, r)] = k[np.ix_(r, r)] - k[np.ix_(r, rel)] @ inv @ k[np.ix_(rel, r)]
            fc = np.zeros(6)
            fc[r] = f[r] - k[np.ix_(r, rel)] @ inv @ f[rel]
            k, f = kc, fc
        t3 = np.array([[cx, cy, 0], [-cy, cx, 0], [0, 0, 1]])
        t = np.zeros((6, 6))
        t[:3, :3] = t3
        t[3:, 3:] = t3
        d = [3 * (i - 1), 3 * (i - 1) + 1, 3 * (i - 1) + 2, 3 * (j - 1), 3 * (j - 1) + 1, 3 * (j - 1) + 2]
        K[np.ix_(d, d)] += t.T @ k @ t
        F[d] += t.T @ f
        keep[eid] = (k, f, t, d, Le)
    free = [3, 4, 5, 6, 7, 8]
    u = np.zeros(12)
    u[free] = np.linalg.solve(K[np.ix_(free, free)], F[free])
    k, f, t, d, Le = keep["B"]
    e = k @ (t @ u[d]) - f  # local end forces of the beam

    if load_type == "uniform":
        M = lambda x: -e[2] + e[1] * x + q * x**2 / 2
    elif load_type == "point":
        a = a_ratio * Le
        M = lambda x: -e[2] + e[1] * x + (py * (x - a) if x >= a else 0.0)
    elif load_type == "triangular":
        M = lambda x: -e[2] + e[1] * x + q * x**3 / (6 * Le)

    return [M(0.0), M(Le / 2), M(Le)]


def _build_fem2d_portal(rel_i, rel_j, include_shear=False):
    mat = ElasticMaterial(E=200e6, nu=0.3)
    sec = Section(0.01, 1e-4, material=mat)
    if include_shear:
        sec.As = 0.01 * 5.0 / 6.0

    s = Structure()
    n1 = Node(1, 0.0, 0.0)
    n2 = Node(2, 0.0, 3.0)
    n3 = Node(3, 8.0, 3.0)
    n4 = Node(4, 8.0, 0.0)
    for n in (n1, n2, n3, n4):
        s.add_node(n)

    n1.set_support(True, True, True)
    n4.set_support(True, True, True)

    c1 = BeamElement("C1", n1, n2, mat, sec, include_shear_deformation=include_shear)
    beam = BeamElement(
        "B",
        n2,
        n3,
        mat,
        sec,
        releases_i=rel_i,
        releases_j=rel_j,
        include_shear_deformation=include_shear,
    )
    c2 = BeamElement("C2", n4, n3, mat, sec, include_shear_deformation=include_shear)

    s.add_element(c1)
    s.add_element(beam)
    s.add_element(c2)

    load = DistributedLoad(beam, wy=-10.0)
    s.add_load(load)
    return s, beam


def _build_fem3d_portal(rel_i, rel_j, include_shear=False):
    # Portal in XZ plane (X horizontal, Z vertical upward)
    mat = ElasticMaterial3D(E=200e6, nu=0.3)
    as_val = 0.01 * 5.0 / 6.0 if include_shear else None
    sec = Section3D(A=0.01, I2=1e-4, I3=1e-4, J=1e-4, As2=as_val, As3=as_val, material=mat)

    s = Structure3D()
    n1 = Node3D(1, 0.0, 0.0, 0.0)
    n2 = Node3D(2, 0.0, 0.0, 3.0)
    n3 = Node3D(3, 8.0, 0.0, 3.0)
    n4 = Node3D(4, 8.0, 0.0, 0.0)
    for n in (n1, n2, n3, n4):
        s.add_node(n)

    # Base supports: fully fixed
    n1.set_support([True, True, True, True, True, True])
    n4.set_support([True, True, True, True, True, True])
    # Restrain out-of-plane DOFs: Uy, Rx, Rz
    n2.set_support([False, True, False, True, False, True])
    n3.set_support([False, True, False, True, False, True])

    # Convert 3-entry 2D release [u1, u2, r3] to 6-entry 3D release [u1, u2, u3, r1, r2, r3]
    rel_i_3d = [rel_i[0], rel_i[1], False, False, False, rel_i[2]] if any(rel_i) else None
    rel_j_3d = [rel_j[0], rel_j[1], False, False, False, rel_j[2]] if any(rel_j) else None

    c1 = Frame3D("C1", n1, n2, mat, sec, include_shear_deformation=include_shear)
    beam = Frame3D(
        "B",
        n2,
        n3,
        mat,
        sec,
        releases_i=rel_i_3d,
        releases_j=rel_j_3d,
        include_shear_deformation=include_shear,
    )
    c2 = Frame3D("C2", n4, n3, mat, sec, include_shear_deformation=include_shear)

    s.add_element(c1)
    s.add_element(beam)
    s.add_element(c2)

    load = DistributedLoad3D(beam, w2=-10.0)
    s.add_load(load)
    return s, beam


def _beam_moments_fem2d(beam):
    # In fem2d: beam spans along X from n2 to n3
    # Internal bending moment M3 along member: M(x) = -f[2] + f[1]*x + q*x^2/2
    f = beam.get_local_forces()
    q = -10.0
    L = 8.0
    M = lambda x: -f[2] + f[1] * x + q * x**2 / 2.0
    return [M(0.0), M(L / 2.0), M(L)]


def _beam_moments_fem3d(beam):
    # In fem3d: f is 12-vector: f1_i, f2_i, f3_i, m1_i, m2_i, m3_i, f1_j, ...
    f = beam.get_local_forces()
    q = -10.0
    L = 8.0
    M = lambda x: -f[5] + f[1] * x + q * x**2 / 2.0
    return [M(0.0), M(L / 2.0), M(L)]


@pytest.mark.parametrize(
    "name, rel_i, rel_j",
    [
        ("None", [False, False, False], [False, False, False]),
        ("Moment hinge at i", [False, False, True], [False, False, False]),
        ("Moment hinge at j only", [False, False, False], [False, False, True]),
        ("Hinges at both ends", [False, False, True], [False, False, True]),
        ("Shear (u2) release at i", [False, True, False], [False, False, False]),
        ("Axial (u1) release at i", [True, False, False], [False, False, False]),
    ],
)
def test_portal_beam_releases_match_reference_solver(name, rel_i, rel_j):
    ref_m = portal_reference(rel_i, rel_j)
    s, beam = _build_fem2d_portal(rel_i, rel_j, include_shear=False)
    s.solve()
    fem2d_m = _beam_moments_fem2d(beam)

    np.testing.assert_allclose(fem2d_m, ref_m, atol=0.05, err_msg=f"Failed on case {name}")


def test_portal_point_load_with_releases():
    """Verify releases with ElementPointLoad at 0.3L against reference solver."""
    rel_i = [False, False, False]
    rel_j = [False, False, True]
    ref_m = portal_reference(rel_i, rel_j, load_type="point", py=-50.0, a_ratio=0.3)

    mat = ElasticMaterial(E=200e6)
    sec = Section(0.01, 1e-4, material=mat)
    s = Structure()
    n1, n2, n3, n4 = Node(1, 0, 0), Node(2, 0, 3), Node(3, 8, 3), Node(4, 8, 0)
    for n in (n1, n2, n3, n4):
        s.add_node(n)
    n1.set_support(True, True, True)
    n4.set_support(True, True, True)

    c1 = BeamElement("C1", n1, n2, mat, sec)
    beam = BeamElement("B", n2, n3, mat, sec, releases_i=rel_i, releases_j=rel_j)
    c2 = BeamElement("C2", n4, n3, mat, sec)
    for el in (c1, beam, c2):
        s.add_element(el)

    s.add_load(ElementPointLoad(beam, py=-50.0, x=2.4))
    s.solve()

    f = beam.get_local_forces()
    M = lambda x: -f[2] + f[1] * x + (-50.0 * (x - 2.4) if x >= 2.4 else 0.0)
    fem2d_m = [M(0.0), M(4.0), M(8.0)]

    np.testing.assert_allclose(fem2d_m, ref_m, atol=0.05)


def test_portal_triangular_load_with_releases():
    """Verify releases with TriangularLoad against reference solver."""
    rel_i = [False, False, True]
    rel_j = [False, False, False]
    ref_m = portal_reference(rel_i, rel_j, load_type="triangular", q=-10.0)

    mat = ElasticMaterial(E=200e6)
    sec = Section(0.01, 1e-4, material=mat)
    s = Structure()
    n1, n2, n3, n4 = Node(1, 0, 0), Node(2, 0, 3), Node(3, 8, 3), Node(4, 8, 0)
    for n in (n1, n2, n3, n4):
        s.add_node(n)
    n1.set_support(True, True, True)
    n4.set_support(True, True, True)

    c1 = BeamElement("C1", n1, n2, mat, sec)
    beam = BeamElement("B", n2, n3, mat, sec, releases_i=rel_i, releases_j=rel_j)
    c2 = BeamElement("C2", n4, n3, mat, sec)
    for el in (c1, beam, c2):
        s.add_element(el)

    s.add_load(TriangularLoad(beam, wy1=0.0, wy2=-10.0))
    s.solve()

    f = beam.get_local_forces()
    M = lambda x: -f[2] + f[1] * x + (-10.0 * x**3 / (6 * 8.0))
    fem2d_m = [M(0.0), M(4.0), M(8.0)]

    np.testing.assert_allclose(fem2d_m, ref_m, atol=0.05)


@pytest.mark.parametrize(
    "rel_i, rel_j",
    [
        ([False, False, False], [False, False, False]),
        ([False, False, True], [False, False, False]),
        ([False, False, False], [False, False, True]),
        ([False, False, True], [False, False, True]),
        ([False, True, False], [False, False, False]),
        ([True, False, False], [False, False, False]),
    ],
)
def test_portal_fem2d_fem3d_parity(rel_i, rel_j):
    s2d, b2d = _build_fem2d_portal(rel_i, rel_j, include_shear=False)
    s2d.solve()
    m2d = _beam_moments_fem2d(b2d)

    s3d, b3d = _build_fem3d_portal(rel_i, rel_j, include_shear=False)
    s3d.solve()
    m3d = _beam_moments_fem3d(b3d)

    np.testing.assert_allclose(m2d, m3d, atol=1e-3)


def test_portal_shear_deformation_fem2d_fem3d_parity():
    rel_i = [False, False, True]
    rel_j = [False, False, False]
    s2d, b2d = _build_fem2d_portal(rel_i, rel_j, include_shear=True)
    s2d.solve()
    m2d = _beam_moments_fem2d(b2d)

    s3d, b3d = _build_fem3d_portal(rel_i, rel_j, include_shear=True)
    s3d.solve()
    m3d = _beam_moments_fem3d(b3d)

    np.testing.assert_allclose(m2d, m3d, atol=1e-3)


def test_releasing_axial_both_ends_raises_mechanism():
    mat = ElasticMaterial(E=200e6)
    sec = Section(0.01, 1e-4, material=mat)
    n1 = Node(1, 0.0, 0.0)
    n2 = Node(2, 5.0, 0.0)
    beam = BeamElement("B", n1, n2, mat, sec, releases_i=[True, False, False], releases_j=[True, False, False])
    with pytest.raises(ValueError, match="mechanism"):
        beam.local_stiffness()
