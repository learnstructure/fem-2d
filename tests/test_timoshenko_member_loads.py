"""
Equivalent nodal loads of member loads on beam elements, with and without shear deformation.

The equivalent nodal loads are the work of the member load on the exact deflected shapes of the element
(``fem2d.loads.bending_shape_functions``: Timoshenko shapes, Hermite cubics without shear deformation). Checks:

1. Fixed-end forces from the flexibility (unit-load) method, which uses neither shape functions nor the stiffness
   matrix: a cantilever under the load, with the end force and moment as redundants, and the flexibilities
   integral(M m / EI) + integral(V v / kGA) by Gauss quadrature on each load segment (exact for these polynomials).
2. A member point force or moment gives the same solution as a mesh with a node at the load point, where the load
   is a nodal load (exact for any element whose stiffness is exact), with and without an end release.
3. The closed-form Euler-Bernoulli fixed-end forces of a concentrated moment M at a (the slopes of the Hermite
   cubics at a): moment reactions M b (2a - b) / L^2 at i and M a (2b - a) / L^2 at j, end shears 6 M a b / L^3.

Kassimali (2021), Matrix Analysis of Structures, 3rd ed., Section 9.7, notes that the usual fixed-end force
expressions do not include shear deformations, and that shear-deformable ones follow from the same derivation.
"""

import numpy as np
import pytest

from fem2d import BeamElement, ElasticMaterial, Node, Results, Section, Structure
from fem2d.loads import DistributedLoad, ElementPointLoad, PointLoad, TriangularLoad, bending_shape_functions

E, NU = 25.0e6, 0.2
L = 2.0          # short, deep member: phi = 0.46
WIDTH, DEPTH = 0.4, 0.8


def _beam(shear=True, releases_j=None, node_i=None, node_j=None, eid=1):
    material = ElasticMaterial(E=E, nu=NU)
    section = Section.from_rectangle(width=WIDTH, depth=DEPTH, material=material)
    return BeamElement(eid, node_i or Node(1, 0.0, 0.0), node_j or Node(2, L, 0.0), material=material,
                       section=section, releases_j=releases_j, include_shear_deformation=shear)


def _eq(load):
    """Local equivalent nodal loads [u_i, v_i, theta_i, u_j, v_j, theta_j] stored by a load on its element."""
    return np.asarray(load.element._raw_eq_load_local)


def fixed_end_loads(length, ei, kga, distributed=(), forces=(), moments=()):
    """Equivalent nodal loads (v_i, theta_i, v_j, theta_j) of a fixed-fixed beam, by the flexibility method.

    ``distributed``: (w_start, w_end, start, end) with relative start and end; ``forces``: (P, a) transverse;
    ``moments``: (m, a) counterclockwise; a is the distance from end i. ``kga = None``: no shear deformation.
    The beam is released at j (a cantilever); the end force R and moment Mj at j follow from v_j = theta_j = 0 with
    the unit-load method. M(x) and V(x) are the moment and force resultants of the actions on [x, L].
    """
    gauss_x, gauss_w = np.polynomial.legendre.leggauss(12)
    breaks = sorted({0.0, length, *[s * length for _, _, s, _ in distributed],
                     *[e * length for _, _, _, e in distributed], *[a for _, a in forces], *[a for _, a in moments]})

    def integrate(f, lo, hi):
        total = 0.0
        pts = [lo] + [b for b in breaks if lo < b < hi] + [hi]
        for a, b in zip(pts[:-1], pts[1:]):
            total += 0.5 * (b - a) * sum(w * f(0.5 * (b - a) * t + 0.5 * (a + b)) for t, w in zip(gauss_x, gauss_w))
        return total

    def q(s):
        total = 0.0
        for w0, w1, st, en in distributed:
            a, b = st * length, en * length
            if a <= s <= b:
                total += w0 + (w1 - w0) * (s - a) / (b - a)
        return total

    def shear0(x):
        return integrate(q, x, length) + sum(p for p, a in forces if a > x)

    def moment0(x):
        return (integrate(lambda s: q(s) * (s - x), x, length) + sum(p * (a - x) for p, a in forces if a > x)
                + sum(m for m, a in moments if a > x))

    flex_v = (lambda x: 0.0) if kga is None else (lambda x: 1.0 / kga)
    f11 = integrate(lambda x: (length - x) ** 2 / ei + flex_v(x), 0.0, length)
    f12 = integrate(lambda x: (length - x) / ei, 0.0, length)
    f22 = length / ei
    d1 = integrate(lambda x: moment0(x) * (length - x) / ei + shear0(x) * flex_v(x), 0.0, length)
    d2 = integrate(lambda x: moment0(x) / ei, 0.0, length)
    r_j, m_j = np.linalg.solve([[f11, f12], [f12, f22]], [-d1, -d2])
    total_force = integrate(q, 0.0, length) + sum(p for p, _ in forces)
    total_moment_i = (integrate(lambda s: q(s) * s, 0.0, length) + sum(p * a for p, a in forces)
                      + sum(m for m, _ in moments))
    r_i = -(total_force + r_j)
    m_i = -(total_moment_i + r_j * length + m_j)
    return np.array([-r_i, -m_i, -r_j, -m_j])  # equivalent nodal loads = - fixed-end reactions


def _stiffness_terms(beam, shear):
    ei = E * beam.section.I
    kga = beam.material.G * beam.section.As if shear else None
    return ei, kga


@pytest.mark.parametrize("shear", [True, False])
@pytest.mark.parametrize("w0, w1, start, end", [(-12.0, -4.0, 0.15, 0.7), (5.0, 5.0, 0.0, 1.0), (-8.0, 0.0, 0.4, 1.0)])
def test_distributed_load_matches_flexibility_method(shear, w0, w1, start, end):
    beam = _beam(shear)
    eq = _eq(DistributedLoad(beam, wx=1.0, wy=w0, wy_end=w1, start=start, end=end))
    ref = fixed_end_loads(L, *_stiffness_terms(beam, shear), distributed=[(w0, w1, start, end)])
    np.testing.assert_allclose(eq[[1, 2, 4, 5]], ref, rtol=1e-10, atol=1e-10)


@pytest.mark.parametrize("shear", [True, False])
def test_triangular_load_matches_flexibility_method(shear):
    beam = _beam(shear)
    eq = _eq(TriangularLoad(beam, wy1=-3.0, wy2=-11.0))
    ref = fixed_end_loads(L, *_stiffness_terms(beam, shear), distributed=[(-3.0, -11.0, 0.0, 1.0)])
    np.testing.assert_allclose(eq[[1, 2, 4, 5]], ref, rtol=1e-10, atol=1e-10)


@pytest.mark.parametrize("shear", [True, False])
@pytest.mark.parametrize("pos, force, moment", [(0.3, -15.0, 7.0), (0.0, 10.0, -4.0), (0.5, 0.0, 9.0), (1.0, 6.0, 2.0)])
def test_point_force_and_moment_match_flexibility_method(shear, pos, force, moment):
    beam = _beam(shear)
    a = pos * L
    eq = _eq(ElementPointLoad(beam, py=force, mz=moment, x=a))
    ref = fixed_end_loads(L, *_stiffness_terms(beam, shear), forces=[(force, a)], moments=[(moment, a)])
    np.testing.assert_allclose(eq[[1, 2, 4, 5]], ref, rtol=1e-10, atol=1e-10)


def test_concentrated_moment_closed_form_fixed_end_forces():
    """Regression: a member moment was split between the end rotations by the lever rule, with no end shears."""
    m, a = 9.0, 0.6
    b = L - a
    eq = _eq(ElementPointLoad(_beam(shear=False), mz=m, x=a))
    expected = [0.0, -6.0 * m * a * b / L**3, m * b * (b - 2.0 * a) / L**2,
                0.0, 6.0 * m * a * b / L**3, m * a * (a - 2.0 * b) / L**2]
    np.testing.assert_allclose(eq, expected, rtol=1e-12, atol=1e-12)


def test_shape_functions_take_unit_end_values():
    for phi in (0.0, 0.46, 3.0):
        start, end = bending_shape_functions(0.0, L, phi), bending_shape_functions(L, L, phi)
        values = np.array([start[0], start[1], end[0], end[1]])  # rows: v_i, theta_i, v_j, theta_j
        np.testing.assert_allclose(values, np.eye(4), atol=1e-14)


def _solve(with_mesh, releases_j, force=-15.0, moment=7.0):
    structure = Structure()
    a, b = Node("A", 0.0, 0.0), Node("B", L, 0.0)
    nodes = [a, b] + ([Node("P", 0.3 * L, 0.0)] if with_mesh else [])
    for n in nodes:
        structure.add_node(n)
    a.set_support(ux_fixed=True, uy_fixed=True, rz_fixed=True)
    b.set_support(uy_fixed=True, rz_fixed=True)
    if with_mesh:
        structure.add_element(_beam(node_i=a, node_j=nodes[2], eid=1))
        structure.add_element(_beam(node_i=nodes[2], node_j=b, eid=2, releases_j=releases_j))
        structure.add_load(PointLoad(nodes[2], fx=2.0, fy=force, mz=moment))
    else:
        beam = _beam(node_i=a, node_j=b, releases_j=releases_j)
        structure.add_element(beam)
        structure.add_load(ElementPointLoad(beam, px=2.0, py=force, mz=moment, x=0.3 * L))
    structure.solve()
    results = Results(structure)
    reactions = results.reactions().set_index("node").loc[["A", "B"], ["Fx", "Fy", "M"]].to_numpy()
    ux_b = results.node_displacements().set_index("node").loc["B", "ux"]
    return reactions, ux_b


@pytest.mark.parametrize("releases_j", [None, [False, False, True]])
def test_member_point_load_equals_nodal_load_on_a_mesh(releases_j):
    reactions_member, ux_member = _solve(False, releases_j)
    reactions_mesh, ux_mesh = _solve(True, releases_j)
    np.testing.assert_allclose(reactions_member, reactions_mesh, rtol=1e-9, atol=1e-9)
    assert ux_member == pytest.approx(ux_mesh, rel=1e-9)
