"""
Example 03: Polymorphism — adapter accepts raw ``fem2d.Structure``,
``SimpleFrame``, and ``Results`` wrappers.

We build a simple truss + beam hybrid directly on a ``fem2d.Structure``
(no ``SimpleFrame`` sugar) and round-trip it through ``struct_core``.
This demonstrates that:

* ``model_to_core`` works on any object exposing a ``.structure`` attr
  (or a bare ``Structure``), and
* Element subclasses (``BeamElement``, ``TrussElement``, ``SpringElement``)
  collapse onto the appropriate ``struct_core`` element class.

Run with::

    python examples/adapters/03_structure_polymorphism.py
"""

from fem2d import (
    ElasticMaterial,
    Results,
    Section,
    SimpleFrame,
    Structure,
    model_from_core,
    model_to_core,
    result_to_core,
)
from fem2d.elements.beam import BeamElement as FemBeamElement
from fem2d.elements.truss import TrussElement
from fem2d.loads import PointLoad as FemPointLoad
from fem2d.nodes import Node
from struct_core import (
    BeamElement as ScBeamElement,
    TrussElement as ScTrussElement,
)


def build_raw_structure() -> Structure:
    """Build a small cantilever with a truss brace directly on Structure."""
    structure = Structure()

    # Materials & sections
    steel = ElasticMaterial(E=200e9, rho=7850.0)
    beam_sec = Section.from_rectangle(width=0.05, depth=0.20)   # 50x200 mm
    # Truss section: pin-ended axial-only. struct_core.GeneralSection
    # requires Iz > 0, so we pass a tiny non-zero moment of inertia.
    truss_sec = Section(area=2.0e-4, moi=1.0e-12)

    # Geometry: a vertical cantilever with a diagonal truss brace
    n_root = Node(nid=1, x=0.0, y=0.0)
    n_tip  = Node(nid=2, x=0.0, y=2.0)
    n_anc  = Node(nid=3, x=1.0, y=0.0)

    n_root.set_support(ux_fixed=True, uy_fixed=True, rz_fixed=True)
    n_anc.set_support(ux_fixed=True, uy_fixed=True, rz_fixed=True)
    # n_tip is free.

    structure.add_node(n_root)
    structure.add_node(n_tip)
    structure.add_node(n_anc)

    beam = FemBeamElement(eid=1, node_i=n_root, node_j=n_tip,
                          material=steel, section=beam_sec)
    truss = TrussElement(eid=2, node_i=n_tip, node_j=n_anc,
                         material=steel, section=truss_sec)

    structure.add_element(beam)
    structure.add_element(truss)

    # Tip load
    structure.add_load(FemPointLoad(node=n_tip, fx=5.0e3, fy=0.0, mz=0.0))
    return structure


def main() -> None:
    # --- 1. Raw Structure --------------------------------------------------
    structure = build_raw_structure()
    structure.solve()
    print(f"Raw Structure solved: {len(structure.elements)} elements, "
          f"{len(structure.nodes)} nodes.")

    sc_model = model_to_core(structure)
    types = sorted({type(e).__name__ for e in sc_model.elements})
    print(f"struct_core element types produced: {types}")
    # We should see both a BeamElement and a TrussElement in the schema.
    assert any(isinstance(e, ScBeamElement) for e in sc_model.elements), \
        "Expected at least one struct_core.BeamElement"
    assert any(isinstance(e, ScTrussElement) for e in sc_model.elements), \
        "Expected at least one struct_core.TrussElement"

    # --- 2. SimpleFrame (polymorphic input) --------------------------------
    # Beam + diagonal truss: the beam gives bending stiffness at the free
    # tip so the K-matrix is non-singular; the truss element demonstrates
    # element-type polymorphism in the exported struct_core schema.
    frame = SimpleFrame()
    frame.add_node(1, 0.0, 0.0)
    frame.add_node(2, 0.0, 2.0)
    frame.add_node(3, 1.0, 0.0)
    frame.add_frame(1, 1, 2, E=200e9, A=0.01,   I=3.333e-5)  # vertical beam
    # struct_core.GeneralSection requires Iz > 0; pass a tiny non-zero I.
    frame.add_truss(2, 2, 3, E=200e9, A=2.0e-4, I=1.0e-12)   # diagonal brace
    frame.add_support(1, [True, True, True])
    frame.add_support(3, [True, True, True])
    frame.add_node_load(2, [5.0e3, 0.0, 0.0])
    frame.solve()

    sc_from_frame = model_to_core(frame)
    print(f"model_to_core(SimpleFrame) -> {len(sc_from_frame.elements)} elements, "
          f"types = {sorted({type(e).__name__ for e in sc_from_frame.elements})}")

    # --- 3. Results wrapper (also accepted by model_to_core) ---------------
    results = Results(frame)
    sc_from_results = model_to_core(results)
    assert sc_from_results.nodes == sc_from_frame.nodes, \
        "Results and SimpleFrame should yield equivalent schemas"
    print("model_to_core(Results) yielded equivalent schema to model_to_core(SimpleFrame).")

    # --- 4. result_to_core accepts all three flavours ----------------------
    res_a = result_to_core(structure, analysis_case_id="Static")
    res_b = result_to_core(frame,     analysis_case_id="Static")
    res_c = result_to_core(results,   analysis_case_id="Static")
    print(f"result_to_core produced AnalysisResult with "
          f"{len(res_a.node_results)} node_results / {len(res_a.element_results)} element_results")
    assert len(res_a.node_results) == len(res_b.node_results) == len(res_c.node_results)

    # --- 5. model_from_core round-trip on raw Structure --------------------
    rt = model_from_core(sc_model)
    rt.solve()
    if not (rt.disp == structure.disp).all():
        raise AssertionError("Round-trip via raw Structure changed displacements.")
    print("Round-trip via raw Structure: disp vectors identical.")

    print("\nPolymorphism verified: Structure / SimpleFrame / Results all accepted.")


if __name__ == "__main__":
    main()
