"""
Example 04: ID stability — ``model_to_core`` assigns ``M1/M2/…`` and
``S1/S2/…`` based on first-encountered uniqueness, and the assignment
must be stable across repeated exports of the same fem2d model.

This matters because downstream consumers (visualizers, code-checkers)
key off those IDs; if the assignment drifts between calls, downstream
artefacts would churn even when the underlying model didn't change.

Run with::

    python examples/adapters/04_id_stability.py
"""

from fem2d import Section, SimpleFrame, model_to_core
from fem2d.elements.truss import TrussElement
from fem2d.materials import ElasticMaterial


def build_frame() -> SimpleFrame:
    """Portal frame with two distinct materials and two distinct sections."""
    frame = SimpleFrame()
    frame.add_node(1, 0.0, 0.0)
    frame.add_node(2, 0.0, 4.0)
    frame.add_node(3, 6.0, 4.0)
    frame.add_node(4, 6.0, 0.0)

    # Columns: A36 steel (200 GPa, lighter column section)
    # Beam:    A992 steel (200 GPa, deeper beam section)
    E_col, A_col, I_col = 200.0e9, 30.0e-4, 30.0e-8
    E_bm,  A_bm,  I_bm  = 200.0e9, 12.5e-4, 22.5e-8

    frame.add_frame(1, 1, 2, E_col, A_col, I_col)  # left column  (M1, S1)
    frame.add_frame(2, 2, 3, E_bm,  A_bm,  I_bm)   # beam        (M1, S2)
    frame.add_frame(3, 3, 4, E_col, A_col, I_col)  # right column (M1, S1)

    frame.add_support(1, [True, True, True])
    frame.add_support(4, [True, True, True])
    frame.add_node_load(2, [20.0e3, 0.0, 0.0])    # lateral wind
    frame.add_node_load(3, [0.0, -30.0e3, 0.0])   # vertical gravity
    return frame


def material_section_ids(model) -> tuple[list[str], list[str]]:
    mats = [m.id for m in model.materials]
    secs = [s.id for s in model.sections]
    return mats, secs


def main() -> None:
    frame = build_frame()
    frame.solve()

    # --- First export ------------------------------------------------------
    sc_model_1 = model_to_core(frame)
    mats_1, secs_1 = material_section_ids(sc_model_1)
    print(f"Export 1 materials: {mats_1}")
    print(f"Export 1 sections : {secs_1}")

    # --- Second export -----------------------------------------------------
    sc_model_2 = model_to_core(frame)
    mats_2, secs_2 = material_section_ids(sc_model_2)
    print(f"Export 2 materials: {mats_2}")
    print(f"Export 2 sections : {secs_2}")

    assert mats_1 == mats_2, f"Material IDs drifted: {mats_1} -> {mats_2}"
    assert secs_1 == secs_2, f"Section IDs drifted:  {secs_1} -> {secs_2}"

    # --- Property assertions ---------------------------------------------
    # Only one material was used (E=200 GPa everywhere), so we expect M1.
    assert mats_1 == ["M1"], f"Expected single material 'M1', got {mats_1}"

    # Two distinct sections were used: column sec first (element 1), beam
    # sec second (element 2). So we expect ['S1', 'S2'].
    assert secs_1 == ["S1", "S2"], f"Expected ['S1', 'S2'], got {secs_1}"

    # And verify that the element-to-section mapping is what we expect:
    elem_to_sec = {e.id: e.section_id for e in sc_model_1.elements}
    print(f"Element -> section_id mapping: {elem_to_sec}")
    assert elem_to_sec[1] == "S1"   # left column
    assert elem_to_sec[2] == "S2"   # beam
    assert elem_to_sec[3] == "S1"   # right column

    # --- Now build a frame with TWO materials to exercise M1/M2 ------------
    # Use the low-level Structure API for the second model so the truss
    # brace can be given a tiny non-zero moment of inertia (struct_core
    # requires Iz > 0). SimpleFrame.add_truss() sets I=0 by construction.
    from fem2d.nodes import Node
    from fem2d.loads import PointLoad
    from fem2d import Structure

    frame2 = Structure()
    for nid, x, y in [(1, 0.0, 0.0), (2, 0.0, 2.0), (3, 2.0, 0.0)]:
        frame2.add_node(Node(nid, x, y))
    n_root, n_tip, n_anc = frame2.nodes[1], frame2.nodes[2], frame2.nodes[3]
    n_root.set_support(ux_fixed=True, uy_fixed=True, rz_fixed=True)
    n_anc.set_support(ux_fixed=True, uy_fixed=True, rz_fixed=True)

    steel     = ElasticMaterial(E=200.0e9)
    aluminium = ElasticMaterial(E=70.0e9)
    col_sec   = Section(area=1.0e-3, moi=1.0e-6)
    # struct_core.GeneralSection requires Iz > 0; pass a tiny positive I.
    truss_sec = Section(area=2.0e-4, moi=1.0e-12)

    from fem2d.elements.beam import BeamElement
    frame2.add_element(BeamElement(1, n_root, n_tip, steel, col_sec))
    frame2.add_element(TrussElement(2, n_tip, n_anc, aluminium, truss_sec))
    frame2.add_load(PointLoad(node=n_tip, fx=1.0e3, fy=0.0, mz=0.0))
    frame2.solve()

    sc_2mat = model_to_core(frame2)
    mats_2mat, _ = material_section_ids(sc_2mat)
    print(f"\nTwo-material frame materials: {mats_2mat}")
    assert mats_2mat == ["M1", "M2"], f"Expected ['M1', 'M2'], got {mats_2mat}"

    # Element 1 (column) should reference the first-encountered material.
    assert sc_2mat.elements[0].material_id == "M1"
    # Element 2 (aluminium truss) should reference the second material.
    assert sc_2mat.elements[1].material_id == "M2"

    print("\nID assignment is stable across exports and follows first-encounter order.")


if __name__ == "__main__":
    main()
