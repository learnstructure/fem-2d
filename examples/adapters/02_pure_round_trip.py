"""
Example 02: Pure in-memory round-trip — no JSON involved.

Build a small cantilever in fem2d, export it via :func:`model_to_core`,
then re-import it via :func:`model_from_core`, solve both models, and
assert the displacement fields are equal.

We also verify the analytical tip-deflection formula::

    delta = P L^3 / (3 E I)

which provides a physics-level sanity check that the round-trip is
preserving both the geometry and the structural properties.

Run with::

    python examples/adapters/02_pure_round_trip.py
"""

from fem2d import SimpleFrame, model_from_core, model_to_core


# Cantilever beam: 2 m long, tip-loaded with P = -1 kN.
L = 2.0            # m
P = -1.0e3         # N (downward)
E = 200.0e9        # Pa
A = 1.0e-3         # m^2
I = 1.0e-6         # m^4


def build_frame() -> SimpleFrame:
    """Build a horizontal cantilever with a tip point load."""
    frame = SimpleFrame()
    frame.add_node(1, 0.0, 0.0)
    frame.add_node(2, L,  0.0)
    frame.add_frame(1, 1, 2, E, A, I)

    frame.add_support(1, [True, True, True])      # fully fixed at root
    frame.add_node_load(2, [0.0, P, 0.0])         # tip load
    return frame


def main() -> None:
    # --- Solve on the original fem2d model ---------------------------------
    frame = build_frame()
    frame.solve()
    original_disp = frame.structure.disp.copy()
    print(f"Original fem2d model solved. disp vector size = {len(original_disp)}")

    # --- Export to struct_core and import back -----------------------------
    sc_model = model_to_core(frame)
    print("Exported fem2d model -> struct_core.StructuralModel "
          f"({len(sc_model.materials)} material(s), {len(sc_model.sections)} section(s), "
          f"{len(sc_model.elements)} element(s))")

    round_tripped = model_from_core(sc_model)
    round_tripped.solve()
    rt_disp = round_tripped.disp
    print(f"Round-tripped fem2d model solved. disp vector size = {len(rt_disp)}")

    # --- Round-trip equality ------------------------------------------------
    if not (original_disp == rt_disp).all():
        diff = (original_disp - rt_disp)
        raise AssertionError(
            f"Round-trip mismatch: max |Δu| = {abs(diff).max():.3e}"
        )
    print("Round-trip equality OK: disp vectors are identical.")

    # --- Physics-level sanity check ----------------------------------------
    # Cantilever tip-deflection: delta = P L^3 / (3 E I)
    # Take the y-DOF of node 2 from the round-tripped model.
    node2 = round_tripped.nodes[2]
    uy_tip = round_tripped.disp[node2.dofs[1]]
    uy_expected = P * L ** 3 / (3.0 * E * I)

    print(f"\nTip Uy (fem2d round-tripped) = {uy_tip:+.6e} m")
    print(f"Tip Uy (analytical P L^3/3EI) = {uy_expected:+.6e} m")

    rel = abs(uy_tip - uy_expected) / abs(uy_expected)
    print(f"Relative error = {rel:.2e}")
    assert rel < 1e-9, f"Tip deflection off by {rel:.2e}"

    print("\nRound-trip is lossless AND numerically correct.")


if __name__ == "__main__":
    main()
