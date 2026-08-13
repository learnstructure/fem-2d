"""
Example 05: Error handling — demonstrate the user-facing error paths
of the adapter functions.

Two kinds of errors are produced by the adapter:

1. ``TypeError`` — when a foreign object (not a ``Structure`` /
   ``SimpleFrame`` / ``Results``) is passed to ``model_to_core`` or
   ``result_to_core``.
2. ``ImportError`` — when ``struct_core`` is not installed. The adapter
   imports ``struct_core`` lazily so the rest of fem2d still works.

Run with::

    python examples/adapters/05_error_handling.py
"""

import sys

from fem2d import SimpleFrame, model_to_core, result_to_core


def build_solved_frame() -> SimpleFrame:
    frame = SimpleFrame()
    frame.add_node(1, 0.0, 0.0)
    frame.add_node(2, 1.0, 0.0)
    frame.add_frame(1, 1, 2, E=200.0e9, A=1.0e-3, I=1.0e-6)
    frame.add_support(1, [True, True, True])
    frame.add_node_load(2, [0.0, -1.0e3, 0.0])
    frame.solve()
    return frame


def demo_typeerror() -> None:
    """Pass a non-Structure to model_to_core; expect a TypeError."""
    print("\n--- TypeError on foreign type ---")
    bogus = {"this": "is not a Structure"}
    try:
        model_to_core(bogus)
    except TypeError as exc:
        msg = str(exc)
        print(f"Caught TypeError: {msg!r}")
        assert msg.startswith("expected a fem2d.Structure"), \
            f"Unexpected TypeError wording: {msg!r}"
        assert "dict" in msg, f"Type name 'dict' should appear in: {msg!r}"
    else:
        raise AssertionError("model_to_core(dict) should have raised TypeError")

    # result_to_core uses the same _resolve_structure helper.
    print("Same TypeError applies to result_to_core:")
    try:
        result_to_core(bogus)
    except TypeError as exc:
        print(f"Caught TypeError: {str(exc)!r}")


def demo_importerror() -> None:
    """Simulate struct_core being uninstalled and confirm the ImportError."""
    print("\n--- ImportError when struct_core is missing ---")

    # Hide struct_core so the adapter's lazy ``from struct_core import``
    # raises ImportError. We must block the *package* entry (setting it to
    # ``None`` is the canonical "this import must fail" sentinel) and pop
    # already-loaded submodules — otherwise ``from struct_core import X``
    # would succeed via the cached submodule entries.
    hidden = []
    for name in list(sys.modules):
        if name == "struct_core" or name.startswith("struct_core."):
            hidden.append((name, sys.modules.pop(name)))
    sys.modules["struct_core"] = None  # force ImportError on next lookup

    try:
        try:
            model_to_core("not even a Structure")
        except ImportError as exc:
            msg = str(exc)
            print(f"Caught ImportError: {msg!r}")
            assert "struct_core is required" in msg, \
                f"Unexpected ImportError wording: {msg!r}"
            assert "fem2d[ecosystem]" in msg, \
                f"Install hint should mention [ecosystem]: {msg!r}"
        else:
            raise AssertionError(
                "model_to_core should have raised ImportError when struct_core is hidden"
            )
    finally:
        # Restore everything so the rest of the script (and other tests)
        # can keep using struct_core normally.
        del sys.modules["struct_core"]
        for name, mod in hidden:
            sys.modules[name] = mod


def main() -> None:
    # Sanity check: a normal call works fine.
    frame = build_solved_frame()
    sc_model = model_to_core(frame)
    print(f"Normal call works: {len(sc_model.elements)} elements exported.")

    demo_typeerror()
    demo_importerror()

    print("\nError-handling demo complete.")


if __name__ == "__main__":
    main()
