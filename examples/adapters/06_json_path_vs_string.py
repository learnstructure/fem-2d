"""
Example 06: JSON round-trip — path-based vs in-memory string-based.

``struct_core`` exposes two equivalent pairs of JSON helpers:

* ``save_json(obj, path)``  / ``load_json(cls, path)``  — reads/writes a file.
* ``to_json(obj)``           / ``from_json(cls, json_str)`` — string-based.

Both pairs produce the same Pydantic model; the choice is purely about
where you want the JSON to live (on disk vs in memory, e.g. when the JSON
travels through an HTTP request body).

This example shows a full fem2d → JSON → fem2d round-trip using both
pairs side-by-side and asserts they yield identical results.

Run with::

    python examples/adapters/06_json_path_vs_string.py
"""

import os
import tempfile

from struct_core import (
    Project,
    load_json,
    save_json,
    from_json,
    to_json,
)

from fem2d import SimpleFrame, model_from_core, model_to_core


def build_frame() -> SimpleFrame:
    frame = SimpleFrame()
    frame.add_node(1, 0.0, 0.0)
    frame.add_node(2, 0.0, 2.0)
    frame.add_frame(1, 1, 2, E=200.0e9, A=1.0e-3, I=1.0e-6)
    frame.add_support(1, [True, True, True])
    frame.add_node_load(2, [5.0e3, 0.0, 0.0])
    frame.solve()
    return frame


def wrap_in_project(model) -> Project:
    project = Project()
    project.metadata.title = "JSON Round-Trip Demo"
    project.metadata.project_name = "fem2d/adapters example 06"
    project.metadata.author = "fem2d examples"
    project.model = model
    return project


def main() -> None:
    frame = build_frame()
    sc_model = model_to_core(frame)
    original_project = wrap_in_project(sc_model)

    # --- Path-based pair ----------------------------------------------------
    tmpdir = tempfile.mkdtemp(prefix="fem2d_adapter_06_")
    path = os.path.join(tmpdir, "model.json")
    print(f"Writing to disk: {path}")
    save_json(original_project, path)
    print(f"File size on disk: {os.path.getsize(path)} bytes")

    reloaded_path = load_json(Project, path)
    print(f"Loaded via load_json(Project, path): {type(reloaded_path).__name__}")

    # --- String-based pair --------------------------------------------------
    json_str = to_json(original_project, indent=2)
    print(f"Serialized to in-memory JSON: {len(json_str)} chars")
    reloaded_str = from_json(Project, json_str)
    print(f"Loaded via from_json(Project, json_str): {type(reloaded_str).__name__}")

    # --- Semantic equality --------------------------------------------------
    # The two reloads should yield Pydantic-equal projects.
    assert reloaded_path.model_dump() == reloaded_str.model_dump(), (
        "Path-based and string-based reloads disagree"
    )
    print("Path-based and string-based reloads are semantically identical.")

    # --- Round-trip into fem2d and re-solve --------------------------------
    structure_from_path = model_from_core(reloaded_path)
    structure_from_path.solve()

    structure_from_str = model_from_core(reloaded_str)
    structure_from_str.solve()

    if not (structure_from_path.disp == structure_from_str.disp).all():
        raise AssertionError("fem2d displacement differs between reload paths.")
    if not (structure_from_path.disp == frame.structure.disp).all():
        raise AssertionError("fem2d displacement differs from original after JSON round-trip.")

    print("Both reload paths yield fem2d displacements equal to the original.")

    # Cleanup
    os.remove(path)
    os.rmdir(tmpdir)
    print(f"Cleaned up {path}.")


if __name__ == "__main__":
    main()
