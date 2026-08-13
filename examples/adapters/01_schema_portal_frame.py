"""
Example 01: Build a portal frame using ``struct_core`` Pydantic models and convert
it to a ``fem2d.Structure`` via :func:`model_from_core`.

This is the most common entry-point for users who design their model in
``struct_core`` and want to hand it off to ``fem2d`` for analysis.

Run with::

    python examples/adapters/01_schema_portal_frame.py
"""

from struct_core import (
    BeamElement,
    DistributedLoad,
    ElasticMaterial,
    LinearStaticAnalysis,
    LoadCase,
    Node,
    PointLoad,
    Project,
    RectangularSection,
    StructuralModel,
    Support,
    from_json,
    to_json,
)

from fem2d import model_from_core


def main() -> None:
    print("=== Step 1: Build a StructuralModel via struct_core ===")

    # Portal frame geometry: 6 m span, 4 m tall.
    model = StructuralModel(
        nodes=[
            Node(id=1, x=0.0, y=0.0),   # base left
            Node(id=2, x=0.0, y=4.0),   # top left
            Node(id=3, x=6.0, y=4.0),   # top right
            Node(id=4, x=6.0, y=0.0),   # base right
        ],
        materials=[
            ElasticMaterial(id="StructuralSteel", E=200e9, rho=7850.0),
        ],
        sections=[
            RectangularSection(id="ColSec", b=0.30, h=0.30),   # columns
            RectangularSection(id="BeamSec", b=0.25, h=0.45),  # beam
        ],
        elements=[
            BeamElement(id=1, start_node=1, end_node=2,
                        material_id="StructuralSteel", section_id="ColSec"),
            BeamElement(id=2, start_node=2, end_node=3,
                        material_id="StructuralSteel", section_id="BeamSec"),
            BeamElement(id=3, start_node=3, end_node=4,
                        material_id="StructuralSteel", section_id="ColSec"),
        ],
        supports=[
            Support.fixed(node_id=1, name="Left Support"),
            Support.fixed(node_id=4, name="Right Support"),
        ],
        load_cases=[
            LoadCase(
                id="GravityAndWind",
                name="Dead Load + Lateral Wind",
                point_loads=[
                    PointLoad(node_id=2, fx=20e3),  # 20 kN lateral wind at top-left
                ],
                element_loads=[
                    DistributedLoad(element_id=2, wy=-15e3),  # 15 kN/m on beam
                ],
            )
        ],
        analysis_cases=[
            LinearStaticAnalysis(id="LinearSolve", load_case_id="GravityAndWind"),
        ],
    )

    # Round-trip through JSON to demonstrate the schema is fully serializable.
    project = Project()
    project.metadata.title = "2D Portal Frame Example"
    project.metadata.project_name = "Industrial Warehouse"
    project.metadata.author = "Structural Engineer AI"
    project.model = model

    json_repr = to_json(project, indent=2)
    print(f"Project serialized to JSON successfully! ({len(json_repr)} chars)")
    reloaded_project = from_json(Project, json_repr)

    print("\n=== Step 2: Convert schema to fem2d.Structure via model_from_core ===")
    structure = model_from_core(reloaded_project)
    print(f"Created Structure with {len(structure.nodes)} nodes "
          f"and {len(structure.elements)} elements.")

    print("\n=== Step 3: Solve linear static analysis ===")
    structure.solve()

    print("\n=== Step 4: Nodal displacements ===")
    for nid, node in sorted(structure.nodes.items()):
        dofs = node.dofs
        ux = structure.disp[dofs[0]]
        uy = structure.disp[dofs[1]]
        rz = structure.disp[dofs[2]]
        print(
            f"Node {nid:2d} (x={node.x:3.1f} m, y={node.y:3.1f} m): "
            f"Ux = {ux * 1e3:7.3f} mm, "
            f"Uy = {uy * 1e3:7.3f} mm, "
            f"Rz = {rz * 1e3:7.3f} mrad/1000"
        )

    print("\nAdapter integrated seamlessly with fem2d.")


if __name__ == "__main__":
    main()
