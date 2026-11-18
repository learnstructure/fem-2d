"""fem2d is a plane-frame solver: a model with shell elements is refused, never analyzed without them."""

import pytest

from struct_core import (
    ConcreteMaterial,
    Node,
    ShellElement,
    ShellSection,
    StructuralModel,
    Support,
)

from fem2d.adapters import model_from_core


def test_model_with_shells_is_refused():
    model = StructuralModel(
        units="kN, m",
        nodes=[Node(id=i + 1, x=x, y=y) for i, (x, y) in enumerate([(0, 0), (1, 0), (1, 1), (0, 1)])],
        materials=[ConcreteMaterial(id="C", E=25e6, fc=30e3)],
        shell_sections=[ShellSection(id="S", material_id="C", thickness=0.2, behavior="membrane")],
        shell_elements=[ShellElement(id="P1", node_ids=[1, 2, 3, 4], section_id="S")],
        supports=[Support(node_id=1, ux=True, uy=True, rz=True)],
    )
    with pytest.raises(ValueError, match="does not support shell elements.*fem3d"):
        model_from_core(model)
