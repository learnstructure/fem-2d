from fem2d.elements.beam import BeamElement
import numpy as np


class BeamWithHingesElement(BeamElement):
    def __init__(
        self, eid, node_i, node_j, material, area, inertia, hinge_i=False, hinge_j=False
    ):
        super().__init__(eid, node_i, node_j, material, area, inertia)
        self.hinge_i = hinge_i
        self.hinge_j = hinge_j

    def local_stiffness(self):
        k_full = super().local_stiffness()
        if not self.hinge_i and not self.hinge_j:
            return k_full
        # Apply releases: zero moment at hinge ends.
        # This requires condensing out the rotational DOFs at released ends.
        # For simplicity, we'll implement the common case of a beam with both ends hinged.
        # The resulting stiffness is then a truss‑like element with axial stiffness only,
        # but that's exactly what a truss element does. So for hinges, you might simply use TrussElement.
        # A more sophisticated implementation would keep the bending terms but allow moment release.
        # Here we'll return a modified matrix (for illustration):
        if self.hinge_i and self.hinge_j:
            # Both ends hinged -> axial only
            k = np.zeros((6, 6))
            k[0, 0] = k[0, 3] = k[3, 0] = k[3, 3] = (
                self.area * self.material.E / self.length
            )
            return k
        # Single hinge: you'd need to condense the rotation at that end.
        # This is more involved and usually done by solving for the rotation in terms of translations.
        # For brevity, we'll skip the full implementation here.
        raise NotImplementedError("Single hinge not yet implemented")
