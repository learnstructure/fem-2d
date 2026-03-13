from fem2d.elements.element import ElementBase
import numpy as np


class BeamElement(ElementBase):
    def __init__(self, eid, node_i, node_j, material, area, inertia):
        super().__init__(eid, node_i, node_j)
        self.material = material
        self.area = area
        self.inertia = inertia

    def local_stiffness(self):
        E = self.material.E
        A = self.area
        I = self.inertia
        L = self.length

        EA_L = E * A / L
        EI_L3 = E * I / L**3
        EI_L2 = E * I / L**2  # = EI / L²
        EI_L = E * I / L

        return np.array(
            [
                [EA_L, 0, 0, -EA_L, 0, 0],
                [0, 12 * EI_L3, 6 * EI_L2, 0, -12 * EI_L3, 6 * EI_L2],
                [0, 6 * EI_L2, 4 * EI_L, 0, -6 * EI_L2, 2 * EI_L],
                [-EA_L, 0, 0, EA_L, 0, 0],
                [0, -12 * EI_L3, -6 * EI_L2, 0, 12 * EI_L3, -6 * EI_L2],
                [0, 6 * EI_L2, 2 * EI_L, 0, -6 * EI_L2, 4 * EI_L],
            ]
        )

    # transformation_matrix() from base class is already 6x6, so global_stiffness() works directly.

    def equivalent_nodal_loads(self):
        # If a distributed load is stored, compute fixed‑end forces.
        # For example, uniform load w in local y‑direction:
        if hasattr(self, "w"):
            w = self.w
            L = self.length
            return np.array([0, w * L / 2, w * L**2 / 12, 0, w * L / 2, -w * L**2 / 12])
        return np.zeros(6)
