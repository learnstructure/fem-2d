import numpy as np


class PointLoad:
    def __init__(self, node, fx=0, fy=0, mz=0):
        self.node = node
        self.fx, self.fy, self.mz = fx, fy, mz


class DistributedLoad:
    """
    Represents a uniformly distributed load acting on an element.
    The load is defined in the element's local coordinate system:
        wx: force per unit length along local x (axial)
        wy: force per unit length along local y (transverse)
    """

    def __init__(self, element, wx=0.0, wy=0.0):
        self.element = element
        self.wx = wx
        self.wy = wy
        self._compute_equivalent_loads()

    def _compute_equivalent_loads(self):
        """Compute equivalent nodal loads in global coordinates and store in element.eq_load."""
        L = self.element.length
        # Standard fixed‑end forces for uniform load (local coordinates)
        eq_local = np.array(
            [
                self.wx * L / 2,  # Fx_i
                self.wy * L / 2,  # Fy_i
                self.wy * L**2 / 12,  # M_i
                self.wx * L / 2,  # Fx_j
                self.wy * L / 2,  # Fy_j
                -self.wy * L**2 / 12,  # M_j
            ]
        )
        # Transform to global using the element's initial orientation
        c = self.element.cos
        s = self.element.sin
        T = np.array(
            [
                [c, -s, 0, 0, 0, 0],
                [s, c, 0, 0, 0, 0],
                [0, 0, 1, 0, 0, 0],
                [0, 0, 0, c, -s, 0],
                [0, 0, 0, s, c, 0],
                [0, 0, 0, 0, 0, 1],
            ]
        )
        eq_global = T @ eq_local
        # eq_global = eq_local
        self.element.eq_load = eq_global
