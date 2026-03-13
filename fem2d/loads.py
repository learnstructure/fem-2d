import numpy as np


class PointLoad:
    def __init__(self, node, fx=0, fy=0, mz=0):
        self.node = node
        self.fx, self.fy, self.mz = fx, fy, mz


class DistributedLoad:
    def __init__(self, element, wx=0, wy=0):
        self.element = element
        self.wx, self.wy = wx, wy  # components in local axes

    def equivalent_nodal_loads(self):
        # Compute fixed‑end forces for uniform wy
        L = self.element.length
        return np.array(
            [
                0,
                self.wy * L / 2,
                self.wy * L**2 / 12,
                0,
                self.wy * L / 2,
                -self.wy * L**2 / 12,
            ]
        )
