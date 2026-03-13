from fem2d.elements.element import ElementBase
import numpy as np

from fem2d.elements.truss import TrussElement


class SpringElement(TrussElement):
    def __init__(self, eid, node_i, node_j, stiffness):
        # We don't need material/area; just pass a dummy material or handle separately.
        self.k = stiffness
        super().__init__(eid, node_i, node_j, material=None, area=None)

    def local_stiffness(self):
        return self.k * np.array(
            [[1, 0, -1, 0], [0, 0, 0, 0], [-1, 0, 1, 0], [0, 0, 0, 0]]
        )
