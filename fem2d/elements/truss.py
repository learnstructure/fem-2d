from fem2d.elements.element import ElementBase
import numpy as np


class TrussElement(ElementBase):
    def __init__(self, eid, node_i, node_j, material, area):
        super().__init__(eid, node_i, node_j)
        self.material = material
        self.area = area

    def local_stiffness(self):
        k = self.area * self.material.E / self.length
        return np.array([[k, 0, -k, 0], [0, 0, 0, 0], [-k, 0, k, 0], [0, 0, 0, 0]])

    def global_stiffness(self):
        # Expand 4x4 to 6x6 (insert axial DOFs at positions 0,1 and 3,4)
        k_local = self.local_stiffness()
        T4 = self.transformation_matrix(dof_per_node=2)  # custom 4x4 T
        k_global_4 = T4.T @ k_local @ T4
        k_global = np.zeros((6, 6))
        # Map: local [u1,v1,u2,v2] -> global indices [0,1,3,4] of 6x6
        dof_map = [0, 1, 3, 4]
        for i, ii in enumerate(dof_map):
            for j, jj in enumerate(dof_map):
                k_global[ii, jj] = k_global_4[i, j]
        return k_global

    def transformation_matrix(self, dof_per_node=2):
        c = self.cos
        s = self.sin
        return np.array([[c, s, 0, 0], [-s, c, 0, 0], [0, 0, c, s], [0, 0, -s, c]])
