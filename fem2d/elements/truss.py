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

    def get_local_forces(self):
        """Return local end forces as a 6‑component array: [Fx_i, Fy_i, M_i, Fx_j, Fy_j, M_j]."""
        # Get global displacements
        u_i = self.structure.disp[self.node_i.dofs]  # [ux, uy, rz]
        u_j = self.structure.disp[self.node_j.dofs]
        u_global = np.concatenate([u_i, u_j])  # 6×1

        # Transform to local (using 4×4 transformation for a truss)
        T4 = self.transformation_matrix()  # 4×4
        # Extract only translational DOFs from global (positions 0,1,3,4)
        u_global_4 = u_global[[0, 1, 3, 4]]  # 4×1
        u_local_4 = T4 @ u_global_4  # 4×1

        # Compute local forces (4×1)
        k_local_4 = self.local_stiffness()  # 4×4
        f_local_4 = k_local_4 @ u_local_4

        # Expand to 6×1 (insert zeros for moments)
        f_local = np.zeros(6)
        f_local[[0, 1, 3, 4]] = f_local_4
        return f_local

    def deformed_shape_points(self, global_disp, n_points=2, scale=1.0):
        """
        For a truss, just return the two displaced end points.
        """
        u_i = global_disp[self.node_i.dofs]
        u_j = global_disp[self.node_j.dofs]

        # Original coordinates
        x_i, y_i = self.node_i.x, self.node_i.y
        x_j, y_j = self.node_j.x, self.node_j.y

        # Displaced coordinates
        x_i_def = x_i + scale * u_i[0]
        y_i_def = y_i + scale * u_i[1]
        x_j_def = x_j + scale * u_j[0]
        y_j_def = y_j + scale * u_j[1]

        return [(x_i_def, y_i_def), (x_j_def, y_j_def)]
