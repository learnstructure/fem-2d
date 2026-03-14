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

    def get_local_forces(self):
        """Return local end forces as a 6‑component array."""
        u_i = self.structure.disp[self.node_i.dofs]
        u_j = self.structure.disp[self.node_j.dofs]
        u_global = np.concatenate([u_i, u_j])
        T = self.transformation_matrix()  # 6×6
        u_local = T @ u_global
        f_local = self.local_stiffness() @ u_local
        # Subtract equivalent nodal loads if any (e.g., from distributed loads)
        if hasattr(self, "eq_load") and self.eq_load is not None:
            f_local -= self.eq_load
        return f_local

    def deformed_shape_points(self, global_disp, n_points=20, scale=1.0):
        """
        Returns a list of (x, y) points in global coordinates representing
        the deformed shape of the beam, scaled by `scale`.
        global_disp: full global displacement vector of the structure.
        """
        # Extract nodal displacements (global)
        u_i = global_disp[self.node_i.dofs]  # [ux_i, uy_i, rz_i]
        u_j = global_disp[self.node_j.dofs]  # [ux_j, uy_j, rz_j]

        # Build 6‑vector of global displacements
        u_global = np.concatenate([u_i, u_j])

        # Transform to local displacements (local = T @ global)
        T = self.transformation_matrix()
        u_local = T @ u_global

        # Local displacements: [u1, v1, θ1, u2, v2, θ2]
        u1, v1, theta1, u2, v2, theta2 = u_local

        L = self.length
        c, s = self.cos, self.sin

        # Generate points along the element
        points = []
        for xi in np.linspace(0, L, n_points):
            # Local coordinate
            x_local = xi
            xi_norm = x_local / L  # ξ

            # Axial displacement (linear interpolation)
            u_axial = u1 * (1 - xi_norm) + u2 * xi_norm

            # Transverse displacement (cubic Hermitian)
            # Shape functions
            N1 = 1 - 3 * xi_norm**2 + 2 * xi_norm**3
            N2 = (xi_norm - 2 * xi_norm**2 + xi_norm**3) * L
            N3 = 3 * xi_norm**2 - 2 * xi_norm**3
            N4 = (-(xi_norm**2) + xi_norm**3) * L

            v_trans = v1 * N1 + theta1 * N2 + v2 * N3 + theta2 * N4

            # Local displacement vector at this point (only translations)
            disp_local = np.array([u_axial, v_trans])

            # Transform back to global
            R = T[:2, :2]
            disp_global = R.T @ disp_local

            # Original global coordinates of this point (linear interpolation)
            x_orig = self.node_i.x + xi_norm * (self.node_j.x - self.node_i.x)
            y_orig = self.node_i.y + xi_norm * (self.node_j.y - self.node_i.y)

            # Deformed coordinates (apply scaled displacement)
            x_def = x_orig + scale * disp_global[0]
            y_def = y_orig + scale * disp_global[1]

            points.append((x_def, y_def))

        return points
