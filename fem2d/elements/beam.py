"""
Beam element module defining 2D elastic Euler-Bernoulli beam elements with mass options.
"""

from fem2d.elements.element import ElementBase
import numpy as np


class BeamElement(ElementBase):
    """
    Elastic 2D Euler-Bernoulli beam element including axial, shear, and bending stiffness.

    Attributes
    ----------
    material : ElasticMaterial
        Material definition for the element.
    area : float
        Cross-sectional area.
    inertia : float
        Moment of inertia of the cross-section.
    extra_mass : float
        Additional distributed mass per unit length (e.g. non-structural mass).
    """

    def __init__(self, eid, node_i, node_j, material, area, inertia, extra_mass=0.0):
        """
        Initialize a BeamElement.

        Parameters
        ----------
        eid : int or str
            Unique identifier of the element.
        node_i : Node
            Start node.
        node_j : Node
            End node.
        material : ElasticMaterial
            Material definition.
        area : float
            Cross-sectional area.
        inertia : float
            Moment of inertia.
        extra_mass : float, optional
            Additional distributed mass per unit length. Defaults to 0.0.
        """
        super().__init__(eid, node_i, node_j)
        self.material = material
        self.area = area
        self.inertia = inertia  # Moment of inertia
        self.extra_mass = extra_mass  # kg/m, additional distributed mass

    def local_stiffness(self):
        """
        Return the 6x6 element stiffness matrix in local coordinates.

        Returns
        -------
        numpy.ndarray
            Stiffness matrix in local coordinate system.
        """
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

    # def equivalent_nodal_loads(self):
    #     """
    #     Return local equivalent nodal loads due to distributed element loads.

    #     Returns
    #     -------
    #     numpy.ndarray
    #         Equivalent nodal load vector (6x1).
    #     """
    #     # If a distributed load is stored, compute fixed‑end forces.
    #     # For example, uniform load w in local y‑direction:
    #     if hasattr(self, "w"):
    #         w = self.w
    #         L = self.length
    #         return np.array([0, w * L / 2, w * L**2 / 12, 0, w * L / 2, -w * L**2 / 12])
    #     return np.zeros(6)

    def get_local_forces(self):
        """
        Compute and return the element internal forces in local coordinates.

        Returns
        -------
        numpy.ndarray
            6-component vector [Fx_i, Fy_i, Mz_i, Fx_j, Fy_j, Mz_j] in local coordinates.
        """
        u_i = self.structure.disp[self.node_i.dofs]
        u_j = self.structure.disp[self.node_j.dofs]
        u_global = np.concatenate([u_i, u_j])
        T = self.transformation_matrix()  # 6×6
        u_local = T @ u_global
        f_local = self.local_stiffness() @ u_local
        # Subtract equivalent nodal loads if any (e.g., from distributed loads)
        if hasattr(self, "eq_load") and self.eq_load is not None:
            eq_load_local = T @ self.eq_load
            f_local -= eq_load_local
        return f_local

    def deformed_shape_points(self, global_disp, n_points=20, scale=1.0):
        """
        Return a list of (x, y) points in global coordinates representing
        the deformed shape of the beam, scaled by `scale`.

        Parameters
        ----------
        global_disp : numpy.ndarray
            Full global displacement vector of the structure.
        n_points : int, optional
            Number of points to generate along the element. Defaults to 20.
        scale : float, optional
            Displacement magnification factor. Defaults to 1.0.

        Returns
        -------
        list of tuple of float
            List of (x, y) coordinates representing the deformed shape.
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

    def mass_matrix(self):
        """
        Return the 6x6 global mass matrix (including rotational inertia and extra mass).

        Returns
        -------
        numpy.ndarray
            6x6 global mass matrix.
        """
        L = self.length
        A = self.area
        I = self.inertia
        rho = self.material.rho
        extra = self.extra_mass

        # Effective density for translational part: material + extra
        rho_eff = rho + (extra / A) if A != 0 else rho

        # ---- Translational mass (m1) ----
        factor_trans = rho_eff * A * L / 420.0
        m1 = np.zeros((6, 6))
        m1[0, 0] = m1[3, 3] = 140.0 * factor_trans
        m1[0, 3] = m1[3, 0] = 70.0 * factor_trans

        m1[1, 1] = 156.0 * factor_trans
        m1[1, 2] = 22.0 * L * factor_trans
        m1[1, 4] = 54.0 * factor_trans
        m1[1, 5] = -13.0 * L * factor_trans

        m1[2, 1] = 22.0 * L * factor_trans
        m1[2, 2] = 4.0 * L**2 * factor_trans
        m1[2, 4] = 13.0 * L * factor_trans
        m1[2, 5] = -3.0 * L**2 * factor_trans

        m1[4, 1] = 54.0 * factor_trans
        m1[4, 2] = 13.0 * L * factor_trans
        m1[4, 4] = 156.0 * factor_trans
        m1[4, 5] = -22.0 * L * factor_trans

        m1[5, 1] = -13.0 * L * factor_trans
        m1[5, 2] = -3.0 * L**2 * factor_trans
        m1[5, 4] = -22.0 * L * factor_trans
        m1[5, 5] = 4.0 * L**2 * factor_trans

        # ---- Rotational mass (m2) ----
        factor_rot = rho * I / (30.0 * L)
        m2 = np.zeros((6, 6))
        m2[1, 1] = 36.0 * factor_rot
        m2[1, 2] = 3.0 * L * factor_rot
        m2[1, 4] = -36.0 * factor_rot
        m2[1, 5] = 3.0 * L * factor_rot

        m2[2, 1] = 3.0 * L * factor_rot
        m2[2, 2] = 4.0 * L**2 * factor_rot
        m2[2, 4] = -3.0 * L * factor_rot
        m2[2, 5] = -(L**2) * factor_rot

        m2[4, 1] = -36.0 * factor_rot
        m2[4, 2] = -3.0 * L * factor_rot
        m2[4, 4] = 36.0 * factor_rot
        m2[4, 5] = -3.0 * L * factor_rot

        m2[5, 1] = 3.0 * L * factor_rot
        m2[5, 2] = -(L**2) * factor_rot
        m2[5, 4] = -3.0 * L * factor_rot
        m2[5, 5] = 4.0 * L**2 * factor_rot

        # total local mass matrix
        m_local = m1 + m2

        # Transform to global
        T = self.transformation_matrix()
        m_global = T.T @ m_local @ T
        return m_global

    def axial_force(self):
        """Return the axial force in the beam element in local coordinates.

        Positive value corresponds to compression for buckling sign convention.
        """
        f_local = self.get_local_forces()
        # f_local[0] is axial force at node i (positive in tension for stiffness)
        # invert sign so compression is positive
        return -f_local[0]

    def geometric_stiffness(self, axial_force):
        """Return the 6x6 geometric stiffness matrix (local->global) for the beam.

        Uses the standard Euler-Bernoulli beam geometric stiffness for axial
        compression P (positive in compression):

            Kg_local = (P / (30*L)) * [
                [0,   0,    0,    0,    0,    0],
                [0,  36,   3L,    0,  -36,   3L],
                [0, 3L,  4L^2,   0,  -3L,  -L^2],
                [0,   0,    0,    0,    0,    0],
                [0, -36,  -3L,   0,   36,  -3L],
                [0,  3L,  -L^2,  0,  -3L,  4L^2]
            ]

        The returned matrix is in global coordinates.
        """
        P = axial_force
        L = self.length
        if L == 0:
            return np.zeros((6, 6))

        factor = P / (30.0 * L)
        kg_local = factor * np.array(
            [
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                [0.0, 36.0, 3.0 * L, 0.0, -36.0, 3.0 * L],
                [0.0, 3.0 * L, 4.0 * L * L, 0.0, -3.0 * L, -1.0 * L * L],
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                [0.0, -36.0, -3.0 * L, 0.0, 36.0, -3.0 * L],
                [0.0, 3.0 * L, -1.0 * L * L, 0.0, -3.0 * L, 4.0 * L * L],
            ]
        )

        T = self.transformation_matrix()
        kg_global = T.T @ kg_local @ T
        return kg_global
