"""
Materially and geometrically non-linear 2D beam element.
"""

import numpy as np
from fem2d.elements.element import ElementBase


class BeamElementMaterialNL(ElementBase):
    """
    Geometrically and materially non-linear 2D beam element.
    Uses a corotational formulation for geometric nonlinearity and numerical
    integration (Gauss-Lobatto) along the length for material nonlinearity.

    Attributes
    ----------
    id : int or str
        Unique identifier of the element.
    node_i : Node
        Start node.
    node_j : Node
        End node.
    section : Section
        Section definition (e.g., FiberSection, MomentCurvatureSection).
    n_int_pts : int
        Number of Gauss-Lobatto integration points along the length (3 or 5).
    """

    def __init__(self, eid, node_i, node_j, section, n_int_pts=3):
        """Initialize the BeamElementMaterialNL."""
        super().__init__(eid, node_i, node_j)
        self.section = section
        self.n_int_pts = n_int_pts

        self.L0 = None  # initial length
        self.alpha0 = None  # initial angle

        # Select Gauss-Lobatto integration points (on [-1, 1]) and weights
        if n_int_pts == 3:
            self.xi_pts = np.array([-1.0, 0.0, 1.0])
            self.weights = np.array([1.0 / 3.0, 4.0 / 3.0, 1.0 / 3.0])
        elif n_int_pts == 5:
            self.xi_pts = np.array([-1.0, -np.sqrt(3.0 / 7.0), 0.0, np.sqrt(3.0 / 7.0), 1.0])
            self.weights = np.array([0.1, 49.0 / 90.0, 32.0 / 45.0, 49.0 / 90.0, 0.1])
        else:
            raise ValueError("n_int_pts must be 3 or 5 for Gauss-Lobatto integration.")

        # Initialize history states for all integration points
        self.committed_states = [None] * self.n_int_pts
        self.trial_states = [None] * self.n_int_pts

        # State outputs (updated during iterations)
        self.f_local = None
        self.F_global = None
        self.k_local = None
        self.K_global = None

    def _update_initial_geometry(self):
        """Compute initial length and chord angle."""
        dx = self.node_j.x - self.node_i.x
        dy = self.node_j.y - self.node_i.y
        self.L0 = np.hypot(dx, dy)
        self.alpha0 = np.arctan2(dy, dx)

    def update_state(self, global_disp):
        """
        Update internal force and tangent stiffness matrices given global displacements.
        """
        if self.L0 is None:
            self._update_initial_geometry()

        # Nodal displacements (global)
        u_i = global_disp[self.node_i.dofs]  # [ux, uy, theta]
        u_j = global_disp[self.node_j.dofs]

        # Current coordinates of nodes
        x_i = self.node_i.x + u_i[0]
        y_i = self.node_i.y + u_i[1]
        x_j = self.node_j.x + u_j[0]
        y_j = self.node_j.y + u_j[1]

        # Current chord length and angle
        dx = x_j - x_i
        dy = y_j - y_i
        L = np.hypot(dx, dy)
        alpha = np.arctan2(dy, dx)

        # Rigid body rotation of the chord
        beta = alpha - self.alpha0

        # Local rotations relative to the chord
        theta1_rel = u_i[2] - beta
        theta2_rel = u_j[2] - beta

        # Axial elongation
        u_axial = L - self.L0

        # Local force and stiffness integration (corotational frame)
        f_corot = np.zeros(3)  # [N_axial, M_i, M_j]
        k_corot = np.zeros((3, 3))

        self.trial_states = []

        for p in range(self.n_int_pts):
            xi = self.xi_pts[p]
            w = self.weights[p]
            W = w * (self.L0 / 2.0)  # Map weight to [0, L0] length

            # Physical coordinate along length (x in [0, L0])
            x = 0.5 * (xi + 1.0) * self.L0

            # Strain displacement matrix (B_corot)
            # eps_axial = u_axial / L0
            # curvature = N''_v2 * theta1_rel + N''_v4 * theta2_rel
            N_pp_v2 = -4.0 / self.L0 + 6.0 * x / (self.L0 ** 2)
            N_pp_v4 = -2.0 / self.L0 + 6.0 * x / (self.L0 ** 2)

            eps_axial = u_axial / self.L0
            kappa = N_pp_v2 * theta1_rel + N_pp_v4 * theta2_rel

            # Query section response
            sec_disp = np.array([eps_axial, kappa])
            sec_force, sec_tangent, sec_state = self.section.get_response(
                sec_disp, self.committed_states[p]
            )
            self.trial_states.append(sec_state)

            # Assemble corotational force vector contribution
            # f_corot = integral( B^T * sec_force )
            f_corot[0] += W * (1.0 / self.L0) * sec_force[0]
            f_corot[1] += W * N_pp_v2 * sec_force[1]
            f_corot[2] += W * N_pp_v4 * sec_force[1]

            # Assemble corotational stiffness matrix contribution
            # k_corot = integral( B^T * sec_tangent * B )
            # Let sec_tangent = [[k11, k12], [k21, k22]]
            k11, k12 = sec_tangent[0, 0], sec_tangent[0, 1]
            k21, k22 = sec_tangent[1, 0], sec_tangent[1, 1]

            # Column 1 (d_axial)
            k_corot[0, 0] += W * (1.0 / self.L0 ** 2) * k11
            k_corot[1, 0] += W * (1.0 / self.L0) * N_pp_v2 * k21
            k_corot[2, 0] += W * (1.0 / self.L0) * N_pp_v4 * k21

            # Column 2 (theta1_rel)
            k_corot[0, 1] += W * (1.0 / self.L0) * N_pp_v2 * k12
            k_corot[1, 1] += W * N_pp_v2 * N_pp_v2 * k22
            k_corot[2, 1] += W * N_pp_v4 * N_pp_v2 * k22

            # Column 3 (theta2_rel)
            k_corot[0, 2] += W * (1.0 / self.L0) * N_pp_v4 * k12
            k_corot[1, 2] += W * N_pp_v2 * N_pp_v4 * k22
            k_corot[2, 2] += W * N_pp_v4 * N_pp_v4 * k22

        # Transform from corotational coordinates (size 3) to local DOFs (size 6)
        T_corot = np.array(
            [
                [-1.0, 0.0, 0.0, 1.0, 0.0, 0.0],
                [0.0, 1.0 / self.L0, 1.0, 0.0, -1.0 / self.L0, 0.0],
                [0.0, 1.0 / self.L0, 0.0, 0.0, -1.0 / self.L0, 1.0],
            ]
        )

        f_local = T_corot.T @ f_corot
        k_local = T_corot.T @ k_corot @ T_corot

        # Add geometric stiffness matrix k_geo (similar to BeamElementNL)
        N_axial = f_corot[0]  # tension force
        k_geo = np.zeros((6, 6))
        factor = N_axial / (30.0 * self.L0)

        k_geo[1, 1] = 36.0
        k_geo[1, 2] = 3.0 * self.L0
        k_geo[1, 4] = -36.0
        k_geo[1, 5] = 3.0 * self.L0

        k_geo[2, 1] = 3.0 * self.L0
        k_geo[2, 2] = 4.0 * self.L0 ** 2
        k_geo[2, 4] = -3.0 * self.L0
        k_geo[2, 5] = -(self.L0 ** 2)

        k_geo[4, 1] = -36.0
        k_geo[4, 2] = -3.0 * self.L0
        k_geo[4, 4] = 36.0
        k_geo[4, 5] = -3.0 * self.L0

        k_geo[5, 1] = 3.0 * self.L0
        k_geo[5, 2] = -(self.L0 ** 2)
        k_geo[5, 4] = -3.0 * self.L0
        k_geo[5, 5] = 4.0 * self.L0 ** 2

        k_geo *= factor
        k_local += k_geo

        # Transform to global coordinates
        c = np.cos(alpha)
        s = np.sin(alpha)
        R = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
        T = np.zeros((6, 6))
        T[:3, :3] = R
        T[3:, 3:] = R

        self.f_local = f_local
        self.F_global = T @ f_local
        self.k_local = k_local
        self.K_global = T @ k_local @ T.T

    def commit_state(self):
        """Commit the history/state variables of all integration points."""
        self.committed_states = self.trial_states

    def get_tangent_stiffness(self):
        """Return global tangent stiffness matrix."""
        return self.K_global

    def get_internal_forces(self):
        """Return global internal forces vector."""
        return self.F_global

    def get_local_forces(self):
        """Return local internal forces vector."""
        return self.f_local
