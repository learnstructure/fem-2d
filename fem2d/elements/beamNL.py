import numpy as np


class BeamElementNL:
    """
    Geometrically nonlinear 2D beam element using corotational formulation
    (Euler-Bernoulli). Follows the interface required by NewtonRaphsonSolver.
    """

    def __init__(self, eid, node_i, node_j, material, area, inertia):
        self.id = eid
        self.node_i = node_i
        self.node_j = node_j
        self.material = material
        self.area = area
        self.inertia = inertia
        self.L0 = None  # initial length
        self.alpha0 = None  # initial chord angle (rad)
        self.structure = None  # set when added to a Structure

        # state variables (updated during iterations)
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
        Compute current internal forces and tangent stiffness
        given the global displacement vector of the whole structure.
        """
        if self.L0 is None:
            self._update_initial_geometry()

        # nodal displacements (global)
        u_i = global_disp[self.node_i.dofs]  # [ux, uy, theta]
        u_j = global_disp[self.node_j.dofs]

        # current nodal coordinates
        x_i = self.node_i.x + u_i[0]
        y_i = self.node_i.y + u_i[1]
        x_j = self.node_j.x + u_j[0]
        y_j = self.node_j.y + u_j[1]

        # current chord length and angle
        dx = x_j - x_i
        dy = y_j - y_i
        L = np.hypot(dx, dy)
        alpha = np.arctan2(dy, dx)

        # rigid body rotation
        beta = alpha - self.alpha0

        # local rotations (relative to chord)
        theta1_bar = u_i[2] - beta
        theta2_bar = u_j[2] - beta

        # axial deformation (elongation)
        u2 = L - self.L0  # positive = tension

        # local displacement vector (order: u1,w1,θ1, u2,w2,θ2)
        d_local = np.array([0.0, 0.0, theta1_bar, u2, 0.0, theta2_bar])

        # material properties
        E = self.material.E
        A = self.area
        I = self.inertia
        L0 = self.L0

        # axial force (tension positive)
        EA_L = E * A / L0
        f_x2 = EA_L * u2

        # ---------- elastic stiffness matrix k1 (Euler-Bernoulli) ----------
        k1 = np.zeros((6, 6))

        # axial
        k1[0, 0] = k1[3, 3] = EA_L
        k1[0, 3] = k1[3, 0] = -EA_L

        # bending
        EI_L3 = E * I / L0**3
        EI_L2 = E * I / L0**2
        EI_L = E * I / L0

        k1[1, 1] = 12 * EI_L3
        k1[1, 2] = 6 * EI_L2
        k1[1, 4] = -12 * EI_L3
        k1[1, 5] = 6 * EI_L2

        k1[2, 1] = 6 * EI_L2
        k1[2, 2] = 4 * EI_L
        k1[2, 4] = -6 * EI_L2
        k1[2, 5] = 2 * EI_L

        k1[4, 1] = -12 * EI_L3
        k1[4, 2] = -6 * EI_L2
        k1[4, 4] = 12 * EI_L3
        k1[4, 5] = -6 * EI_L2

        k1[5, 1] = 6 * EI_L2
        k1[5, 2] = 2 * EI_L
        k1[5, 4] = -6 * EI_L2
        k1[5, 5] = 4 * EI_L

        # ---------- geometric stiffness matrix k2 (from paper) ----------
        # only rows/cols 1,2,4,5 (0‑based indices) are non‑zero
        factor = f_x2 / (30.0 * L0)
        k2 = np.zeros((6, 6))

        k2[1, 1] = 36
        k2[1, 2] = 3 * L0
        k2[1, 4] = -36
        k2[1, 5] = 3 * L0

        k2[2, 1] = 3 * L0
        k2[2, 2] = 4 * L0**2
        k2[2, 4] = -3 * L0
        k2[2, 5] = -(L0**2)

        k2[4, 1] = -36
        k2[4, 2] = -3 * L0
        k2[4, 4] = 36
        k2[4, 5] = -3 * L0

        k2[5, 1] = 3 * L0
        k2[5, 2] = -(L0**2)
        k2[5, 4] = -3 * L0
        k2[5, 5] = 4 * L0**2

        k2 *= factor

        # total local stiffness
        k_local = k1 + k2

        # local internal forces
        f_local = k_local @ d_local

        # ---------- transformation to global ----------
        c = np.cos(alpha)
        s = np.sin(alpha)
        R = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])
        T = np.zeros((6, 6))
        T[:3, :3] = R
        T[3:, 3:] = R

        # global internal forces
        F_global = T @ f_local

        # global tangent stiffness
        K_global = T @ k_local @ T.T

        # store for later queries
        self.f_local = f_local
        self.F_global = F_global
        self.k_local = k_local
        self.K_global = K_global

    def get_tangent_stiffness(self):
        return self.K_global

    def get_internal_forces(self):
        return self.F_global

    def get_local_forces(self):
        """Return local end forces (6 components) for post‑processing."""
        return self.f_local
