"""
Non-linear truss element module implementing a geometrically non-linear 2D truss element.
"""

import numpy as np
from fem2d.elements.element import ElementBase


class TrussElementNL(ElementBase):
    """
    Geometrically nonlinear truss element using corotational formulation.
    Nodes have 2 translational DOFs: [ux, uy], with rotation index mapped to zero.

    Attributes
    ----------
    material : ElasticMaterial
        Material definition for the element.
    area : float
        Cross-sectional area.
    L0 : float
        Initial (undeformed) length.
    Q : float
        Axial force (compression positive, tension negative).
    cx : float
        Current direction cosine along the global x-axis.
    cy : float
        Current direction cosine along the global y-axis.
    L : float
        Current deformed length.
    ke : numpy.ndarray or None
        Elastic stiffness matrix in global coordinates.
    kg : numpy.ndarray or None
        Geometric stiffness matrix in global coordinates.
    k_t : numpy.ndarray or None
        Tangent stiffness matrix in global coordinates.
    F_global : numpy.ndarray or None
        Internal force vector in global coordinates.
    """

    def __init__(self, eid, node_i, node_j, material, area):
        """
        Initialize a TrussElementNL.

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
        """
        super().__init__(eid, node_i, node_j)
        self.material = material
        self.area = area
        # Store initial length (undeformed)
        self.L0 = self.length  # computed by _update_geometry() in base class
        # State variables
        self.Q = 0.0  # axial force (compression positive, as in text)
        self.cx = self.cos  # current direction cosines (initial)
        self.cy = self.sin
        self.L = self.L0  # current length
        self.ke = None  # elastic stiffness matrix (global)
        self.kg = None  # geometric stiffness matrix (global)
        self.k_t = None  # tangent stiffness (global)

    def update_state(self, global_disp):
        """
        Given global displacement vector of the whole structure,
        extract this element's end displacements and update internal state (internal force, tangent stiffness).

        Parameters
        ----------
        global_disp : numpy.ndarray
            Full global displacement vector of the structure.
        """
        # End displacements in global coordinates (2 per node)
        u_i = global_disp[self.node_i.dofs]  # [ux_i, uy_i]
        u_j = global_disp[self.node_j.dofs]  # [ux_j, uy_j]

        # Current coordinates of nodes
        x_i = self.node_i.x + u_i[0]
        y_i = self.node_i.y + u_i[1]
        x_j = self.node_j.x + u_j[0]
        y_j = self.node_j.y + u_j[1]

        # Current length and direction cosines (Eqs. 10.20–10.22)
        dx = x_j - x_i
        dy = y_j - y_i
        self.L = np.hypot(dx, dy)
        if self.L == 0:
            self.cx = 1.0
            self.cy = 0.0
        else:
            self.cx = dx / self.L
            self.cy = dy / self.L

        # Axial deformation (shortening positive) – Eq. (10.15)
        u = self.L0 - self.L

        # Axial force (compression positive) – Eq. (10.18)
        self.Q = (self.material.E * self.area / self.L0) * u

        # Global internal force vector (size 4) – Eq. (10.24)
        F = np.array(
            [self.cx * self.Q, self.cy * self.Q, -self.cx * self.Q, -self.cy * self.Q]
        )

        # Expand to 6x1 for assembly (mapping to [ux_i, uy_i, 0, ux_j, uy_j, 0])
        self.F_global = np.zeros(6)
        self.F_global[[0, 1, 3, 4]] = F

        # Tangent stiffness matrix (global, 4x4) – Eq. (10.31)
        EA_L = self.material.E * self.area / self.L0
        Q_L = self.Q / self.L
        c = self.cx
        s = self.cy

        # Elastic part (EA/L * T^T T)
        ke_local = EA_L * np.array(
            [
                [c * c, c * s, -c * c, -c * s],
                [c * s, s * s, -c * s, -s * s],
                [-c * c, -c * s, c * c, c * s],
                [-c * s, -s * s, c * s, s * s],
            ]
        )

        # Geometric part (Q/L * g) with g from Eq. (10.33)
        kg_local = Q_L * np.array(
            [
                [-s * s, c * s, s * s, -c * s],
                [c * s, -c * c, -c * s, c * c],
                [s * s, -c * s, -s * s, c * s],
                [-c * s, c * c, c * s, -c * c],
            ]
        )

        k_t_local = ke_local + kg_local

        # Expand to 6x6 (insert zeros for rotational DOFs)
        self.k_t = np.zeros((6, 6))
        idx = [0, 1, 3, 4]
        for i, ii in enumerate(idx):
            for j, jj in enumerate(idx):
                self.k_t[ii, jj] = k_t_local[i, j]

    def get_tangent_stiffness(self):
        """
        Return the global tangent stiffness matrix of the element.

        Returns
        -------
        numpy.ndarray
            6x6 tangent stiffness matrix in global coordinates.
        """
        return self.k_t

    def get_internal_forces(self):
        """
        Return the global internal force vector of the element.

        Returns
        -------
        numpy.ndarray
            6x1 global internal force vector.
        """
        return self.F_global

    def get_local_forces(self):
        """
        Return local end forces as a 6‑component array:
        [Fx_i, Fy_i, Mz_i, Fx_j, Fy_j, Mz_j]
        Axial force Q is positive in compression.

        Returns
        -------
        numpy.ndarray
            6x1 force vector in local coordinates.
        """
        return np.array([self.Q, 0.0, 0.0, -self.Q, 0.0, 0.0])
