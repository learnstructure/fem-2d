"""
Truss element module defining 2D elastic truss (bar) elements with axial stiffness.
"""

from fem2d.elements.element import ElementBase
import numpy as np


class TrussElement(ElementBase):
    """
    Elastic 2D truss (pin-jointed bar) element with only axial stiffness.

    Attributes
    ----------
    material : ElasticMaterial or None
        Material definition for the element.
    area : float or None
        Cross-sectional area.
    extra_mass : float
        Additional distributed mass per unit length.
    """

    def __init__(
        self,
        eid,
        node_i,
        node_j,
        material,
        section,
        extra_mass=0.0,
        geometric_nonlinear=False,
    ):
        """
        Initialize a TrussElement.

        Parameters
        ----------
        eid : int or str
            Unique identifier of the element.
        node_i : Node
            Start node.
        node_j : Node
            End node.
        material : ElasticMaterial or UniaxialMaterial
            Material definition.
        section : Section
            Section definition (containing area A).
        extra_mass : float, optional
            Additional distributed mass per unit length. Defaults to 0.0.
        geometric_nonlinear : bool, optional
            Whether to use corotational geometric non-linear formulation. Defaults to False.
        """
        super().__init__(eid, node_i, node_j)

        if hasattr(section, "A"):
            self.section = section
            self.area = section.A
            self.material = material
        elif hasattr(material, "A"):
            self.section = material
            self.area = material.A
            self.material = getattr(material, "material", None)
            if isinstance(section, (int, float)):
                extra_mass = float(section)
        elif isinstance(section, (int, float)):
            from fem2d.sections.section import Section
            area_val = float(section)
            self.section = Section(area_val)
            self.area = area_val
            self.material = material
        else:
            self.material = material
            self.section = section
            self.area = getattr(section, "A", None)

        self.extra_mass = extra_mass
        self.geometric_nonlinear = geometric_nonlinear

        # Initial length and state variables
        self.L0 = self.length
        self.Q = 0.0
        self.cx = self.cos
        self.cy = self.sin
        self.L = self.L0
        self.k_t = None
        self.F_global = None

    def local_stiffness(self):
        """
        Return the 4x4 element stiffness matrix in local coordinates.

        Returns
        -------
        numpy.ndarray
            4x4 stiffness matrix.
        """
        E = self.material.E if self.material is not None else 1.0
        k = self.area * E / self.length
        return np.array([[k, 0, -k, 0], [0, 0, 0, 0], [-k, 0, k, 0], [0, 0, 0, 0]])

    def global_stiffness(self):
        """
        Assemble and return the 6x6 global stiffness matrix.
        Expands the 4x4 matrix by inserting zeros for rotational DOFs.

        Returns
        -------
        numpy.ndarray
            6x6 global stiffness matrix.
        """
        k_local = self.local_stiffness()
        T4 = self.transformation_matrix(dof_per_node=2)
        k_global_4 = T4.T @ k_local @ T4
        k_global = np.zeros((6, 6))
        dof_map = [0, 1, 3, 4]
        for i, ii in enumerate(dof_map):
            for j, jj in enumerate(dof_map):
                k_global[ii, jj] = k_global_4[i, j]
        return k_global

    def transformation_matrix(self, dof_per_node=2):
        """
        Return the 4x4 transformation matrix from local to global coordinates.

        Returns
        -------
        numpy.ndarray
            4x4 transformation matrix.
        """
        c = self.cos
        s = self.sin
        return np.array([[c, s, 0, 0], [-s, c, 0, 0], [0, 0, c, s], [0, 0, -s, c]])

    def update_state(self, global_disp):
        """
        Update internal force and stiffness matrix given global displacement vector.

        Parameters
        ----------
        global_disp : numpy.ndarray
            Full global displacement vector of the structure.
        """
        if self.geometric_nonlinear:
            # Corotational geometrically non-linear formulation
            u_i = global_disp[self.node_i.dofs]  # [ux, uy, rz]
            u_j = global_disp[self.node_j.dofs]

            x_i = self.node_i.x + u_i[0]
            y_i = self.node_i.y + u_i[1]
            x_j = self.node_j.x + u_j[0]
            y_j = self.node_j.y + u_j[1]

            dx = x_j - x_i
            dy = y_j - y_i
            self.L = np.hypot(dx, dy)
            if self.L == 0:
                self.cx = 1.0
                self.cy = 0.0
            else:
                self.cx = dx / self.L
                self.cy = dy / self.L

            u = self.L0 - self.L
            E = self.material.E if self.material is not None else 1.0
            self.Q = (E * self.area / self.L0) * u

            F = np.array(
                [self.cx * self.Q, self.cy * self.Q, -self.cx * self.Q, -self.cy * self.Q]
            )

            self.F_global = np.zeros(6)
            self.F_global[[0, 1, 3, 4]] = F

            EA_L = E * self.area / self.L0
            Q_L = self.Q / self.L
            c = self.cx
            s = self.cy

            ke_local = EA_L * np.array(
                [
                    [c * c, c * s, -c * c, -c * s],
                    [c * s, s * s, -c * s, -s * s],
                    [-c * c, -c * s, c * c, c * s],
                    [-c * s, -s * s, c * s, s * s],
                ]
            )

            kg_local = Q_L * np.array(
                [
                    [-s * s, c * s, s * s, -c * s],
                    [c * s, -c * c, -c * s, c * c],
                    [s * s, -c * s, -s * s, c * s],
                    [-c * s, c * c, c * s, -c * c],
                ]
            )

            k_t_local = ke_local + kg_local

            self.k_t = np.zeros((6, 6))
            idx = [0, 1, 3, 4]
            for i, ii in enumerate(idx):
                for j, jj in enumerate(idx):
                    self.k_t[ii, jj] = k_t_local[i, j]
        else:
            # Linear state update
            u_i = global_disp[self.node_i.dofs]
            u_j = global_disp[self.node_j.dofs]
            u_global = np.concatenate([u_i, u_j])
            K_g = self.global_stiffness()
            self.F_global = K_g @ u_global
            self.k_t = K_g

    def get_tangent_stiffness(self):
        """Return global tangent stiffness matrix of the element."""
        if hasattr(self, "k_t") and self.k_t is not None:
            return self.k_t
        return self.global_stiffness()

    def get_internal_forces(self):
        """Return global internal forces vector of the element."""
        if hasattr(self, "F_global") and self.F_global is not None:
            return self.F_global
        return np.zeros(6)

    def get_local_forces(self):
        """Compute and return local end forces."""
        if self.geometric_nonlinear and hasattr(self, "Q"):
            return np.array([self.Q, 0.0, 0.0, -self.Q, 0.0, 0.0])

        if hasattr(self, "structure") and self.structure is not None and self.structure.disp is not None:
            u_i = self.structure.disp[self.node_i.dofs]
            u_j = self.structure.disp[self.node_j.dofs]
            u_global = np.concatenate([u_i, u_j])
            T4 = self.transformation_matrix()
            u_global_4 = u_global[[0, 1, 3, 4]]
            u_local_4 = T4 @ u_global_4
            f_local_4 = self.local_stiffness() @ u_local_4
            f_local = np.zeros(6)
            f_local[[0, 1, 3, 4]] = f_local_4
            return f_local
        return np.zeros(6)

    def axial_force(self):
        """Return the axial force in the truss element in local coordinates.

        Positive values correspond to compression.
        """
        local_forces = self.get_local_forces()
        return local_forces[0]      #compression positive for buckling analysis

    def geometric_stiffness(self, axial_force):
        """Return the 6x6 geometric stiffness matrix for the truss element.

        Parameters
        ----------
        axial_force : float
            Axial force in the element (positive in compression).

        Returns
        -------
        numpy.ndarray
            6x6 geometric stiffness matrix in global coordinates.
        """
        N = axial_force
        if self.length == 0:
            return np.zeros((6, 6))

        kg_local = (N/self.length) * np.array([
    [0, 0, 0, 0],
    [0, 1, 0,-1],
    [0, 0, 0, 0],
    [0,-1, 0, 1]
])
        # kg_local = (N / self.length) * np.array([[1, 0, -1, 0], [0, 0, 0, 0], [-1, 0, 1, 0], [0, 0, 0, 0]])

        T4 = self.transformation_matrix(dof_per_node=2)
        kg_global_4 = T4.T @ kg_local @ T4
        kg_global = np.zeros((6, 6))
        dof_map = [0, 1, 3, 4]
        for i, ii in enumerate(dof_map):
            for j, jj in enumerate(dof_map):
                kg_global[ii, jj] = kg_global_4[i, j]
        return kg_global

    def deformed_shape_points(self, global_disp, n_points=2, scale=1.0):
        """
        Return end point coordinates representing the deformed shape of the truss.

        Parameters
        ----------
        global_disp : numpy.ndarray
            Full global displacement vector.
        n_points : int, optional
            Number of points along the element. Defaults to 2.
        scale : float, optional
            Displacement magnification factor. Defaults to 1.0.

        Returns
        -------
        list of tuple of float
            End point coordinates (start and end).
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

    def mass_matrix(self):
        """
        Return the 6x6 global mass matrix.
        Lumps the element and extra mass equally to translational DOFs at ends.

        Returns
        -------
        numpy.ndarray
            6x6 global mass matrix.
        """
        L = self.length
        A = self.area
        rho = self.material.rho
        extra = self.extra_mass

        total_mass = (rho * A + extra) * L
        # Lumped at nodes: half at each node, translational only
        m_local = np.zeros((4, 4))
        m_local[0, 0] = m_local[1, 1] = m_local[2, 2] = m_local[3, 3] = total_mass / 2.0

        # Expand to 6x6 (zero for rotational DOFs)
        m_global = np.zeros((6, 6))
        idx = [0, 1, 3, 4]  # mapping from local (4) to global (6) DOFs
        for i, ii in enumerate(idx):
            for j, jj in enumerate(idx):
                m_global[ii, jj] = m_local[i, j]

        # Transform to global orientation if needed (since truss can be rotated)
        T4 = self.transformation_matrix()
        T_full = np.zeros((6, 6))
        T_full[0:2, 0:2] = T4[0:2, 0:2]
        T_full[3:5, 3:5] = T4[2:4, 2:4]
        m_global = T_full.T @ m_global @ T_full
        return m_global
