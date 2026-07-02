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

    def __init__(self, eid, node_i, node_j, material, area, extra_mass=0.0):
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
        material : ElasticMaterial or None
            Material definition.
        area : float or None
            Cross-sectional area.
        extra_mass : float, optional
            Additional distributed mass per unit length. Defaults to 0.0.
        """
        super().__init__(eid, node_i, node_j)
        self.material = material
        self.area = area
        self.extra_mass = extra_mass  # kg/m additional distributed mass

    def local_stiffness(self):
        """
        Return the 4x4 element stiffness matrix in local coordinates.

        Returns
        -------
        numpy.ndarray
            4x4 stiffness matrix.
        """
        k = self.area * self.material.E / self.length
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
        """
        Return the transformation matrix from local to global coordinates.

        Parameters
        ----------
        dof_per_node : int, optional
            Number of degrees of freedom per node. Defaults to 2.

        Returns
        -------
        numpy.ndarray
            4x4 transformation matrix.
        """
        c = self.cos
        s = self.sin
        return np.array([[c, s, 0, 0], [-s, c, 0, 0], [0, 0, c, s], [0, 0, -s, c]])

    def get_local_forces(self):
        """
        Compute and return the local member end forces.

        Returns
        -------
        numpy.ndarray
            6-component vector [Fx_i, Fy_i, Mz_i, Fx_j, Fy_j, Mz_j] in local coordinates.
        """
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

    def axial_force(self):
        """Return the axial force in the truss element in local coordinates.

        Positive values correspond to compression.
        """
        local_forces = self.get_local_forces()
        # get_local_forces() returns local axial force at node i in the first entry.
        # In a standard truss sign convention, positive value indicates tension.
        # For geometric stiffness we need compression positive, so invert the sign.
        return -local_forces[0]

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
