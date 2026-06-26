"""
ElementBase module defining the base class for all structural finite elements.
"""

import numpy as np


class ElementBase:
    """
    Base class for structural elements in a 2D finite element model.

    Attributes
    ----------
    id : int or str
        Unique identifier of the element.
    node_i : Node
        Start node of the element.
    node_j : Node
        End node of the element.
    structure : Structure or None
        The parent structure containing this element.
    length : float
        Length of the element.
    cos : float
        Cosine of the element orientation angle with respect to global x-axis.
    sin : float
        Sine of the element orientation angle with respect to global x-axis.
    """

    def __init__(self, eid, node_i, node_j):
        """
        Initialize an ElementBase object.

        Parameters
        ----------
        eid : int or str
            Unique identifier of the element.
        node_i : Node
            Start node of the element.
        node_j : Node
            End node of the element.
        """
        self.id = eid
        self.node_i = node_i
        self.node_j = node_j
        self.structure = None
        self._update_geometry()

    def _update_geometry(self):
        """Compute length, sine, cosine from node coordinates."""
        dx = self.node_j.x - self.node_i.x
        dy = self.node_j.y - self.node_i.y
        self.length = np.hypot(dx, dy)
        self.cos = dx / self.length
        self.sin = dy / self.length

    def transformation_matrix(self, dof_per_node=3):
        """
        Return the transformation matrix from local to global coordinates.

        Parameters
        ----------
        dof_per_node : int, optional
            Number of degrees of freedom per node. Defaults to 3.

        Returns
        -------
        numpy.ndarray
            Transformation matrix of size (2*dof_per_node, 2*dof_per_node).
        """
        c = self.cos
        s = self.sin
        T = np.array(
            [
                [c, s, 0, 0, 0, 0],
                [-s, c, 0, 0, 0, 0],
                [0, 0, 1, 0, 0, 0],
                [0, 0, 0, c, s, 0],
                [0, 0, 0, -s, c, 0],
                [0, 0, 0, 0, 0, 1],
            ]
        )
        # If the element uses fewer DOFs (e.g., truss), subclasses can slice.
        return T

    def local_stiffness(self):
        """
        Return element stiffness matrix in local coordinates.

        Returns
        -------
        numpy.ndarray
            Local stiffness matrix.

        Raises
        ------
        NotImplementedError
            Must be implemented by subclasses.
        """
        raise NotImplementedError

    def global_stiffness(self):
        """
        Assemble global stiffness matrix (expanded to global system size).

        Returns
        -------
        numpy.ndarray
            Global stiffness matrix (typically 6x6).
        """
        k_local = self.local_stiffness()
        T = self.transformation_matrix()
        # For elements with fewer than 6 DOF, subclasses must expand k_local and T.
        return T.T @ k_local @ T

    def equivalent_nodal_loads(self):
        """
        Return equivalent nodal loads due to element loads in global coordinates.

        Returns
        -------
        numpy.ndarray
            Global equivalent load vector (size 6).
        """
        return np.zeros(6)  # override if needed

    def internal_forces(self, displacements):
        """
        Given nodal displacements (global, 6‑vector), return local end forces.

        Parameters
        ----------
        displacements : numpy.ndarray
            Global displacement vector for the element nodes.

        Returns
        -------
        numpy.ndarray
            Local end forces vector.

        Raises
        ------
        NotImplementedError
            Must be implemented by subclasses.
        """
        # Used for post‑processing and non‑linear iterations.
        raise NotImplementedError
