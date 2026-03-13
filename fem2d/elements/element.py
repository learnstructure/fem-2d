import numpy as np


class ElementBase:
    def __init__(self, eid, node_i, node_j):
        self.id = eid
        self.node_i = node_i
        self.node_j = node_j
        self._update_geometry()

    def _update_geometry(self):
        """Compute length, sine, cosine from node coordinates."""
        dx = self.node_j.x - self.node_i.x
        dy = self.node_j.y - self.node_i.y
        self.length = np.hypot(dx, dy)
        self.cos = dx / self.length
        self.sin = dy / self.length

    def transformation_matrix(self, dof_per_node=3):
        """Return 6x6 (or smaller) transformation matrix from local to global."""
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
        """Return element stiffness matrix in local coordinates (size depends on element)."""
        raise NotImplementedError

    def global_stiffness(self):
        """Assemble global stiffness matrix (expanded to 6x6 if necessary)."""
        k_local = self.local_stiffness()
        T = self.transformation_matrix()
        # For elements with fewer than 6 DOF, subclasses must expand k_local and T.
        return T.T @ k_local @ T

    def equivalent_nodal_loads(self):
        """Return equivalent nodal loads due to element loads (global 6‑vector)."""
        return np.zeros(6)  # override if needed

    def internal_forces(self, displacements):
        """Given nodal displacements (global, 6‑vector), return local end forces."""
        # Used for post‑processing and non‑linear iterations.
        raise NotImplementedError
