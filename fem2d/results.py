import numpy as np


class Results:
    def __init__(self, structure):
        self.structure = structure

    def node_displacement(self, node):
        return self.structure.disp[node.dofs]

    def element_forces(self, element):
        # Get global displacements of element nodes
        u_i = self.structure.disp[element.node_i.dofs]
        u_j = self.structure.disp[element.node_j.dofs]
        u_global = np.concatenate([u_i, u_j])
        # Transform to local
        T = element.transformation_matrix()
        u_local = T @ u_global
        # Compute local forces (k_local * u_local - fixed_end_forces)
        f_local = element.local_stiffness() @ u_local
        # Subtract equivalent loads if any (for consistent sign)
        # ... (simplified)
        return f_local
