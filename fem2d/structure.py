import numpy as np
from .loads import DistributedLoad


class Structure:
    def __init__(self):
        self.nodes = {}  # id -> Node
        self.elements = {}  # id -> ElementBase
        self.loads = []  # list of load objects
        self.K = None  # global stiffness
        self.F = None  # global load vector
        self.disp = None  # displacement vector
        self.reactions = None

    def add_node(self, node):
        self.nodes[node.id] = node

    def add_element(self, element):
        self.elements[element.id] = element

    def add_load(self, load):
        self.loads.append(load)

    def number_dofs(self):
        """Assign DOF indices to each node (3 per node)."""
        dof = 0
        for nid in sorted(self.nodes):
            self.nodes[nid].dofs = [dof, dof + 1, dof + 2]
            dof += 3
        self.neq = dof

    def assemble_stiffness(self):
        self.K = np.zeros((self.neq, self.neq))
        for el in self.elements.values():
            k_global = el.global_stiffness()
            dofs = el.node_i.dofs + el.node_j.dofs
            for i, ii in enumerate(dofs):
                for j, jj in enumerate(dofs):
                    self.K[ii, jj] += k_global[i, j]

    def assemble_loads(self):
        self.F = np.zeros(self.neq)
        # nodal loads
        for node in self.nodes.values():
            if any(node.load):
                self.F[node.dofs] += node.load
        # element loads (distributed, etc.)
        for load in self.loads:
            if isinstance(load, DistributedLoad):
                eq = load.equivalent_nodal_loads()
                dofs = load.element.node_i.dofs + load.element.node_j.dofs
                self.F[dofs] += eq
            # point loads already handled via node.load

    def apply_boundary_conditions(self):
        """Determine free and fixed DOFs, partition matrices."""
        fixed = []
        free = []
        for node in self.nodes.values():
            for i, fixed_flag in enumerate(node.support):
                if fixed_flag:
                    fixed.append(node.dofs[i])
                else:
                    free.append(node.dofs[i])
        self.fixed_dofs = fixed
        self.free_dofs = free

    def solve(self):
        self.number_dofs()
        self.assemble_stiffness()
        self.assemble_loads()
        self.apply_boundary_conditions()

        # Partition
        K_ff = self.K[np.ix_(self.free_dofs, self.free_dofs)]
        F_f = self.F[self.free_dofs]

        # Solve for free displacements
        self.disp = np.zeros(self.neq)
        self.disp[self.free_dofs] = np.linalg.solve(K_ff, F_f)

        # Compute reactions
        self.reactions = np.zeros(self.neq)
        if self.fixed_dofs:
            self.reactions[self.fixed_dofs] = (
                self.K[np.ix_(self.fixed_dofs, self.free_dofs)]
                @ self.disp[self.free_dofs]
                - self.F[self.fixed_dofs]
            )
