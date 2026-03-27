import numpy as np

from fem2d.elements.beam import BeamElement
from fem2d.elements.spring import SpringElement
from .loads import DistributedLoad
from .solver import NewtonRaphsonSolver


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
        element.structure = self
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
        for el in self.elements.values():
            if hasattr(el, "eq_load") and el.eq_load is not None:
                dofs = el.node_i.dofs + el.node_j.dofs
                self.F[dofs] += el.eq_load

    def apply_boundary_conditions(self):
        """Determine free and fixed DOFs, partition matrices."""
        self._auto_fix_unstable_dofs()  # ensure rotations are fixed at truss-only nodes
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

    def _auto_fix_unstable_dofs(self):
        """Automatically fix DOFs at nodes connected only to certain element types."""
        for node in self.nodes.values():
            connected_elements = [
                el for el in self.elements.values() if node in (el.node_i, el.node_j)
            ]

            is_beam_connected = any(
                isinstance(el, BeamElement) for el in connected_elements
            )
            if is_beam_connected:
                continue

            is_spring_only = all(
                isinstance(el, SpringElement) for el in connected_elements
            )

            if is_spring_only and connected_elements:
                # Fix rotation and perpendicular force for spring-only nodes
                node.support[1] = True  # Fix perpendicular force
                node.support[2] = True  # Fix rotation
            else:
                # Fix rotation for truss-only nodes (original logic)
                node.support[2] = True

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

    def solve_nonlinear(self, tolerance=1e-8, max_iter=30):
        """
        Perform geometrically nonlinear analysis using Newton‑Raphson.
        Assumes all elements implement the nonlinear interface
        (update_state, get_tangent_stiffness, get_internal_forces).
        """
        # Ensure DOFs are numbered and boundary conditions are known
        self.number_dofs()
        self.apply_boundary_conditions()
        # Assemble external load vector from nodal and element loads
        self.assemble_loads()
        P_ext = self.F.copy()  # external loads (full vector)

        # Create solver and run
        from .solver import NewtonRaphsonSolver  # adjust import as needed

        solver = NewtonRaphsonSolver(self, tolerance, max_iter)
        self.disp = solver.solve(P_ext)

        # Compute reactions after convergence
        self._compute_reactions_nonlinear()

    def _compute_reactions_nonlinear(self):
        """Assemble internal forces at fixed DOFs to obtain reactions."""
        # Sum internal forces from all elements
        f_int = np.zeros(self.neq)
        for el in self.elements.values():
            dofs = el.node_i.dofs + el.node_j.dofs
            f_int[dofs] += el.get_internal_forces()
        self.reactions = np.zeros(self.neq)
        if self.fixed_dofs:
            self.reactions[self.fixed_dofs] = f_int[self.fixed_dofs]

    def assemble_mass_matrix(self):
        """Assemble the global mass matrix."""
        if not hasattr(self, "neq"):
            self.number_dofs()
        self.M = np.zeros((self.neq, self.neq))

        # Element contributions
        for el in self.elements.values():
            m_el = el.mass_matrix()
            dofs = el.node_i.dofs + el.node_j.dofs
            for i, ii in enumerate(dofs):
                for j, jj in enumerate(dofs):
                    self.M[ii, jj] += m_el[i, j]

        # Node lumped masses
        for node in self.nodes.values():
            if node.mass != 0.0:
                self.M[node.dofs[0], node.dofs[0]] += node.mass
                self.M[node.dofs[1], node.dofs[1]] += node.mass
            if node.inertia != 0.0:
                self.M[node.dofs[2], node.dofs[2]] += node.inertia

    def get_reduced_matrices(self):
        """
        Return (K_ff, M_ff) reduced to free degrees of freedom.
        Assumes assemble_stiffness() and assemble_mass_matrix() have been called,
        and apply_boundary_conditions() has been called.
        """
        if self.free_dofs is None:
            raise ValueError(
                "Boundary conditions not applied. Call apply_boundary_conditions() first."
            )
        if self.K is None:
            self.assemble_stiffness()
        if self.M is None:
            self.assemble_mass_matrix()
        K_ff = self.K[np.ix_(self.free_dofs, self.free_dofs)]
        M_ff = self.M[np.ix_(self.free_dofs, self.free_dofs)]
        return K_ff, M_ff
