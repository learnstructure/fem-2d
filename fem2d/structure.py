import numpy as np

from fem2d.elements.beam import BeamElement
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
        for load in self.loads:
            if isinstance(load, DistributedLoad):
                eq = load.equivalent_nodal_loads()
                dofs = load.element.node_i.dofs + load.element.node_j.dofs
                self.F[dofs] += eq
            # point loads already handled via node.load

    def apply_boundary_conditions(self):
        """Determine free and fixed DOFs, partition matrices."""
        self._auto_fix_rotations()  # ensure rotations are fixed at truss-only nodes
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

    def _auto_fix_rotations(self):
        """Automatically fix rotational DOF at nodes that have only truss elements."""
        # First, identify nodes that have at least one beam element
        nodes_with_beam = set()
        for el in self.elements.values():
            # Check if element is a beam (you need a way to identify element type)
            if hasattr(el, "inertia") or isinstance(
                el, BeamElement
            ):  # adjust as needed
                nodes_with_beam.add(el.node_i)
                nodes_with_beam.add(el.node_j)

        # For all other nodes, set rotation as fixed
        for node in self.nodes.values():
            if node not in nodes_with_beam:
                # Override user's setting for rotational DOF
                node.support[2] = True  # fix rotation

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
