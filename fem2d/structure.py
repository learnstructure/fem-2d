"""
Structure module defining the main model container and global assembly/solver processes.
"""

from typing import Any, Optional, Union
import numpy as np

from fem2d.elements.beam import BeamElement
from fem2d.elements.spring import SpringElement
from .loads import DistributedLoad, PointLoad
from .solver import NewtonRaphsonSolver


class Structure:
    """
    Represents a 2D finite element model container.

    Manages nodes, elements, and loads, and handles degrees of freedom numbering,
    global stiffness/mass matrix assembly, boundary conditions, and solution procedures.

    Attributes
    ----------
    nodes : dict of {int/str: Node}
        Dictionary of nodes mapped by their IDs.
    elements : dict of {int/str: ElementBase}
        Dictionary of elements mapped by their IDs.
    loads : list
        List of load objects applied to the structure.
    K : numpy.ndarray or None
        Global stiffness matrix.
    F : numpy.ndarray or None
        Global force vector.
    M : numpy.ndarray or None
        Global mass matrix.
    disp : numpy.ndarray or None
        Global displacement vector.
    reactions : numpy.ndarray or None
        Global reaction force vector at fixed degrees of freedom.
    neq : int
        Total number of degrees of freedom.
    fixed_dofs : list of int
        List of fixed degrees of freedom indices.
    free_dofs : list of int
        List of free degrees of freedom indices.
    """

    def __init__(self, units: Optional[Union[Any, str]] = None):
        """Initialize an empty Structure."""
        self.units = units
        self.nodes = {}  # id -> Node
        self.elements = {}  # id -> ElementBase
        self.loads = []  # list of load objects
        self.K = None  # global stiffness
        self.F = None  # global load vector
        self.M = None  # global mass matrix
        self.disp = None  # displacement vector
        self.reactions = None
        self.modal_results = None  # stores (omega, phi) from modal analysis

    def add_node(self, node):
        """
        Add a Node to the structure.

        Parameters
        ----------
        node : Node
            The Node object to add.
        """
        self.nodes[node.id] = node

    def add_element(self, element):
        """
        Add an element to the structure.

        Parameters
        ----------
        element : ElementBase
            The Element object to add.
        """
        element.structure = self
        self.elements[element.id] = element

    def add_load(self, load):
        """
        Add a load to the structure.

        Parameters
        ----------
        load : PointLoad or ElementLoad
            The load object to add.
        """
        self.loads.append(load)
    def set_modal_results(self, omega, phi):
        """
        Register modal analysis results on the structure for later use by plotting utilities.

        Once registered, ``plot_mode_shape`` (and ``DrawStructure.draw_mode_shape``) can be
        called without explicitly passing ``phi`` or ``period`` — they are read from here.

        Parameters
        ----------
        omega : numpy.ndarray
            Array of natural circular frequencies (rad/s), one per mode.
        phi : numpy.ndarray
            Mode shape matrix of shape ``(n_free_dofs, n_modes)``.  Each column is a mode.
        """
        import numpy as np
        self.modal_results = {"omega": np.asarray(omega), "phi": np.asarray(phi)}

    def number_dofs(self):
        """Assign DOF indices to each node (3 DOFs per node: [ux, uy, rz])."""
        dof = 0
        for nid in sorted(self.nodes):
            self.nodes[nid].dofs = [dof, dof + 1, dof + 2]
            dof += 3
        self.neq = dof

    def assemble_stiffness(self):
        """Assemble the global stiffness matrix from all element stiffness contributions."""
        self.K = np.zeros((self.neq, self.neq))
        for el in self.elements.values():
            k_global = el.global_stiffness()
            dofs = el.node_i.dofs + el.node_j.dofs
            for i, ii in enumerate(dofs):
                for j, jj in enumerate(dofs):
                    self.K[ii, jj] += k_global[i, j]

    def assemble_loads(self):
        """Assemble the global force vector from nodal forces, point loads, and element equivalent loads."""
        self.F = np.zeros(self.neq)
        # nodal loads stored on nodes
        for node in self.nodes.values():
            if any(node.load):
                self.F[node.dofs] += node.load

        # point loads stored in self.loads
        for load in self.loads:
            if isinstance(load, PointLoad):
                dofs = load.node.dofs
                self.F[dofs] += np.array([load.fx, load.fy, load.mz])

        # element loads (distributed, etc.)
        for el in self.elements.values():
            if hasattr(el, "eq_load") and el.eq_load is not None:
                dofs = el.node_i.dofs + el.node_j.dofs
                self.F[dofs] += el.eq_load

    def apply_boundary_conditions(self):
        """
        Determine free and fixed DOFs, partitioning matrices accordingly.
        Also automatically flags unstable DOFs before partitioning.
        """
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
        """
        Automatically fix rotational and unstable translation DOFs at nodes
        that are connected only to certain element types (e.g., truss or spring elements).
        """
        for node in self.nodes.values():
            connected_elements = [
                el for el in self.elements.values() if node in (el.node_i, el.node_j)
            ]

            is_beam_connected = any(
                "Beam" in type(el).__name__ for el in connected_elements
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

    def solve(
        self,
        geometric_nonlinear=False,
        material_nonlinear=False,
        n_steps=1,
        tolerance=1e-8,
        max_iter=30,
    ):
        """
        Perform static analysis on the structure (linear or non-linear).

        Parameters
        ----------
        geometric_nonlinear : bool, optional
            Whether to include large displacement geometric non-linear effects. Defaults to False.
        material_nonlinear : bool, optional
            Whether material non-linear behavior is present. Defaults to False.
        n_steps : int, optional
            Number of load increments for non-linear analysis. Defaults to 1.
        tolerance : float, optional
            Convergence tolerance for non-linear solver. Defaults to 1e-8.
        max_iter : int, optional
            Maximum iterations per load step in non-linear solver. Defaults to 30.
        """
        if geometric_nonlinear:
            for el in self.elements.values():
                if hasattr(el, "geometric_nonlinear"):
                    el.geometric_nonlinear = True

        is_nonlinear = (
            geometric_nonlinear
            or material_nonlinear
            or n_steps > 1
            or any(
                getattr(el, "geometric_nonlinear", False)
                for el in self.elements.values()
            )
            or any(
                hasattr(el, "section") and hasattr(el.section, "fibers")
                for el in self.elements.values()
            )
        )

        if is_nonlinear:
            self.solve_nonlinear(
                tolerance=tolerance, max_iter=max_iter, n_steps=n_steps
            )
        else:
            self._solve_linear()

    def _solve_linear(self):
        """
        Perform linear static analysis.
        Numbers DOFs, assembles matrices, partitions, solves for free displacements,
        and computes support reaction forces.
        """
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

    def solve_nonlinear(self, tolerance=1e-8, max_iter=30, n_steps=1):
        """
        Perform non-linear static analysis using Newton-Raphson.

        Supports geometric and/or material non-linear behavior. For path-dependent
        material models, the load is applied incrementally in n_steps.

        Parameters
        ----------
        tolerance : float, optional
            Convergence tolerance. Defaults to 1e-8.
        max_iter : int, optional
            Maximum iterations. Defaults to 30.
        n_steps : int, optional
            Number of load increments. Defaults to 1.
        """
        # Ensure DOFs are numbered and boundary conditions are known
        self.number_dofs()
        self.apply_boundary_conditions()
        # Assemble external load vector from nodal and element loads
        self.assemble_loads()
        P_ext_full = self.F.copy()  # external loads (full vector)

        # Create solver
        from .solver import NewtonRaphsonSolver

        self.disp = np.zeros(self.neq)

        # Incremental load stepping
        for step in range(n_steps):
            factor = (step + 1) / n_steps
            P_ext_step = P_ext_full * factor
            if n_steps > 1:
                print(f"--- Load Step {step+1}/{n_steps} (Load Factor: {factor:.2f}) ---")

            solver = NewtonRaphsonSolver(self, tolerance, max_iter)
            self.disp = solver.solve(P_ext_step, d_init=self.disp)

            # Commit element/material state variables for this step
            self.commit_state()

        # Compute reactions after convergence of the final step
        self._compute_reactions_nonlinear()

    def commit_state(self):
        """Commit history/state variables of all elements (e.g. plastic strains)."""
        for el in self.elements.values():
            if hasattr(el, "commit_state"):
                el.commit_state()

    def _compute_reactions_nonlinear(self):
        """Assemble internal forces at fixed DOFs to obtain reaction forces."""
        # Sum internal forces from all elements
        f_int = np.zeros(self.neq)
        for el in self.elements.values():
            dofs = el.node_i.dofs + el.node_j.dofs
            f_int[dofs] += el.get_internal_forces()
        self.reactions = np.zeros(self.neq)
        if self.fixed_dofs:
            self.reactions[self.fixed_dofs] = f_int[self.fixed_dofs]

    def assemble_mass_matrix(self, lumped: bool = False, include_rotational_inertia: bool = False):
        """Assemble the global mass matrix from element and nodal mass contributions."""
        if not hasattr(self, "neq"):
            self.number_dofs()
        self.M = np.zeros((self.neq, self.neq))
        self._mass_lumped = lumped
        self._mass_rotational = include_rotational_inertia

        # Element contributions
        for el in self.elements.values():
            m_el = el.mass_matrix(lumped=lumped, include_rotational_inertia=include_rotational_inertia)
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

    def get_reduced_matrices(self, lumped: Optional[bool] = None, include_rotational_inertia: bool = False):
        """
        Return (K_ff, M_ff) reduced to free degrees of freedom.

        Assumes assemble_stiffness() and assemble_mass_matrix() have been called,
        and apply_boundary_conditions() has been called.

        Returns
        -------
        K_ff : numpy.ndarray
            Reduced global stiffness matrix for free degrees of freedom.
        M_ff : numpy.ndarray
            Reduced global mass matrix for free degrees of freedom.

        Raises
        ------
        ValueError
            If boundary conditions have not been applied yet.
        """
        if self.free_dofs is None:
            raise ValueError(
                "Boundary conditions not applied. Call apply_boundary_conditions() first."
            )
        if self.K is None:
            self.assemble_stiffness()
        if lumped is not None:
            if (
                self.M is None
                or getattr(self, "_mass_lumped", None) != lumped
                or getattr(self, "_mass_rotational", None) != include_rotational_inertia
            ):
                self.assemble_mass_matrix(lumped=lumped, include_rotational_inertia=include_rotational_inertia)
        elif self.M is None:
            self.assemble_mass_matrix(lumped=False, include_rotational_inertia=include_rotational_inertia)
        K_ff = self.K[np.ix_(self.free_dofs, self.free_dofs)]
        M_ff = self.M[np.ix_(self.free_dofs, self.free_dofs)]
        return K_ff, M_ff

    def modal_analysis(
        self,
        num_modes: int = 6,
        lumped: bool = False,
        include_rotational_inertia: bool = False,
    ):
        """
        Perform 2D undamped eigenvalue modal analysis.

        Solves the generalized eigenvalue problem:
        (K_ff - omega^2 * M_ff) * phi = 0
        delegating matrix regularization and eigensolution to structdyn.ModalAnalysis.

        Parameters
        ----------
        num_modes : int, optional
            Number of vibration modes to compute. Defaults to 6.
        lumped : bool, optional
            Whether to use lumped mass formulation. Defaults to False.
        include_rotational_inertia : bool, optional
            Whether to include frame rotational lumped inertia when lumped=True. Defaults to False.

        Returns
        -------
        dict
            Dictionary containing natural frequencies, periods, and mode shapes:
            - 'omega': circular frequencies in rad/s (1D array)
            - 'frequencies': cyclic frequencies in Hz (1D array)
            - 'periods': natural periods in seconds (1D array)
            - 'modes': mode shapes expanded to full global DOF space (neq x n_modes), normalized to peak 1.0
            - 'phi': reduced mode shape matrix on free DOFs (n_free x n_modes)
        """
        from structdyn.mdf.mdf import MDF

        if not hasattr(self, "neq") or self.neq == 0:
            self.number_dofs()
        self.apply_boundary_conditions()

        n_free = len(self.free_dofs)
        modes_to_calc = min(num_modes, n_free)
        if modes_to_calc <= 0:
            self.modal_results = {
                "omega": np.array([], dtype=float),
                "frequencies": np.array([], dtype=float),
                "periods": np.array([], dtype=float),
                "modes": np.zeros((self.neq, 0), dtype=float),
                "phi": np.zeros((n_free, 0), dtype=float),
            }
            return self.modal_results

        mdf = MDF.from_fem2d(
            self,
            lumped=lumped,
            include_rotational_inertia=include_rotational_inertia,
        )
        omega, phi = mdf.modal.modal_analysis(n_modes=modes_to_calc)

        frequencies = omega / (2.0 * np.pi)
        periods = np.zeros_like(frequencies)
        non_zero = frequencies > 1e-10
        periods[non_zero] = 1.0 / frequencies[non_zero]

        full_modes = np.zeros((self.neq, len(omega)), dtype=float)
        for m_idx in range(len(omega)):
            full_modes[self.free_dofs, m_idx] = phi[:, m_idx]
            max_val = np.max(np.abs(full_modes[:, m_idx]))
            if max_val > 1e-12:
                full_modes[:, m_idx] /= max_val

        self.modal_results = {
            "omega": omega,
            "frequencies": frequencies,
            "periods": periods,
            "modes": full_modes,
            "phi": phi,
        }
        return self.modal_results

    def response_spectrum_analysis(
        self,
        params=None,
        spectrum_curve=None,
        gm=None,
        gm_component=None,
        Sa=None,
        Sd=None,
        n_modes=None,
        ji=None,
        direction="x",
        scale_factor=1.0,
        combination="CQC",
        lumped: bool = True,
        include_rotational_inertia: bool = False,
    ):
        """
        Perform linear dynamic Response Spectrum Analysis (RSA) on the structure.

        Parameters
        ----------
        params : SeismicParameters, optional
            SeismicParameters object from struct_load providing design response spectrum.
        spectrum_curve : callable or tuple or dict, optional
            Continuous response spectrum curve Sa(T).
        gm : GroundMotion or str, optional
            Ground motion object, built-in event name (e.g. 'elcentro', 'northridge',
            'loma_prieta', 'san_fernando', 'elcentro_chopra'), or composite string (e.g. 'northridge:hor2').
        gm_component : str, optional
            Specific component when gm is a string ('hor1', 'hor2', 'up').
            If None, automatically inferred from excitation direction:
            - direction='x' -> 'hor1'
            - direction='y' -> 'hor2'
        Sa : array-like, optional
            Array of modal spectral accelerations.
        Sd : array-like, optional
            Array of modal spectral displacements.
        n_modes : int, optional
            Number of modes to include. Defaults to all available free DOFs.
        ji : float or array-like, optional
            Modal damping ratio(s). If None, defaults to params.damping_ratio or 0.05 (5%).
        direction : str, optional
            Excitation direction ('x' or 'y'). Defaults to 'x'.
        scale_factor : float, optional
            Scale factor for spectral values (e.g. 9.80665 or 386.09 if Sa is in g).
        combination : str, optional
            Modal combination method: 'CQC', 'SRSS', or 'ABS'. Defaults to 'CQC'.
        lumped : bool, optional
            Whether to use lumped mass formulation. Defaults to True (commercial standard for RSA).
        include_rotational_inertia : bool, optional
            Whether to include rotational mass inertia when lumped=True. Defaults to False.

        Returns
        -------
        dict
            Dictionary containing rsa object, modal and combined displacements,
            element forces, base shear, and mass participation metrics.
        """
        from structdyn.mdf.mdf import MDF
        from structdyn.mdf.analytical_methods.response_spectrum_analysis import (
            ResponseSpectrumAnalysis,
        )

        if not hasattr(self, "neq") or self.neq == 0:
            self.number_dofs()
        self.apply_boundary_conditions()

        if ji is None:
            ji = getattr(params, "damping_ratio", 0.05) if params is not None else 0.05

        if spectrum_curve is None and params is not None and getattr(params, "sa_spectrum", None) is not None:
            spectrum_curve = params.sa_spectrum

        mdf = MDF.from_fem2d(
            self,
            zeta=ji,
            n_modes=n_modes,
            lumped=lumped,
            include_rotational_inertia=include_rotational_inertia,
        )
        rsa = ResponseSpectrumAnalysis(
            mdf=mdf,
            ji=ji,
            n_modes=n_modes,
            gm=gm,
            gm_component=gm_component,
            Sd=Sd,
            Sa=Sa,
            spectrum_curve=spectrum_curve,
            params=params,
            direction=direction,
            scale_factor=scale_factor,
        )

        modal_u = rsa.modal_displacements()
        combined_u = rsa.combine_modal_responses(modal_u, method=combination)

        d_comb = np.zeros(self.neq, dtype=float)
        d_comb[self.free_dofs] = combined_u
        self.disp = d_comb

        elem_forces = rsa.compute_modal_element_forces(method=combination)
        modal_vb = rsa.modal_base_shear()
        combined_vb = rsa.combine_modal_responses(modal_vb, method=combination)

        self.set_modal_results(rsa.omega, rsa.phi)

        self.rsa_results = {
            "rsa": rsa,
            "omega": rsa.omega,
            "periods": rsa.T,
            "phi": rsa.phi,
            "participation_factors": rsa.lambda_,
            "effective_mass": rsa.effective_mass,
            "mass_participation_ratio": rsa.mass_participation_ratio,
            "modal_displacements": modal_u,
            "combined_displacements": combined_u,
            "modal_element_forces": elem_forces["modal"],
            "combined_element_forces": elem_forces["combined"],
            "modal_base_shear": modal_vb,
            "combined_base_shear": combined_vb,
            "combination_method": combination,
        }
        return self.rsa_results

    def apply_self_weight(
        self,
        direction: str = "-y",
        gravity_accel: Optional[float] = None,
        units: Optional[Any] = None,
    ) -> "Structure":
        """
        Compute and apply element self-weight distributed loads directly to the structure.

        Parameters
        ----------
        direction : str, optional
            Gravity direction in global coordinates ('-y' or '-x'). Defaults to '-y'.
        gravity_accel : float, optional
            Acceleration due to gravity. If None, derived from structure.units.
        units : str or Units, optional
            Units for calculation. If None, uses structure.units.

        Returns
        -------
        Structure
            self (for chaining).
        """
        from struct_load import apply_self_weight as _apply_sw
        return _apply_sw(self, direction=direction, gravity_accel=gravity_accel, units=units)

