"""
Beam element module defining 2D elastic Euler-Bernoulli beam elements with mass options.
"""

from typing import Optional, List, Union
from fem2d.elements.element import ElementBase
import numpy as np


class BeamElement(ElementBase):
    """
    Elastic 2D Euler-Bernoulli beam element.

    The element formulation is the classical Euler-Bernoulli beam (shear
    deformation is neglected), so only axial and flexural degrees of freedom
    are considered. Optional distributed ``extra_mass`` and a corotational
    geometrically-nonlinear formulation are supported.

    End releases (hinges) can be specified either via boolean flags:
    - ``hinge_i=True``, ``hinge_j=True`` (moment releases / hinges)
    or via general 3-DOF release vectors:
    - ``releases_i=[u1, u2, theta]``, ``releases_j=[u1, u2, theta]``
    Released DOFs are condensed out via static condensation.

    Attributes
    ----------
    material : ElasticMaterial
        Material definition for the element.
    area : float
        Cross-sectional area.
    inertia : float
        Moment of inertia of the cross-section.
    extra_mass : float
        Additional distributed mass per unit length (e.g. non-structural mass).
    hinge_i : bool
        Whether end i has a moment release (hinge).
    hinge_j : bool
        Whether end j has a moment release (hinge).
    releases_i : list of bool
        End release flags at start node [u1, u2, theta].
    releases_j : list of bool
        End release flags at end node [u1, u2, theta].
    """

    def __init__(
        self,
        eid: Union[int, str],
        node_i,
        node_j,
        material,
        section,
        extra_mass: float = 0.0,
        geometric_nonlinear: bool = False,
        hinge_i: bool = False,
        hinge_j: bool = False,
        releases_i: Optional[List[bool]] = None,
        releases_j: Optional[List[bool]] = None,
        include_shear_deformation: bool = False,
    ):
        """
        Initialize a BeamElement.

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
            Section definition (containing A and I).
        extra_mass : float, optional
            Additional distributed mass per unit length. Defaults to 0.0.
        geometric_nonlinear : bool, optional
            Whether to use corotational geometric non-linear formulation. Defaults to False.
        hinge_i : bool, optional
            Whether start node i has a moment release (hinge). Defaults to False.
        hinge_j : bool, optional
            Whether end node j has a moment release (hinge). Defaults to False.
        releases_i : list of bool, optional
            End release fixity at start node [u1, u2, theta]. True indicates released.
        releases_j : list of bool, optional
            End release fixity at end node [u1, u2, theta]. True indicates released.
        include_shear_deformation : bool, optional
            Whether to include Timoshenko shear deformation in the local stiffness matrix.
            Requires shear area `As` (or shape factor `fs`) on Section and shear modulus `G`
            (or `nu`) on ElasticMaterial. Defaults to False.
        """
        super().__init__(eid, node_i, node_j)

        # Handle Section object or numeric area/inertia fallback
        if hasattr(section, "A") and hasattr(section, "I"):
            self.section = section
            self.area = section.A
            self.inertia = section.I
            self.material = material
        elif hasattr(material, "A") and hasattr(material, "I"):
            self.section = material
            self.area = material.A
            self.inertia = material.I
            self.material = getattr(material, "material", None)
            if isinstance(section, (int, float)):
                extra_mass = float(section)
        elif isinstance(section, (int, float)):
            from fem2d.sections.section import Section
            area_val = float(section)
            moi_val = float(extra_mass) if isinstance(extra_mass, (int, float)) else 0.0
            self.section = Section(area_val, moi_val)
            self.area = area_val
            self.inertia = moi_val
            self.material = material
            extra_mass = 0.0
        else:
            self.material = material
            self.section = section
            self.area = getattr(section, "A", None)
            self.inertia = getattr(section, "I", None)

        self.extra_mass = extra_mass
        self.geometric_nonlinear = geometric_nonlinear
        self.include_shear_deformation = bool(include_shear_deformation)

        # Releases: 3 DOFs per node [u1, u2, theta] (axial, transverse, rotation)
        self.releases_i = list(releases_i) if releases_i is not None else [False, False, False]
        self.releases_j = list(releases_j) if releases_j is not None else [False, False, False]
        if hinge_i:
            self.releases_i[2] = True
        if hinge_j:
            self.releases_j[2] = True

        self.hinge_i = bool(self.releases_i[2])
        self.hinge_j = bool(self.releases_j[2])

        # State variables for non-linear updates
        self.L0 = None
        self.alpha0 = None
        self.f_local = None
        self.F_global = None
        self.k_local = None
        self.K_global = None

    def _update_initial_geometry(self):
        """Compute initial length and chord angle."""
        dx = self.node_j.x - self.node_i.x
        dy = self.node_j.y - self.node_i.y
        self.L0 = np.hypot(dx, dy)
        self.alpha0 = np.arctan2(dy, dx)

    def _apply_releases(self, k: np.ndarray) -> np.ndarray:
        """
        Condense out released degrees of freedom via static condensation.

        Parameters
        ----------
        k : numpy.ndarray
            6x6 unreleased element stiffness matrix in local coordinates.

        Returns
        -------
        numpy.ndarray
            6x6 condensed stiffness matrix in local coordinates.
        """
        rel_flags = [False] * 6
        if self.releases_i:
            for i, r in enumerate(self.releases_i):
                if r:
                    rel_flags[i] = True
        if self.releases_j:
            for i, r in enumerate(self.releases_j):
                if r:
                    rel_flags[3 + i] = True

        released_dofs = [idx for idx, flag in enumerate(rel_flags) if flag]
        retained_dofs = [idx for idx, flag in enumerate(rel_flags) if not flag]

        if not released_dofs:
            return k

        # Static condensation:
        # [ K_rr  K_rc ] [ u_r ] = [ F_r ]
        # [ K_cr  K_cc ] [ u_c ]   [ 0   ]
        # K_condensed = K_rr - K_rc @ inv(K_cc) @ K_cr
        K_rr = k[np.ix_(retained_dofs, retained_dofs)]
        K_rc = k[np.ix_(retained_dofs, released_dofs)]
        K_cr = k[np.ix_(released_dofs, retained_dofs)]
        K_cc = k[np.ix_(released_dofs, released_dofs)]

        K_cc_inv = np.linalg.pinv(K_cc)
        K_r_star = K_rr - K_rc @ K_cc_inv @ K_cr

        k_new = np.zeros((6, 6), dtype=float)
        for i_idx, r_i in enumerate(retained_dofs):
            for j_idx, r_j in enumerate(retained_dofs):
                k_new[r_i, r_j] = K_r_star[i_idx, j_idx]
        return k_new

    def local_stiffness(self):
        """
        Return the 6x6 element stiffness matrix in local coordinates,
        accounting for shear deformations (if enabled) and any released
        degrees of freedom.

        Reference:
        Kassimali A. (2022). Matrix Analysis of Structures, Chapter 9,
        Eq. (9.60) & Eq. (9.63), pp. 546-552.

        Returns
        -------
        numpy.ndarray
            6x6 local stiffness matrix.
        """
        E = self.material.E if self.material is not None else 1.0
        A = self.area
        I = self.inertia
        L = self.length

        EA_L = E * A / L

        # Dimensionless shear deformation constant beta_s = 12*E*I / (G*As*L^2) (Eq. 9.60)
        beta_s = 0.0
        if self.include_shear_deformation and L > 0.0:
            as_val = getattr(self.section, "As", None)
            if as_val is None:
                fs_val = getattr(self.section, "fs", None)
                if fs_val is not None and fs_val > 0.0 and A is not None and A > 0.0:
                    as_val = A / fs_val
            g_val = getattr(self.material, "G", None)
            if g_val is None and self.material is not None and hasattr(self.material, "E"):
                nu_val = getattr(self.material, "nu", 0.3)
                g_val = E / (2.0 * (1.0 + nu_val)) if nu_val != -1.0 else 0.0

            if as_val is not None and as_val > 0.0 and g_val is not None and g_val > 0.0:
                beta_s = (12.0 * E * I) / (g_val * as_val * L**2)

        factor = 1.0 + beta_s
        k12 = 12.0 * E * I / (L**3 * factor)
        k6 = 6.0 * E * I / (L**2 * factor)
        k4 = (4.0 + beta_s) * E * I / (L * factor)
        k2 = (2.0 - beta_s) * E * I / (L * factor)

        k = np.array(
            [
                [EA_L, 0.0, 0.0, -EA_L, 0.0, 0.0],
                [0.0, k12, k6, 0.0, -k12, k6],
                [0.0, k6, k4, 0.0, -k6, k2],
                [-EA_L, 0.0, 0.0, EA_L, 0.0, 0.0],
                [0.0, -k12, -k6, 0.0, k12, -k6],
                [0.0, k6, k2, 0.0, -k6, k4],
            ]
        )

        if any(self.releases_i) or any(self.releases_j):
            k = self._apply_releases(k)

        return k

    def update_state(self, global_disp):
        """
        Update internal forces and stiffness matrices given global displacement vector.

        Parameters
        ----------
        global_disp : numpy.ndarray
            Full global displacement vector of the structure.
        """
        if self.L0 is None:
            self._update_initial_geometry()

        if self.geometric_nonlinear:
            # Corotational geometrically non-linear formulation
            u_i = global_disp[self.node_i.dofs]
            u_j = global_disp[self.node_j.dofs]

            x_i = self.node_i.x + u_i[0]
            y_i = self.node_i.y + u_i[1]
            x_j = self.node_j.x + u_j[0]
            y_j = self.node_j.y + u_j[1]

            dx = x_j - x_i
            dy = y_j - y_i
            L = np.hypot(dx, dy)
            alpha = np.arctan2(dy, dx)

            beta = alpha - self.alpha0

            theta1_bar = u_i[2] - beta
            theta2_bar = u_j[2] - beta
            u2 = L - self.L0

            d_local = np.array([0.0, 0.0, theta1_bar, u2, 0.0, theta2_bar])

            E = self.material.E if self.material is not None else 1.0
            A = self.area
            I = self.inertia
            L0 = self.L0

            EA_L = E * A / L0
            f_x2 = EA_L * u2

            k1 = np.zeros((6, 6))
            k1[0, 0] = k1[3, 3] = EA_L
            k1[0, 3] = k1[3, 0] = -EA_L

            EI_L3 = E * I / L0**3
            EI_L2 = E * I / L0**2
            EI_L = E * I / L0

            k1[1, 1] = 12 * EI_L3
            k1[1, 2] = 6 * EI_L2
            k1[1, 4] = -12 * EI_L3
            k1[1, 5] = 6 * EI_L2

            k1[2, 1] = 6 * EI_L2
            k1[2, 2] = 4 * EI_L
            k1[2, 4] = -6 * EI_L2
            k1[2, 5] = 2 * EI_L

            k1[4, 1] = -12 * EI_L3
            k1[4, 2] = -6 * EI_L2
            k1[4, 4] = 12 * EI_L3
            k1[4, 5] = -6 * EI_L2

            k1[5, 1] = 6 * EI_L2
            k1[5, 2] = 2 * EI_L
            k1[5, 4] = -6 * EI_L2
            k1[5, 5] = 4 * EI_L

            factor = f_x2 / (30.0 * L0)
            k2 = np.zeros((6, 6))

            k2[1, 1] = 36
            k2[1, 2] = 3 * L0
            k2[1, 4] = -36
            k2[1, 5] = 3 * L0

            k2[2, 1] = 3 * L0
            k2[2, 2] = 4 * L0**2
            k2[2, 4] = -3 * L0
            k2[2, 5] = -(L0**2)

            k2[4, 1] = -36
            k2[4, 2] = -3 * L0
            k2[4, 4] = 36
            k2[4, 5] = -3 * L0

            k2[5, 1] = 3 * L0
            k2[5, 2] = -(L0**2)
            k2[5, 4] = -3 * L0
            k2[5, 5] = 4 * L0**2

            k2 *= factor
            k_local = k1 + k2
            f_local = k_local @ d_local

            c = np.cos(alpha)
            s = np.sin(alpha)
            R = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])
            T = np.zeros((6, 6))
            T[:3, :3] = R
            T[3:, 3:] = R

            self.f_local = f_local
            self.F_global = T @ f_local
            self.k_local = k_local
            self.K_global = T @ k_local @ T.T
        else:
            # Linear static state update
            u_i = global_disp[self.node_i.dofs]
            u_j = global_disp[self.node_j.dofs]
            u_global = np.concatenate([u_i, u_j])
            T = self.transformation_matrix()
            u_local = T @ u_global
            self.k_local = self.local_stiffness()
            self.f_local = self.k_local @ u_local
            if hasattr(self, "eq_load") and self.eq_load is not None:
                eq_load_local = T @ self.eq_load
                self.f_local -= eq_load_local
            if self.releases_i:
                for i_dof, rel in enumerate(self.releases_i):
                    if rel:
                        self.f_local[i_dof] = 0.0
            if self.releases_j:
                for j_dof, rel in enumerate(self.releases_j):
                    if rel:
                        self.f_local[3 + j_dof] = 0.0
            self.K_global = T.T @ self.k_local @ T
            self.F_global = self.K_global @ u_global

    def get_tangent_stiffness(self):
        """Return global tangent stiffness matrix of the element."""
        if hasattr(self, "K_global") and self.K_global is not None:
            return self.K_global
        return self.global_stiffness()

    def get_internal_forces(self):
        """Return global internal forces vector of the element."""
        if hasattr(self, "F_global") and self.F_global is not None:
            return self.F_global
        return np.zeros(6)

    def get_local_forces(self):
        """Return local end forces vector of the element."""
        if hasattr(self, "f_local") and self.f_local is not None:
            return self.f_local
        if hasattr(self, "structure") and self.structure is not None and self.structure.disp is not None:
            u_i = self.structure.disp[self.node_i.dofs]
            u_j = self.structure.disp[self.node_j.dofs]
            u_global = np.concatenate([u_i, u_j])
            T = self.transformation_matrix()
            u_local = T @ u_global
            f_loc = self.local_stiffness() @ u_local
            if hasattr(self, "eq_load") and self.eq_load is not None:
                eq_load_local = T @ self.eq_load
                f_loc -= eq_load_local
            if self.releases_i:
                for i_dof, rel in enumerate(self.releases_i):
                    if rel:
                        f_loc[i_dof] = 0.0
            if self.releases_j:
                for j_dof, rel in enumerate(self.releases_j):
                    if rel:
                        f_loc[3 + j_dof] = 0.0
            return f_loc
        return np.zeros(6)

    def deformed_shape_points(self, global_disp, n_points=20, scale=1.0):
        """
        Return a list of (x, y) points in global coordinates representing
        the deformed shape of the beam, scaled by `scale`.

        Parameters
        ----------
        global_disp : numpy.ndarray
            Full global displacement vector of the structure.
        n_points : int, optional
            Number of points to generate along the element. Defaults to 20.
        scale : float, optional
            Displacement magnification factor. Defaults to 1.0.

        Returns
        -------
        list of tuple of float
            List of (x, y) coordinates representing the deformed shape.
        """
        # Extract nodal displacements (global)
        u_i = global_disp[self.node_i.dofs]  # [ux_i, uy_i, rz_i]
        u_j = global_disp[self.node_j.dofs]  # [ux_j, uy_j, rz_j]

        # Build 6‑vector of global displacements
        u_global = np.concatenate([u_i, u_j])

        # Transform to local displacements (local = T @ global)
        T = self.transformation_matrix()
        u_local = T @ u_global

        # Local displacements: [u1, v1, θ1, u2, v2, θ2]
        u1, v1, theta1, u2, v2, theta2 = u_local

        L = self.length
        c, s = self.cos, self.sin

        # Generate points along the element
        points = []
        for xi in np.linspace(0, L, n_points):
            # Local coordinate
            x_local = xi
            xi_norm = x_local / L  # ξ

            # Axial displacement (linear interpolation)
            u_axial = u1 * (1 - xi_norm) + u2 * xi_norm

            # Transverse displacement (cubic Hermitian)
            # Shape functions
            N1 = 1 - 3 * xi_norm**2 + 2 * xi_norm**3
            N2 = (xi_norm - 2 * xi_norm**2 + xi_norm**3) * L
            N3 = 3 * xi_norm**2 - 2 * xi_norm**3
            N4 = (-(xi_norm**2) + xi_norm**3) * L

            v_trans = v1 * N1 + theta1 * N2 + v2 * N3 + theta2 * N4

            # Local displacement vector at this point (only translations)
            disp_local = np.array([u_axial, v_trans])

            # Transform back to global
            R = T[:2, :2]
            disp_global = R.T @ disp_local

            # Original global coordinates of this point (linear interpolation)
            x_orig = self.node_i.x + xi_norm * (self.node_j.x - self.node_i.x)
            y_orig = self.node_i.y + xi_norm * (self.node_j.y - self.node_i.y)

            # Deformed coordinates (apply scaled displacement)
            x_def = x_orig + scale * disp_global[0]
            y_def = y_orig + scale * disp_global[1]

            points.append((x_def, y_def))

        return points

    def mass_matrix(self):
        """
        Return the 6x6 global mass matrix (including rotational inertia and extra mass).

        Returns
        -------
        numpy.ndarray
            6x6 global mass matrix.
        """
        L = self.length
        A = self.area
        I = self.inertia
        rho = self.material.rho
        extra = self.extra_mass

        # Effective density for translational part: material + extra
        rho_eff = rho + (extra / A) if A != 0 else rho

        # ---- Translational mass (m1) ----
        factor_trans = rho_eff * A * L / 420.0
        m1 = np.zeros((6, 6))
        m1[0, 0] = m1[3, 3] = 140.0 * factor_trans
        m1[0, 3] = m1[3, 0] = 70.0 * factor_trans

        m1[1, 1] = 156.0 * factor_trans
        m1[1, 2] = 22.0 * L * factor_trans
        m1[1, 4] = 54.0 * factor_trans
        m1[1, 5] = -13.0 * L * factor_trans

        m1[2, 1] = 22.0 * L * factor_trans
        m1[2, 2] = 4.0 * L**2 * factor_trans
        m1[2, 4] = 13.0 * L * factor_trans
        m1[2, 5] = -3.0 * L**2 * factor_trans

        m1[4, 1] = 54.0 * factor_trans
        m1[4, 2] = 13.0 * L * factor_trans
        m1[4, 4] = 156.0 * factor_trans
        m1[4, 5] = -22.0 * L * factor_trans

        m1[5, 1] = -13.0 * L * factor_trans
        m1[5, 2] = -3.0 * L**2 * factor_trans
        m1[5, 4] = -22.0 * L * factor_trans
        m1[5, 5] = 4.0 * L**2 * factor_trans

        # ---- Rotational mass (m2) ----
        factor_rot = rho * I / (30.0 * L)
        m2 = np.zeros((6, 6))
        m2[1, 1] = 36.0 * factor_rot
        m2[1, 2] = 3.0 * L * factor_rot
        m2[1, 4] = -36.0 * factor_rot
        m2[1, 5] = 3.0 * L * factor_rot

        m2[2, 1] = 3.0 * L * factor_rot
        m2[2, 2] = 4.0 * L**2 * factor_rot
        m2[2, 4] = -3.0 * L * factor_rot
        m2[2, 5] = -(L**2) * factor_rot

        m2[4, 1] = -36.0 * factor_rot
        m2[4, 2] = -3.0 * L * factor_rot
        m2[4, 4] = 36.0 * factor_rot
        m2[4, 5] = -3.0 * L * factor_rot

        m2[5, 1] = 3.0 * L * factor_rot
        m2[5, 2] = -(L**2) * factor_rot
        m2[5, 4] = -3.0 * L * factor_rot
        m2[5, 5] = 4.0 * L**2 * factor_rot

        # total local mass matrix
        m_local = m1 + m2

        # Transform to global
        T = self.transformation_matrix()
        m_global = T.T @ m_local @ T
        return m_global

    def axial_force(self):
        """Return the axial force in the beam element in local coordinates.

        Positive value corresponds to compression for buckling sign convention.
        """
        f_local = self.get_local_forces()
        return f_local[0]

    def geometric_stiffness(self, axial_force):
        """Return the 6x6 geometric stiffness matrix (local->global) for the beam.

        Uses the standard Euler-Bernoulli beam geometric stiffness for axial
        compression P (positive in compression):

            Kg_local = (P / (30*L)) * [
                [0,   0,    0,    0,    0,    0],
                [0,  36,   3L,    0,  -36,   3L],
                [0, 3L,  4L^2,   0,  -3L,  -L^2],
                [0,   0,    0,    0,    0,    0],
                [0, -36,  -3L,   0,   36,  -3L],
                [0,  3L,  -L^2,  0,  -3L,  4L^2]
            ]

        The returned matrix is in global coordinates.
        """
        P = axial_force
        L = self.length
        if L == 0:
            return np.zeros((6, 6))

        factor = P / (30.0 * L)
        kg_local = factor * np.array(
            [
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                [0.0, 36.0, 3.0 * L, 0.0, -36.0, 3.0 * L],
                [0.0, 3.0 * L, 4.0 * L * L, 0.0, -3.0 * L, -1.0 * L * L],
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                [0.0, -36.0, -3.0 * L, 0.0, 36.0, -3.0 * L],
                [0.0, 3.0 * L, -1.0 * L * L, 0.0, -3.0 * L, 4.0 * L * L],
            ]
        )

        T = self.transformation_matrix()
        kg_global = T.T @ kg_local @ T
        return kg_global
