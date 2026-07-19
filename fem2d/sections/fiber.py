"""
Fiber section discretization for material nonlinear analysis.
"""

import numpy as np
from fem2d.sections.section import Section


class FiberSection(Section):
    """
    Represents a cross-section discretized into 1D fibers (layers).
    Each fiber has a y-coordinate, cross-sectional area, and uniaxial material model.
    """

    def __init__(self, fibers):
        """
        Initialize the FiberSection.

        Parameters
        ----------
        fibers : list of dict
            List of fibers, where each fiber is a dictionary:
            {"y": float, "area": float, "material": UniaxialMaterial}
        """
        self.fibers = fibers

        # Compute elastic equivalent area and moment of inertia for compatibility
        area = sum(f["area"] for f in fibers)
        moi = sum(f["area"] * (f["y"] ** 2) for f in fibers)

        super().__init__(area, moi)

    def get_response(self, deformation, state=None):
        """
        Compute section force resultants (axial force N, bending moment M)
        and section tangent stiffness matrix by integrating fiber responses.

        Parameters
        ----------
        deformation : numpy.ndarray
            Deformation vector [axial_strain, curvature] (size 2).
        state : list of dict, optional
            Committed history variables for each fiber material. If None,
            initializes default states.

        Returns
        -------
        forces : numpy.ndarray
            Section force vector [N, M] (size 2).
        tangent : numpy.ndarray
            Section tangent stiffness matrix (size 2x2).
        new_state : list of dict
            Updated history states for all fibers.
        """
        eps_0, kappa = deformation

        N = 0.0
        M = 0.0
        k11 = 0.0
        k12 = 0.0
        k22 = 0.0

        if state is None:
            state = [None] * len(self.fibers)

        new_state = []

        for i, fiber in enumerate(self.fibers):
            y = fiber["y"]
            area = fiber["area"]
            material = fiber["material"]

            # Strain at fiber coordinate: eps = eps_0 - y * kappa
            eps_fiber = eps_0 - y * kappa

            # Query material response
            stress, E_tan, fiber_state = material.get_response(eps_fiber, state[i])
            new_state.append(fiber_state)

            # Integrate forces
            N += stress * area
            M += -y * stress * area

            # Integrate stiffness
            k11 += E_tan * area
            k12 += -y * E_tan * area
            k22 += (y ** 2) * E_tan * area

        forces = np.array([N, M])
        tangent = np.array([[k11, k12], [k12, k22]])

        return forces, tangent, new_state

    @classmethod
    def from_rectangular_fibers(cls, width, depth, material, n_fibers=20):
        """
        Create a rectangular fiber section.

        Parameters
        ----------
        width : float
            Width of the rectangular cross-section.
        depth : float
            Depth of the rectangular cross-section.
        material : UniaxialMaterial
            Material definition for the fibers.
        n_fibers : int, optional
            Number of layers to divide the depth into. Defaults to 20.

        Returns
        -------
        FiberSection
            A FiberSection object with fibers distributed along the depth.
        """
        fibers = []
        dy = depth / n_fibers
        area_fiber = width * dy

        # Discretize from -depth/2 to depth/2
        # Centroids of the fibers:
        y_coords = np.linspace(-depth / 2 + dy / 2, depth / 2 - dy / 2, n_fibers)

        for y in y_coords:
            fibers.append({"y": y, "area": area_fiber, "material": material})

        return cls(fibers)
