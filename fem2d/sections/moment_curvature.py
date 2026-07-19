"""
Moment-curvature section model for direct action-deformation modeling of beam bending.
"""

import numpy as np
from fem2d.sections.section import Section


class MomentCurvatureSection(Section):
    """
    Represents a cross-section characterized directly by a moment-curvature (M-kappa)
    material model for bending, and a linear elastic relation for axial behavior.
    """

    def __init__(self, E, area, moment_curvature_material):
        """
        Initialize the MomentCurvatureSection.

        Parameters
        ----------
        E : float
            Young's Modulus used for the elastic axial stiffness (EA).
        area : float
            Cross-sectional area.
        moment_curvature_material : UniaxialMaterial
            A material model where the 'strain' input is curvature (kappa) and the
            'stress' output is bending moment (M), and the tangent is EI_tan.
        """
        self.E = E
        self.moment_curvature_material = moment_curvature_material

        # For compatibility, moment of inertia I = EI_initial / E
        EI_init = moment_curvature_material.get_initial_tangent()
        moi = EI_init / E if E != 0.0 else 1.0

        super().__init__(area, moi)

    def get_response(self, deformation, state=None):
        """
        Compute section forces [N, M] and tangent stiffness matrix.

        Parameters
        ----------
        deformation : numpy.ndarray
            Deformation vector [axial_strain, curvature] (size 2).
        state : dict, optional
            Committed state of the moment-curvature material.

        Returns
        -------
        forces : numpy.ndarray
            Section forces [N, M] (size 2).
        tangent : numpy.ndarray
            Tangent stiffness matrix (size 2x2).
        new_state : dict
            Updated material state.
        """
        eps_0, kappa = deformation

        # Axial behavior: linear elastic
        N = self.E * self.A * eps_0
        k11 = self.E * self.A

        # Bending behavior: query the moment-curvature material
        M, EI_tan, new_state = self.moment_curvature_material.get_response(kappa, state)

        forces = np.array([N, M])
        tangent = np.array([[k11, 0.0], [0.0, EI_tan]])

        return forces, tangent, new_state
