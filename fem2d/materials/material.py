"""
Base class for all uniaxial (1D) constitutive material models.
"""


class UniaxialMaterial:
    """
    Abstract base class for uniaxial material models.
    Can represent stress-strain (sigma-epsilon) or moment-curvature (M-kappa) relationships.
    """

    def get_response(self, deformation, state=None):
        """
        Compute the force/stress and tangent stiffness at a given deformation/strain.

        Parameters
        ----------
        deformation : float
            Current strain (or curvature).
        state : dict, optional
            A dictionary containing history variables for path-dependent materials.
            If None, a default initialized state is used.

        Returns
        -------
        force : float
            Current stress (or moment).
        tangent : float
            Current tangent modulus (or tangent stiffness).
        new_state : dict
            Updated history variables.
        """
        raise NotImplementedError

    def get_initial_tangent(self):
        """
        Return the initial elastic tangent stiffness.

        Returns
        -------
        float
            Initial stiffness/modulus.
        """
        raise NotImplementedError
