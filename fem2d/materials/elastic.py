"""
Elastic material model.
"""

from fem2d.materials.material import UniaxialMaterial


class ElasticMaterial(UniaxialMaterial):
    """
    Represents an elastic material with stiffness and density properties.

    Attributes
    ----------
    E : float
        Young's Modulus of the material.
    rho : float, optional
        Mass density (mass per unit volume) of the material. Defaults to 0.0.
    """

    def __init__(self, E, rho=0.0, nu=0.3, G=None):
        """
        Initialize an ElasticMaterial object.

        Parameters
        ----------
        E : float
            Young's Modulus of the material.
        rho : float, optional
            Mass density (mass per unit volume) of the material. Defaults to 0.0.
        nu : float, optional
            Poisson's ratio. Defaults to 0.3.
        G : float, optional
            Shear modulus. If None, derived as G = E / (2 * (1 + nu)).
        """
        self.E = float(E)
        self.rho = float(rho)
        self.nu = float(nu)
        if G is not None:
            self.G = float(G)
        else:
            self.G = self.E / (2.0 * (1.0 + self.nu)) if self.nu != -1.0 else 0.0

    def get_response(self, deformation, state=None):
        """
        Compute the force/stress and tangent stiffness at a given deformation/strain.

        For elastic material, response is linear: force = E * strain.
        """
        return self.E * deformation, self.E, state

    def get_initial_tangent(self):
        """Return the initial elastic tangent (Young's Modulus E)."""
        return self.E
