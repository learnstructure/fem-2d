"""
Materials module defining material properties for structural elements.
"""


class ElasticMaterial:
    """
    Represents an elastic material with stiffness and density properties.

    Attributes
    ----------
    E : float
        Young's Modulus of the material.
    rho : float, optional
        Mass density (mass per unit volume) of the material. Defaults to 0.0.
    """

    def __init__(self, E, rho=0.0):
        """
        Initialize an ElasticMaterial object.

        Parameters
        ----------
        E : float
            Young's Modulus of the material.
        rho : float, optional
            Mass density (mass per unit volume) of the material. Defaults to 0.0.
        """
        self.E = E
        self.rho = rho  # mass density (mass per unit volume)

