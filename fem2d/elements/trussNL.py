"""
Non-linear truss element wrapper for geometrically non-linear 2D truss elements.
"""

from fem2d.elements.truss import TrussElement


class TrussElementNL(TrussElement):
    """
    Geometrically nonlinear 2D truss element using corotational formulation.
    Inherits from TrussElement with geometric_nonlinear=True by default for backward compatibility.
    """

    def __init__(self, eid, node_i, node_j, material, section, extra_mass=0.0):
        """
        Initialize a TrussElementNL object.
        """
        super().__init__(
            eid=eid,
            node_i=node_i,
            node_j=node_j,
            material=material,
            section=section,
            extra_mass=extra_mass,
            geometric_nonlinear=True,
        )
