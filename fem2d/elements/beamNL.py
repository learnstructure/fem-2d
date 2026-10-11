"""
Non-linear beam element wrapper for geometrically non-linear 2D beam elements.
"""

from fem2d.elements.beam import BeamElement


class BeamElementNL(BeamElement):
    """
    Geometrically nonlinear 2D beam element using corotational formulation.
    Inherits from BeamElement with geometric_nonlinear=True by default for backward compatibility.
    """

    def __init__(self, eid, node_i, node_j, material, section, mass=0.0):
        """
        Initialize a BeamElementNL object.
        """
        super().__init__(
            eid=eid,
            node_i=node_i,
            node_j=node_j,
            material=material,
            section=section,
            mass=mass,
            geometric_nonlinear=True,
        )
