"""
Spring element module defining 2D spring elements with axial stiffness.
"""

from fem2d.elements.truss import TrussElement
import numpy as np


class SpringElement(TrussElement):
    """
    Represents a 1D elastic spring element with specified axial stiffness.

    Attributes
    ----------
    k : float
        Axial stiffness coefficient of the spring.
    """

    def __init__(self, eid, node_i, node_j, stiffness):
        """
        Initialize a SpringElement.

        Parameters
        ----------
        eid : int or str
            Unique identifier of the element.
        node_i : Node
            Start node.
        node_j : Node
            End node.
        stiffness : float
            Axial stiffness coefficient of the spring.
        """
        self.k = stiffness
        # A spring has no physical section: its stiffness is k. The truss base class gets a unit, massless section
        # it never uses for stiffness.
        from fem2d.materials.elastic import ElasticMaterial
        from fem2d.sections.section import Section
        super().__init__(eid, node_i, node_j, Section(area=1.0, material=ElasticMaterial(E=1.0)))

    def local_stiffness(self):
        """
        Return the 4x4 spring stiffness matrix in local coordinates.

        Returns
        -------
        numpy.ndarray
            4x4 stiffness matrix.
        """
        return self.k * np.array(
            [[1, 0, -1, 0], [0, 0, 0, 0], [-1, 0, 1, 0], [0, 0, 0, 0]]
        )
