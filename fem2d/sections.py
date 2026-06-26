"""
Sections module defining cross-sectional properties of structural elements.
"""

import math


class Section:
    """
    Represents a cross-section of a structural element with area and moment of inertia.

    Attributes
    ----------
    A : float
        Cross-sectional area.
    I : float
        Moment of inertia (second moment of area).
    """

    def __init__(self, area, moi):
        """
        Initialize a Section object.

        Parameters
        ----------
        area : float
            Cross-sectional area.
        moi : float
            Moment of inertia (second moment of area).
        """
        self.A = area
        self.I = moi  # moment of inertia


    @classmethod
    def from_rectangle(cls, width, depth):
        """
        Create a rectangular section.
        Parameters
        ----------
        width : float
            Width (horizontal dimension) in consistent units.
        depth : float
            Depth (vertical dimension) in consistent units.
        Returns
        -------
        Section
            A Section object with area = width*depth and
            moment of inertia = width * depth**3 / 12.
        """
        area = width * depth
        moi = (width * depth**3) / 12
        return cls(area, moi)

    @classmethod
    def from_circle(cls, diameter):
        """
        Create a circular section.
        Parameters
        ----------
        diameter : float
            Diameter of the circle in consistent units.
        Returns
        -------
        Section
            A Section object with area = π*(d/2)^2 and
            moment of inertia = π*(d/2)^4 / 4.
        """
        radius = diameter / 2
        area = math.pi * radius**2
        moi = math.pi * radius**4 / 4
        return cls(area, moi)
