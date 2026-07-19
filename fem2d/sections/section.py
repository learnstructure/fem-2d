"""
Base Section class defining cross-sectional properties.
"""

import math
import numpy as np


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
        self.I = moi

    def get_response(self, deformation, state=None):
        """
        Compute section force resultants (axial force N, bending moment M)
        and section tangent stiffness matrix for a given section deformation.

        By default, returns linear elastic response.

        Parameters
        ----------
        deformation : numpy.ndarray
            Deformation vector [axial_strain, curvature] (size 2).
        state : dict, optional
            Committed history variables of the section.

        Returns
        -------
        forces : numpy.ndarray
            Section force vector [N, M] (size 2).
        tangent : numpy.ndarray
            Section tangent stiffness matrix (size 2x2).
        new_state : dict
            Updated section state.
        """
        eps, kappa = deformation
        # Linear elastic: N = E*A*eps, M = E*I*kappa
        # For a pure section (no material associated), we assume unit E or raise Error.
        # But for backward compatibility and general interface:
        # In general, sections used in nonlinear analysis will override this.
        raise NotImplementedError("get_response must be implemented by subclasses for nonlinear analysis.")

    @classmethod
    def from_rectangle(cls, width, depth):
        """Create a rectangular section."""
        area = width * depth
        moi = (width * depth**3) / 12
        return cls(area, moi)

    @classmethod
    def from_circle(cls, diameter):
        """Create a circular section."""
        radius = diameter / 2
        area = math.pi * radius**2
        moi = math.pi * radius**4 / 4
        return cls(area, moi)
