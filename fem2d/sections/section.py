"""
Base Section class defining cross-sectional properties.
"""

import math
import numpy as np


class Section:
    """
    Represents a cross-section of a structural element: its area, moment of inertia, and material.

    Elements take a section and read their material from it (``section.material``), so the same shape in two
    materials is two sections.

    Attributes
    ----------
    A : float
        Cross-sectional area.
    I : float
        Moment of inertia (second moment of area).
    material : ElasticMaterial or UniaxialMaterial or None
        Material of the section.
    """

    def __init__(self, area, moi=0.0, As=None, fs=None, material=None):
        """
        Initialize a Section object.

        Parameters
        ----------
        area : float
            Cross-sectional area.
        moi : float, optional
            Moment of inertia (second moment of area). Defaults to 0.0.
        As : float, optional
            Effective shear area. Defaults to None.
        fs : float, optional
            Dimensionless shear shape factor (e.g. 1.2 for rectangle, 10/9 for circle).
            Used to compute As = A / fs if As is not specified.
        material : ElasticMaterial or UniaxialMaterial, optional
            Material of the section. Elements built on the section require it.
        """
        self.material = material
        self.A = float(area)
        self.I = float(moi)
        if As is not None:
            self.As = float(As)
            self.fs = float(fs) if fs is not None else (self.A / self.As if self.As > 0 else 1.0)
        elif fs is not None:
            self.fs = float(fs)
            self.As = self.A / self.fs if self.fs > 0 else None
        else:
            self.As = None
            self.fs = None

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
    def from_rectangle(cls, width, depth, material=None):
        """
        Create a rectangular section with exact inertia and shear area (fs = 1.2), of ``material``.
        """
        area = width * depth
        moi = (width * depth**3) / 12.0
        # For rectangular cross-sections, shape factor fs = 1.2 => As = A / 1.2 = 5/6 A
        as_val = (5.0 / 6.0) * area
        return cls(area, moi, As=as_val, fs=1.2, material=material)

    @classmethod
    def from_circle(cls, diameter, material=None):
        """
        Create a circular section with exact inertia and shear area (fs = 10/9), of ``material``.
        """
        radius = diameter / 2.0
        area = math.pi * radius**2
        moi = math.pi * radius**4 / 4.0
        # For circular cross-sections, shape factor fs = 10/9 => As = 0.9 A
        as_val = 0.9 * area
        return cls(area, moi, As=as_val, fs=10.0 / 9.0, material=material)
