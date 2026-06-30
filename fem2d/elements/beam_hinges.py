"""
Beam with hinges module defining beam elements with internal moment releases.
"""

from fem2d.elements.beam import BeamElement
import numpy as np



class BeamWithHingesElement(BeamElement):
    """
    Beam element that allows specifying internal moment releases (hinges) at ends i and/or j.

    Attributes
    ----------
    hinge_i : bool
        Whether end i has a moment release (hinge).
    hinge_j : bool
        Whether end j has a moment release (hinge).
    """

    def __init__(
        self, eid, node_i, node_j, material, area, inertia, hinge_i=False, hinge_j=False
    ):
        """
        Initialize a BeamWithHingesElement.

        Parameters
        ----------
        eid : int or str
            Unique identifier of the element.
        node_i : Node
            Start node.
        node_j : Node
            End node.
        material : ElasticMaterial
            Material definition.
        area : float
            Cross-sectional area.
        inertia : float
            Moment of inertia.
        hinge_i : bool, optional
            Whether end i has a hinge. Defaults to False.
        hinge_j : bool, optional
            Whether end j has a hinge. Defaults to False.
        """
        super().__init__(eid, node_i, node_j, material, area, inertia)
        self.hinge_i = hinge_i
        self.hinge_j = hinge_j

    def local_stiffness(self):
        """
        Compute local stiffness matrix taking moment releases into account.

        Returns
        -------
        numpy.ndarray
            6x6 local stiffness matrix.
        """
        E = self.material.E
        A = self.area
        I = self.inertia
        L = self.length

        EA_L = E * A / L
        EI_L3 = E * I / L**3
        EI_L2 = E * I / L**2
        EI_L = E * I / L

        k_beam = super().local_stiffness()
        if not self.hinge_i and not self.hinge_j:       #if no hinges, return standard beam stiffness
            return k_beam

        if self.hinge_i and self.hinge_j:               #if both hinges, return zero for shear and moment DOFs
            k_beam[1:3, :] = 0.0
            k_beam[:, 1:3] = 0.0
            k_beam[4:6, :] = 0.0
            k_beam[:, 4:6] = 0.0
            return k_beam

        k_hinge_i = np.zeros((6, 6))
        k_hinge_i[0, 0] = EA_L
        k_hinge_i[0, 3] = -EA_L
        k_hinge_i[3, 0] = -EA_L
        k_hinge_i[3, 3] = EA_L

        k_hinge_i[1, 1] = 3.0 * EI_L3
        k_hinge_i[1, 4] = -3.0 * EI_L3
        k_hinge_i[1, 5] = 3.0 * EI_L2

        k_hinge_i[4, 1] = -3.0 * EI_L3
        k_hinge_i[4, 4] = 3.0 * EI_L3
        k_hinge_i[4, 5] = -3.0 * EI_L2

        k_hinge_i[5, 1] = 3.0 * EI_L2
        k_hinge_i[5, 4] = -3.0 * EI_L2
        k_hinge_i[5, 5] = 3.0 * EI_L

        if self.hinge_i:                                #if hinge_i only, return k_hinge_i
            return k_hinge_i

        # hinge_j only -> swap node DOFs from hinge_i form
        if self.hinge_j:
            k_hinge_j = np.zeros((6, 6))
            k_hinge_j[0, 0] = EA_L
            k_hinge_j[0, 3] = -EA_L
            k_hinge_j[3, 0] = -EA_L
            k_hinge_j[3, 3] = EA_L

            k_hinge_j[1, 1] = 3.0 * EI_L3
            k_hinge_j[1, 2] = 3.0 * EI_L2
            k_hinge_j[1, 4] = -3.0 * EI_L3

            k_hinge_j[2, 1] = 3.0 * EI_L3
            k_hinge_j[2, 2] = 3.0 * EI_L
            k_hinge_j[2, 4] = -3.0 * EI_L2  

            k_hinge_j[4, 1] = -3.0 * EI_L3
            k_hinge_j[4, 2] = -3.0 * EI_L2
            k_hinge_j[4, 4] = 3.0 * EI_L3
            
            return k_hinge_j

    # def equivalent_nodal_loads(self):
    #     """
    #     Return local equivalent nodal loads due to a uniform transverse load
    #     accounting for moment releases at hinged ends.

    #     Returns
    #     -------
    #     numpy.ndarray
    #         Equivalent nodal load vector (6x1) in local coordinates.
    #     """
    #     if not hasattr(self, "w"):
    #         return np.zeros(6)

    #     w = self.w
    #     L = self.length
    #     eq = np.zeros(6)
    #     eq[1] = w * L / 2
    #     eq[4] = w * L / 2

    #     if self.hinge_i and self.hinge_j:
    #         eq[2] = 0.0
    #         eq[5] = 0.0
    #     elif self.hinge_i:
    #         eq[2] = 0.0
    #         eq[5] = -w * L**2 / 8
    #     elif self.hinge_j:
    #         eq[2] = w * L**2 / 8
    #         eq[5] = 0.0
    #     else:
    #         eq[2] = w * L**2 / 12
    #         eq[5] = -w * L**2 / 12

    #     return eq
