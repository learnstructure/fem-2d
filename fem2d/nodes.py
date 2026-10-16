"""
Nodes module representing points in the structural model.
"""


class Node:
    """
    Represents a node in a 2D finite element model.

    Attributes
    ----------
    id : int or str
        Unique identifier of the node.
    x : float
        x-coordinate of the node.
    y : float
        y-coordinate of the node.
    support : list of bool
        Fixity condition of the node's degrees of freedom [ux, uy, rz].
    load : list of float
        Nodal force and moment values [Fx, Fy, Mz].
    dofs : list of int or None
        Global degrees of freedom indices assigned to this node.
    mass : float
        Translational lumped mass of the node.
    inertia : float
        Rotational lumped inertia of the node.
    """

    def __init__(self, nid, x, y):
        """
        Initialize a Node object.

        Parameters
        ----------
        nid : int or str
            Unique identifier of the node.
        x : float
            x-coordinate of the node.
        y : float
            y-coordinate of the node.
        """
        self.id = nid
        self.x = x
        self.y = y
        self.support = [False, False, False]  # [ux, uy, rz] fixed?
        self.load = [0.0, 0.0, 0.0]  # nodal load [Fx, Fy, Mz]
        self.dofs = None  # will be set by Structure
        self.mass = 0.0  # translational lumped mass
        self.inertia = 0.0  # rotational lumped inertia (optional)

    def set_support(self, ux_fixed=False, uy_fixed=False, rz_fixed=False):
        """
        Set support fixity conditions for the degrees of freedom.

        Can be called either with keyword/positional boolean flags::

            node.set_support(ux_fixed=True, uy_fixed=True, rz_fixed=True)
            node.set_support([1, 1, 1])

        or more compactly with three positional ``0``/``1`` (or
        ``True``/``False``) flags ordered as ``(ux, uy, rz)``::

            node.set_support(1, 1, 1)
            node.set_support(0, 1, 0)

        Parameters
        ----------
        ux_fixed : bool, int, or list/tuple, optional
            Whether translation along the x-axis is fixed, or a 3-element list [ux, uy, rz].
        uy_fixed : bool or int, optional
            Whether translation along the y-axis is fixed. Defaults to False.
        rz_fixed : bool or int, optional
            Whether rotation about the z-axis is fixed. Defaults to False.
        """
        if isinstance(ux_fixed, (list, tuple)) and len(ux_fixed) >= 3:
            self.support = [bool(ux_fixed[0]), bool(ux_fixed[1]), bool(ux_fixed[2])]
        else:
            self.support = [bool(ux_fixed), bool(uy_fixed), bool(rz_fixed)]

    def set_load(self, fx=0.0, fy=0.0, mz=0.0):
        """
        Apply forces and moment directly to the node.

        Parameters
        ----------
        fx : float, optional
            Concentrated load in the global x-direction. Defaults to 0.0.
        fy : float, optional
            Concentrated load in the global y-direction. Defaults to 0.0.
        mz : float, optional
            Concentrated moment about the global z-axis. Defaults to 0.0.
        """
        self.load = [fx, fy, mz]

    def set_mass(self, mass=0.0, inertia=0.0):
        """
        Add lumped mass (translational) and optionally rotational inertia.

        Parameters
        ----------
        mass : float, optional
            Translational mass lumped at the node. Defaults to 0.0.
        inertia : float, optional
            Rotational inertia lumped at the node. Defaults to 0.0.
        """
        self.mass = mass
        self.inertia = inertia

