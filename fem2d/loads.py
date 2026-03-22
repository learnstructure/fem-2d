import numpy as np


class PointLoad:
    """Represents a load acting on a node."""

    def __init__(self, node, fx=0, fy=0, mz=0):
        self.node = node
        self.fx, self.fy, self.mz = fx, fy, mz


class ElementLoad:
    """Base class for loads acting on an element."""

    def __init__(self, element):
        self.element = element

    def _compute_equivalent_loads(self):
        raise NotImplementedError("Subclasses must implement this method")

    def _transform_and_store_equivalent_loads(self, eq_local):
        """
        Transform local equivalent nodal loads to global coordinates and
        store them in the element.
        """
        c = self.element.cos
        s = self.element.sin
        T = np.array(
            [
                [c, -s, 0, 0, 0, 0],
                [s, c, 0, 0, 0, 0],
                [0, 0, 1, 0, 0, 0],
                [0, 0, 0, c, -s, 0],
                [0, 0, 0, s, c, 0],
                [0, 0, 0, 0, 0, 1],
            ]
        )
        eq_global = T @ eq_local
        if not hasattr(self.element, "eq_load"):
            self.element.eq_load = np.zeros(6)
        self.element.eq_load += eq_global


class DistributedLoad(ElementLoad):
    """
    Represents a uniformly distributed load acting on an element.
    The load is defined in the element's local coordinate system:
        wx: force per unit length along local x (axial)
        wy: force per unit length along local y (transverse)
    """

    def __init__(self, element, wx=0.0, wy=0.0):
        super().__init__(element)
        self.wx = wx
        self.wy = wy
        self._compute_equivalent_loads()

    def _compute_equivalent_loads(self):
        L = self.element.length
        # Standard fixed-end forces for uniform load (local coordinates)
        eq_local = np.array(
            [
                self.wx * L / 2,
                self.wy * L / 2,
                self.wy * L**2 / 12,
                self.wx * L / 2,
                self.wy * L / 2,
                -self.wy * L**2 / 12,
            ]
        )
        self._transform_and_store_equivalent_loads(eq_local)


class ElementPointLoad(ElementLoad):
    """
    Represents a point load acting on an element at a specific distance from its start node.
    The load is defined in the element's local coordinate system:
        px: force along local x (axial)
        py: force along local y (transverse)
        mz: moment around local z
        x: distance from the start node of the element.
    """

    def __init__(self, element, px=0.0, py=0.0, mz=0.0, x=0.0):
        super().__init__(element)
        self.px = px
        self.py = py
        self.mz = mz
        if x > element.length or x < 0:
            raise ValueError("Load position must be within the element length.")
        self.x = x
        self._compute_equivalent_loads()

    def _compute_equivalent_loads(self):
        L = self.element.length
        a = self.x
        b = L - a

        # Fixed-end forces for a point load in local coordinates
        FAx = (self.px * b) / L
        FBx = (self.px * a) / L
        FAy = (self.py * b**2 * (3 * a + b)) / L**3
        FBy = (self.py * a**2 * (a + 3 * b)) / L**3
        MAz = (self.py * a * b**2) / L**2
        MBz = -(self.py * a**2 * b) / L**2

        # Add moment contributions at nodes
        MAz += (self.mz * b) / L
        MBz += (self.mz * a) / L

        eq_local = np.array([FAx, FAy, MAz, FBx, FBy, MBz])
        self._transform_and_store_equivalent_loads(eq_local)


class TriangularLoad(ElementLoad):
    """
    Represents a trapezoidally distributed load (UVL) on an element.
    This can be used for triangular loads by setting one of w1 or w2 to zero.
    The load is defined in the element's local coordinate system.
    w1 and w2 are load intensities at the start and end of the element.
    Can be axial (wx1, wx2) or transverse (wy1, wy2).
    """

    def __init__(self, element, wx1=0.0, wx2=0.0, wy1=0.0, wy2=0.0):
        super().__init__(element)
        self.wx1 = wx1
        self.wx2 = wx2
        self.wy1 = wy1
        self.wy2 = wy2
        self._compute_equivalent_loads()

    def _compute_equivalent_loads(self):
        L = self.element.length
        eq_local = np.zeros(6)

        # Axial forces (trapezoidal)
        eq_local[0] += (L / 6) * (2 * self.wx1 + self.wx2)
        eq_local[3] += (L / 6) * (self.wx1 + 2 * self.wx2)

        # Transverse forces (trapezoidal), using superposition
        # 1. Uniformly distributed load part (wy1)
        eq_local[1] += self.wy1 * L / 2
        eq_local[2] += self.wy1 * L**2 / 12
        eq_local[4] += self.wy1 * L / 2
        eq_local[5] -= self.wy1 * L**2 / 12

        # 2. Uniformly varying load part (delta from wy1 to wy2)
        wy_delta = self.wy2 - self.wy1
        eq_local[1] += (3 * wy_delta * L) / 20
        eq_local[2] += (wy_delta * L**2) / 30
        eq_local[4] += (7 * wy_delta * L) / 20
        eq_local[5] -= (wy_delta * L**2) / 20

        self._transform_and_store_equivalent_loads(eq_local)
