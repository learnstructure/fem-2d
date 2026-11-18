"""
Loads module defining point, distributed, and varying loads for nodes and elements.
"""

import numpy as np

_GAUSS_POINTS = (-np.sqrt(0.6), 0.0, np.sqrt(0.6))
_GAUSS_WEIGHTS = (5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0)


def bending_shape_functions(x, length, phi=0.0):
    """Exact interpolation of a two-node Timoshenko beam at distance ``x`` from node i.

    Returns ``(Nv, Npsi)``: the transverse displacement ``v(x)`` and the cross-section rotation ``psi(x)`` caused by
    a unit value of each end DOF ``(v_i, theta_i, v_j, theta_j)`` with the others zero. They solve the unloaded
    Timoshenko beam (EI psi'' + kGA (v' - psi) = 0, (v' - psi)' = 0), so they are the exact deflected shapes of the
    element whose stiffness uses the same shear parameter ``phi = 12 EI / (G As L^2)``
    (``BeamElement.shear_parameter``; Friedman and Kosmatka 1993, "An improved two-node Timoshenko beam finite
    element", Computers & Structures 47(3)). With ``phi = 0`` they are the Hermite cubics and their slopes.

    By the reciprocal theorem the equivalent nodal loads (negative fixed-end reactions) of a member load are the
    work of the load on these shapes: the integral of w(x) Nv(x) for transverse forces, m Npsi(a) for a moment at a.
    """
    xi = x / length
    d = 1.0 + phi
    nv = (
        (2.0 * xi**3 - 3.0 * xi**2 - phi * xi + 1.0 + phi) / d,
        length * (xi**3 - (2.0 + 0.5 * phi) * xi**2 + (1.0 + 0.5 * phi) * xi) / d,
        (-2.0 * xi**3 + 3.0 * xi**2 + phi * xi) / d,
        length * (xi**3 - (1.0 - 0.5 * phi) * xi**2 - 0.5 * phi * xi) / d,
    )
    npsi = (
        6.0 * (xi**2 - xi) / (d * length),
        (3.0 * xi**2 - (4.0 + phi) * xi + 1.0 + phi) / d,
        -6.0 * (xi**2 - xi) / (d * length),
        (3.0 * xi**2 - (2.0 - phi) * xi) / d,
    )
    return nv, npsi


def _shear_parameter(element):
    """Shear parameter of a beam element; zero for elements without one (trusses, nonlinear beams)."""
    get = getattr(element, "shear_parameter", None)
    return get() if get is not None else 0.0


def _linear_load_equivalent(element, wx, wy, wx_end, wy_end, start, end):
    """Work-equivalent local nodal loads [u_i, v_i, theta_i, u_j, v_j, theta_j] of a linearly varying axial and
    transverse load on the relative span [start, end]. The transverse shapes are cubic and the load linear, so
    3-point Gauss quadrature is exact."""
    L = element.length
    phi = _shear_parameter(element)
    eq_local = np.zeros(6, dtype=float)
    a = start * L
    b = end * L
    half_len = (b - a) / 2.0
    mid = (a + b) / 2.0
    for xi, w in zip(_GAUSS_POINTS, _GAUSS_WEIGHTS):
        x = mid + half_len * xi
        alpha = (xi + 1.0) / 2.0
        cur_wx = wx + alpha * (wx_end - wx)
        cur_wy = wy + alpha * (wy_end - wy)
        nv, _ = bending_shape_functions(x, L, phi)
        wt = w * half_len
        eq_local[0] += wt * (1.0 - x / L) * cur_wx
        eq_local[3] += wt * (x / L) * cur_wx
        for k, dof in enumerate((1, 2, 4, 5)):
            eq_local[dof] += wt * nv[k] * cur_wy
    return eq_local


class PointLoad:
    """
    Represents a concentrated force/moment load acting on a node.

    Attributes
    ----------
    node : Node
        The Node object to which the load is applied.
    fx : float
        Force along the global x-axis.
    fy : float
        Force along the global y-axis.
    mz : float
        Moment about the global z-axis.
    """

    def __init__(self, node, fx=0.0, fy=0.0, mz=0.0):
        """
        Initialize a PointLoad.

        Parameters
        ----------
        node : Node
            The Node object to which the load is applied.
        fx : float, optional
            Force along the global x-axis. Defaults to 0.0.
        fy : float, optional
            Force along the global y-axis. Defaults to 0.0.
        mz : float, optional
            Moment about the global z-axis. Defaults to 0.0.
        """
        self.node = node
        self.fx, self.fy, self.mz = fx, fy, mz


class ElementLoad:
    """
    Base class for loads acting along the length of an element.

    Attributes
    ----------
    element : ElementBase
        The element to which the load is applied.
    """

    def __init__(self, element):
        """
        Initialize an ElementLoad.

        Parameters
        ----------
        element : ElementBase
            The element to which the load is applied.
        """
        self.element = element

    def _compute_equivalent_loads(self):
        """
        Compute the equivalent local nodal forces and moments for this load type.

        Raises
        ------
        NotImplementedError
            If not implemented by a subclass.
        """
        raise NotImplementedError("Subclasses must implement this method")

    def _transform_and_store_equivalent_loads(self, eq_local):
        """
        Transform local equivalent nodal loads to global coordinates and
        store them in the element.

        Parameters
        ----------
        eq_local : numpy.ndarray
            6-component vector of equivalent loads in the element's local system.
        """
        if hasattr(self.element, "add_equivalent_load"):
            self.element.add_equivalent_load(eq_local)
            return

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
    Represents a uniformly or linearly varying distributed load acting on an element
    along a full or partial relative span [start, end].
    The load is defined in the element's local coordinate system:
        wx: force per unit length along local x (axial) at start
        wy: force per unit length along local y (transverse) at start
        wx_end: axial force intensity per unit length at end (defaults to wx)
        wy_end: transverse force intensity per unit length at end (defaults to wy)
        start: relative start position [0.0, 1.0] (defaults to 0.0)
        end: relative end position [0.0, 1.0] (defaults to 1.0)
    """

    def __init__(
        self,
        element,
        wx: float = 0.0,
        wy: float = 0.0,
        wx_end: float = None,
        wy_end: float = None,
        start: float = 0.0,
        end: float = 1.0,
    ):
        """
        Initialize a DistributedLoad.

        Parameters
        ----------
        element : ElementBase
            The element to which the load is applied.
        wx : float, optional
            Axial force intensity per unit length at start. Defaults to 0.0.
        wy : float, optional
            Transverse force intensity per unit length at start. Defaults to 0.0.
        wx_end : float, optional
            Axial force intensity at end. Defaults to wx.
        wy_end : float, optional
            Transverse force intensity at end. Defaults to wy.
        start : float, optional
            Relative start position along span [0, 1]. Defaults to 0.0.
        end : float, optional
            Relative end position along span [0, 1]. Defaults to 1.0.
        """
        super().__init__(element)
        if start < 0.0 or end > 1.0 or start >= end:
            raise ValueError(f"DistributedLoad relative span must satisfy 0 <= start < end <= 1, got start={start}, end={end}")

        self.wx = float(wx)
        self.wy = float(wy)
        self.wx_end = float(wx if wx_end is None else wx_end)
        self.wy_end = float(wy if wy_end is None else wy_end)
        self.start = float(start)
        self.end = float(end)
        self._compute_equivalent_loads()

    def _compute_equivalent_loads(self):
        """Work-equivalent nodal loads with the element's exact (Timoshenko or Hermite) shapes."""
        eq_local = _linear_load_equivalent(
            self.element, self.wx, self.wy, self.wx_end, self.wy_end, self.start, self.end
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
        """
        Initialize an ElementPointLoad.

        Parameters
        ----------
        element : ElementBase
            The element to which the load is applied.
        px : float, optional
            Local axial force intensity. Defaults to 0.0.
        py : float, optional
            Local transverse force intensity. Defaults to 0.0.
        mz : float, optional
            Local moment. Defaults to 0.0.
        x : float, optional
            Distance from the start node (node_i) along the element length. Defaults to 0.0.

        Raises
        ------
        ValueError
            If the distance `x` is out of bounds (less than 0 or greater than element length).
        """
        super().__init__(element)
        self.px = px
        self.py = py
        self.mz = mz
        if x > element.length or x < 0:
            raise ValueError("Load position must be within the element length.")
        self.x = x
        self._compute_equivalent_loads()

    def _compute_equivalent_loads(self):
        """Work-equivalent nodal loads in local coordinates: the force does work on the transverse shapes, the
        moment on the cross-section rotations of :func:`bending_shape_functions` (exact for Timoshenko beams)."""
        L = self.element.length
        a = self.x
        b = L - a
        nv, npsi = bending_shape_functions(a, L, _shear_parameter(self.element))
        eq_local = np.zeros(6, dtype=float)
        eq_local[0] = (self.px * b) / L
        eq_local[3] = (self.px * a) / L
        for k, dof in enumerate((1, 2, 4, 5)):
            eq_local[dof] = self.py * nv[k] + self.mz * npsi[k]
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
        """
        Initialize a TriangularLoad/TrapezoidalLoad.

        Parameters
        ----------
        element : ElementBase
            The element to which the load is applied.
        wx1 : float, optional
            Axial force intensity at the start node. Defaults to 0.0.
        wx2 : float, optional
            Axial force intensity at the end node. Defaults to 0.0.
        wy1 : float, optional
            Transverse force intensity at the start node. Defaults to 0.0.
        wy2 : float, optional
            Transverse force intensity at the end node. Defaults to 0.0.
        """
        super().__init__(element)
        self.wx1 = wx1
        self.wx2 = wx2
        self.wy1 = wy1
        self.wy2 = wy2
        self._compute_equivalent_loads()

    def _compute_equivalent_loads(self):
        """Work-equivalent nodal loads with the element's exact (Timoshenko or Hermite) shapes, over the full span."""
        eq_local = _linear_load_equivalent(self.element, self.wx1, self.wy1, self.wx2, self.wy2, 0.0, 1.0)
        self._transform_and_store_equivalent_loads(eq_local)
