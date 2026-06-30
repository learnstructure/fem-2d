"""
DrawStructure module implementing visualization of 2D frame structures using Matplotlib.
Supports drawing support symbols, point loads, moments, distributed loads, and deformed shapes.
"""

import numpy as np

try:
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon, Circle, FancyArrowPatch, Arc
    from matplotlib.collections import LineCollection
except ImportError:  # pragma: no cover - optional dependency for docs builds
    plt = None
    Polygon = Circle = FancyArrowPatch = Arc = LineCollection = None


class DrawStructure:
    """
    Visualizes the structural model geometry, boundary conditions, applied loads, and deformed shape.

    Attributes
    ----------
    structure : Structure
        The Structure object to draw.
    scale : float
        Displacement magnification factor for drawing deformed shapes.
    arrow_scale : float
        Scaling factor for loading arrows relative to the characteristic span of the structure.
    """

    def __init__(self, structure, scale=1.0, arrow_scale=0.1):
        """
        Initialize the structure plotter.

        Parameters
        ----------
        structure : Structure
            The Structure object to draw.
        scale : float, optional
            Displacement magnification factor. Defaults to 1.0.
        arrow_scale : float, optional
            Loading arrows scale factor. Defaults to 0.1.
        """
        self.structure = structure
        self.scale = scale
        self.arrow_scale = arrow_scale  # fraction of structure span

    def _get_span(self):
        """
        Return the characteristic length (span) of the structure for scaling symbols.

        Returns
        -------
        float
            The characteristic length.
        """
        x = [node.x for node in self.structure.nodes.values()]
        y = [node.y for node in self.structure.nodes.values()]
        span = max(max(x) - min(x), max(y) - min(y))
        if span == 0:
            span = 1.0
        return span

    def _draw_support(self, node, support_type, span):
        """
        Draw boundary condition (support) symbol at a node based on fixity.

        Parameters
        ----------
        node : Node
            The Node object where support is located.
        support_type : list of bool
            Fixity condition list [ux, uy, rz].
        span : float
            Characteristic structure span for scaling the support symbol.
        """
        # Support size relative to span
        size = 0.03 * span
        x, y = node.x, node.y

        # Determine support type from node.support list
        ux, uy, rz = node.support

        # Fixed (all DOFs fixed) → filled triangle at base
        if ux and uy and rz:
            # Triangle pointing downward from node
            triangle = Polygon(
                [(x - size, y), (x + size, y), (x, y - 1.5 * size)],
                closed=True,
                facecolor="lightgray",
                edgecolor="black",
                linewidth=1,
            )
            plt.gca().add_patch(triangle)

        # Pinned (ux, uy fixed, rz free) → triangle + circle
        elif ux and uy and not rz:
            # Triangle
            triangle = Polygon(
                [(x - size, y), (x + size, y), (x, y - 1.5 * size)],
                closed=True,
                facecolor="lightgray",
                edgecolor="black",
                linewidth=1,
            )
            plt.gca().add_patch(triangle)
            # Circle at top
            circle = Circle(
                (x, y - 0.1 * size),
                radius=0.2 * size,
                facecolor="white",
                edgecolor="black",
                linewidth=1,
            )
            plt.gca().add_patch(circle)

        # Roller (only one translational DOF fixed, e.g., uy fixed) → two small circles under node
        else:
            # Two wheels (circles)
            wheel_radius = 0.15 * size
            circle1 = Circle(
                (x - 0.3 * size, y - wheel_radius),
                radius=wheel_radius,
                facecolor="gray",
                edgecolor="black",
            )
            circle2 = Circle(
                (x + 0.3 * size, y - wheel_radius),
                radius=wheel_radius,
                facecolor="gray",
                edgecolor="black",
            )
            plt.gca().add_patch(circle1)
            plt.gca().add_patch(circle2)
            # Small triangle to connect to node
            plt.plot([x, x - 0.3 * size], [y, y - wheel_radius], "k-", linewidth=1)
            plt.plot([x, x + 0.3 * size], [y, y - wheel_radius], "k-", linewidth=1)

    def _draw_loads(self, span):
        """
        Draw point loads, moments, and distributed loads acting on the structure.

        Parameters
        ----------
        span : float
            Characteristic structure span for scaling the load arrows and symbols.
        """
        from fem2d.loads import PointLoad, ElementPointLoad, DistributedLoad
        
        # Determine max force for scaling arrows
        max_force = 0
        
        # Check node loads and structure loads
        for node in self.structure.nodes.values():
            fx, fy, mz = node.load
            max_force = max(max_force, abs(fx), abs(fy))
        
        for load in self.structure.loads:
            if isinstance(load, PointLoad):
                max_force = max(max_force, abs(load.fx), abs(load.fy))
            elif isinstance(load, ElementPointLoad):
                max_force = max(max_force, abs(load.px), abs(load.py))
            elif isinstance(load, DistributedLoad):
                el = load.element
                if load.wx != 0:
                    max_force = max(max_force, abs(load.wx * el.length))
                if load.wy != 0:
                    max_force = max(max_force, abs(load.wy * el.length))
        
        # Draw distributed loads (UDLs) - transform from local to global coordinates
        for load in self.structure.loads:
            if isinstance(load, DistributedLoad):
                el = load.element
                c, s = el.cos, el.sin
                n_arrows = 10  # Reduced from 11 for better clarity
                x_i, y_i = el.node_i.x, el.node_i.y
                x_j, y_j = el.node_j.x, el.node_j.y
                
                # Collect arrow tip positions for drawing connecting line
                arrow_tips_x = []
                arrow_tips_y = []
                arrow_head_extension = 0.03 * span  # Extension to account for arrow head
                
                for t in np.linspace(0, 1, n_arrows):
                    # Position along element
                    x = x_i + t * (x_j - x_i)
                    y = y_i + t * (y_j - y_i)
                    
                    # Transform local load (wx, wy) to global (fx, gy)
                    # Global: fx = wx*c - wy*s,  fy = wx*s + wy*c
                    fx_local = load.wx
                    fy_local = load.wy
                    fx_global = fx_local * c - fy_local * s
                    fy_global = fx_local * s + fy_local * c
                    
                    # Compute arrow length
                    magnitude = np.hypot(fx_global, fy_global)
                    if magnitude > 0:
                        length = self.arrow_scale * span * 0.5  # Adjusted scale
                        dx = fx_global / magnitude * length
                        dy = fy_global / magnitude * length
                    else:
                        dx = dy = 0
                    
                    if magnitude > 0:
                        # Account for arrow head extending beyond endpoint
                        direction_x = fx_global / magnitude
                        direction_y = fy_global / magnitude
                        # Visual tip includes the arrow head extension
                        arrow_tip_x = x + dx + direction_x * arrow_head_extension
                        arrow_tip_y = y + dy + direction_y * arrow_head_extension
                        arrow_tips_x.append(arrow_tip_x)
                        arrow_tips_y.append(arrow_tip_y)
                        
                        plt.arrow(
                            x, y, dx, dy,
                            head_width=0.025 * span,
                            head_length=0.025 * span,
                            fc="green",
                            ec="green",
                            alpha=0.8,
                            linewidth=1.0
                        )
                
                # Draw connecting line through arrow tips to show distributed load
                if len(arrow_tips_x) > 1:
                    plt.plot(arrow_tips_x, arrow_tips_y, 'g--', linewidth=1.0, alpha=0.8)
        
        # Draw point loads at nodes (from structure.loads PointLoad objects)
        arrow_len = self.arrow_scale * span
        
        for load in self.structure.loads:
            if isinstance(load, PointLoad):
                # PointLoad is already in global coordinates
                fx, fy, mz = load.fx, load.fy, load.mz
                node = load.node
                
                # Draw force
                if fx != 0 or fy != 0:
                    magnitude = np.hypot(fx, fy)
                    min_length = 1 * arrow_len  # Minimum arrow length for visibility
                    length = max(
                        min_length,
                        arrow_len * (magnitude / max_force) if max_force != 0 else arrow_len
                    )
                    dx = fx / magnitude * length
                    dy = fy / magnitude * length
                    
                    # Draw arrow (includes line + head)
                    arrow = FancyArrowPatch(
                        (node.x, node.y),
                        (node.x + dx, node.y + dy),
                        arrowstyle="->",
                        color="red",
                        linewidth=3,
                        mutation_scale=25
                    )
                    plt.gca().add_patch(arrow)
                
                # Draw moment
                if mz != 0:
                    r = 0.05 * span
                    if mz > 0:
                        start_angle, end_angle = 0, 90
                    else:
                        start_angle, end_angle = 90, 0
                    arc = Arc(
                        (node.x, node.y),
                        2 * r, 2 * r,
                        theta1=start_angle,
                        theta2=end_angle,
                        edgecolor="red",
                        linewidth=2
                    )
                    plt.gca().add_patch(arc)
                    angle = np.radians(end_angle)
                    tip_x = node.x + r * np.cos(angle)
                    tip_y = node.y + r * np.sin(angle)
                    arrowhead = Polygon(
                        [
                            (tip_x, tip_y),
                            (tip_x - 0.01 * span, tip_y - 0.01 * span),
                            (tip_x + 0.01 * span, tip_y - 0.01 * span),
                        ],
                        closed=True,
                        facecolor="red",
                    )
                    plt.gca().add_patch(arrowhead)
        
        # Draw element point loads - transform from local to global coordinates
        for load in self.structure.loads:
            if isinstance(load, ElementPointLoad):
                el = load.element
                c, s = el.cos, el.sin
                
                # Position of load on element
                t = load.x / el.length
                x = el.node_i.x + t * (el.node_j.x - el.node_i.x)
                y = el.node_i.y + t * (el.node_j.y - el.node_i.y)
                
                # Transform local load (px, py) to global (fx, fy)
                px_local = load.px
                py_local = load.py
                fx_global = px_local * c - py_local * s
                fy_global = px_local * s + py_local * c
                
                # Draw force
                if fx_global != 0 or fy_global != 0:
                    magnitude = np.hypot(fx_global, fy_global)
                    length = (
                        arrow_len * (magnitude / max_force) if max_force != 0 else arrow_len
                    )
                    dx = fx_global / magnitude * length
                    dy = fy_global / magnitude * length
                    
                    # Draw arrow (includes line + head)
                    arrow = FancyArrowPatch(
                        (x, y),
                        (x + dx, y + dy),
                        arrowstyle="->",
                        color="orange",
                        linewidth=3,
                        mutation_scale=25
                    )
                    plt.gca().add_patch(arrow)
                
                # Draw moment (if any)
                if load.mz != 0:
                    r = 0.05 * span
                    if load.mz > 0:
                        start_angle, end_angle = 0, 90
                    else:
                        start_angle, end_angle = 90, 0
                    arc = Arc(
                        (x, y),
                        2 * r, 2 * r,
                        theta1=start_angle,
                        theta2=end_angle,
                        edgecolor="orange",
                        linewidth=2
                    )
                    plt.gca().add_patch(arc)
                    angle = np.radians(end_angle)
                    tip_x = x + r * np.cos(angle)
                    tip_y = y + r * np.sin(angle)
                    arrowhead = Polygon(
                        [
                            (tip_x, tip_y),
                            (tip_x - 0.01 * span, tip_y - 0.01 * span),
                            (tip_x + 0.01 * span, tip_y - 0.01 * span),
                        ],
                        closed=True,
                        facecolor="orange",
                    )
                    plt.gca().add_patch(arrowhead)
        
        # Draw node loads from node.load attribute
        for node in self.structure.nodes.values():
            fx, fy, mz = node.load
            
            # Force
            if fx != 0 or fy != 0:
                magnitude = np.hypot(fx, fy)
                length = (
                    arrow_len * (magnitude / max_force) if max_force != 0 else arrow_len
                )
                dx = fx / magnitude * length
                dy = fy / magnitude * length
                
                # Draw arrow (includes line + head)
                arrow = FancyArrowPatch(
                    (node.x, node.y),
                    (node.x + dx, node.y + dy),
                    arrowstyle="->",
                    color="red",
                    linewidth=3,
                    mutation_scale=25
                )
                plt.gca().add_patch(arrow)
            
            # Moment
            if mz != 0:
                r = 0.05 * span
                if mz > 0:
                    start_angle, end_angle = 0, 90
                else:
                    start_angle, end_angle = 90, 0
                arc = Arc(
                    (node.x, node.y),
                    2 * r, 2 * r,
                    theta1=start_angle,
                    theta2=end_angle,
                    edgecolor="red",
                    linewidth=2
                )
                plt.gca().add_patch(arc)
                angle = np.radians(end_angle)
                tip_x = node.x + r * np.cos(angle)
                tip_y = node.y + r * np.sin(angle)
                arrowhead = Polygon(
                    [
                        (tip_x, tip_y),
                        (tip_x - 0.01 * span, tip_y - 0.01 * span),
                        (tip_x + 0.01 * span, tip_y - 0.01 * span),
                    ],
                    closed=True,
                    facecolor="red",
                )
                plt.gca().add_patch(arrowhead)

    def draw(self, show_undeformed=True, show_deformed=True, n_points=20):
        """
        Draw the structure plot (geometry, loads, boundary conditions, displacements) using Matplotlib.

        Parameters
        ----------
        show_undeformed : bool, optional
            Whether to draw the undeformed structure. Defaults to True.
        show_deformed : bool, optional
            Whether to draw the deformed structure. Defaults to True.
        n_points : int, optional
            Number of points for plotting deformed shape curves. Defaults to 20.
        """
        if plt is None:
            raise ImportError("matplotlib is required to draw structures")

        plt.figure(figsize=(10, 8))
        ax = plt.gca()
        span = self._get_span()

        # Undeformed (blue)
        if show_undeformed:
            for el in self.structure.elements.values():
                x = [el.node_i.x, el.node_j.x]
                y = [el.node_i.y, el.node_j.y]
                plt.plot(
                    x, y, "b-", linewidth=1, label="Undeformed" if el.id == 1 else ""
                )

        # Deformed (red)
        if show_deformed and self.structure.disp is not None:
            for el in self.structure.elements.values():
                points = el.deformed_shape_points(
                    self.structure.disp, n_points, self.scale
                )
                x_def = [p[0] for p in points]
                y_def = [p[1] for p in points]
                plt.plot(
                    x_def,
                    y_def,
                    "r-",
                    linewidth=1.5,
                    label="Deformed" if el.id == 1 else "",
                )

        # Nodes
        for node in self.structure.nodes.values():
            plt.plot(node.x, node.y, "ko", markersize=4)
            plt.text(node.x, node.y, f" {node.id}", ha="left", va="bottom", fontsize=8)

        # Supports
        for node in self.structure.nodes.values():
            if any(node.support):
                self._draw_support(node, node.support, span)

        # Loads
        self._draw_loads(span)

        # Legend
        if show_undeformed and show_deformed:
            plt.legend()

        plt.xlabel("x")
        plt.ylabel("y")
        plt.title("Structure: Undeformed (blue) and Deformed (red)")
        plt.axis("equal")
        plt.grid(True)
        plt.show()
