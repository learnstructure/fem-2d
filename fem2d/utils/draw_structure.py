import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Polygon, Circle, FancyArrowPatch, Arc
from matplotlib.collections import LineCollection


class DrawStructure:
    def __init__(self, structure, scale=1.0, arrow_scale=0.1):
        self.structure = structure
        self.scale = scale
        self.arrow_scale = arrow_scale  # fraction of structure span

    def _get_span(self):
        """Return characteristic length for scaling symbols."""
        x = [node.x for node in self.structure.nodes.values()]
        y = [node.y for node in self.structure.nodes.values()]
        span = max(max(x) - min(x), max(y) - min(y))
        if span == 0:
            span = 1.0
        return span

    def _draw_support(self, node, support_type, span):
        """Draw support symbol at node based on fixity."""
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
        """Draw point loads, moments, and distributed loads."""
        # Determine max force for scaling arrows
        max_force = 0
        for node in self.structure.nodes.values():
            fx, fy, mz = node.load
            max_force = max(max_force, abs(fx), abs(fy))
        # For UDLs, we need to iterate over loads in structure.loads
        # Assuming UDLs are stored in structure.loads as DistributedLoad objects
        for load in self.structure.loads:
            if hasattr(load, "wy") and load.wy != 0:  # vertical UDL
                el = load.element
                # Draw small arrows along element
                n_arrows = 11
                x_i, y_i = el.node_i.x, el.node_i.y
                x_j, y_j = el.node_j.x, el.node_j.y
                for t in np.linspace(0, 1, n_arrows):
                    x = x_i + t * (x_j - x_i)
                    y = y_i + t * (y_j - y_i)
                    # Arrow direction depends on load sign (positive wy = upward)
                    dy = 0.05 * span * (1 if load.wy > 0 else -1)
                    plt.arrow(
                        x,
                        y,
                        0,
                        dy,
                        head_width=0.02 * span,
                        head_length=0.02 * span,
                        fc="green",
                        ec="green",
                        alpha=0.5,
                    )

        # Point loads at nodes
        arrow_len = self.arrow_scale * span
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
                arrow = FancyArrowPatch(
                    (node.x, node.y),
                    (node.x + dx, node.y + dy),
                    arrowstyle="->",
                    color="red",
                    linewidth=2,
                )
                plt.gca().add_patch(arrow)

            # Moment
            if mz != 0:
                # Curved arrow radius
                r = 0.05 * span
                # Determine direction: positive -> counter‑clockwise
                if mz > 0:
                    start_angle, end_angle = 0, 90
                else:
                    start_angle, end_angle = 90, 0
                arc = Arc(
                    (node.x, node.y),
                    2 * r,
                    2 * r,
                    theta1=start_angle,
                    theta2=end_angle,
                    edgecolor="red",
                    linewidth=2,
                )
                plt.gca().add_patch(arc)
                # Add arrowhead at the end
                angle = np.radians(end_angle)
                tip_x = node.x + r * np.cos(angle)
                tip_y = node.y + r * np.sin(angle)
                # Small triangle
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
