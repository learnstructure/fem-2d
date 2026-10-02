"""
DrawStructure module implementing visualization of 2D frame structures using Matplotlib.
Supports drawing realistic support symbols (fixed, pin, roller, box), point loads, moments,
distributed loads, member hinges, and deformed shapes.
"""

import numpy as np

try:
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon, Circle, FancyArrowPatch, Arc, Rectangle
    from matplotlib.collections import LineCollection
except ImportError:  # pragma: no cover - optional dependency for docs builds
    plt = None
    Polygon = Circle = FancyArrowPatch = Arc = Rectangle = LineCollection = None


class DrawStructure:
    """
    Visualizes structural model geometry, boundary conditions, applied loads,
    internal releases, and deformed shapes using professional engineering aesthetics.

    Attributes
    ----------
    structure : Structure
        The Structure object to draw.
    scale : float
        Displacement magnification factor for drawing deformed shapes.
    arrow_scale : float
        Scaling factor for loading arrows relative to the characteristic span of the structure.
    """

    def __init__(self, structure, scale=1.0, arrow_scale=0.10):
        """
        Initialize the structure plotter.

        Parameters
        ----------
        structure : Structure
            The Structure object to draw.
        scale : float, optional
            Displacement magnification factor. Defaults to 1.0.
        arrow_scale : float, optional
            Loading arrows scale factor relative to structure span. Defaults to 0.10.
        """
        self.structure = structure
        self.scale = scale
        self.arrow_scale = arrow_scale

    def _get_span(self):
        """
        Return the characteristic length (span) of the structure for scaling symbols.

        Returns
        -------
        float
            The characteristic length.
        """
        if not self.structure.nodes:
            return 1.0
        x = [node.x for node in self.structure.nodes.values()]
        y = [node.y for node in self.structure.nodes.values()]
        dx = max(x) - min(x)
        dy = max(y) - min(y)
        span = max(dx, dy)
        return span if span > 1e-6 else 1.0

    def _get_support_normal(self, node):
        """
        Determine the outward normal vector pointing away from the structure into the support.
        This enables intelligent orientation (e.g., ground beneath vertical columns,
        wall behind horizontal cantilever ends, ground beneath continuous beams).

        Parameters
        ----------
        node : Node
            The node being inspected.

        Returns
        -------
        tuple of float
            Unit normal vector (nx, ny) pointing into the support/foundation.
        """
        vx, vy = 0.0, 0.0
        count = 0
        purely_horizontal = True
        purely_vertical = True

        for el in self.structure.elements.values():
            if el.node_i.id == node.id:
                dx = el.node_j.x - node.x
                dy = el.node_j.y - node.y
                length = np.hypot(dx, dy)
                if length > 1e-9:
                    vx += dx / length
                    vy += dy / length
                    count += 1
                    if abs(dy) > 1e-4:
                        purely_horizontal = False
                    if abs(dx) > 1e-4:
                        purely_vertical = False
            elif el.node_j.id == node.id:
                dx = el.node_i.x - node.x
                dy = el.node_i.y - node.y
                length = np.hypot(dx, dy)
                if length > 1e-9:
                    vx += dx / length
                    vy += dy / length
                    count += 1
                    if abs(dy) > 1e-4:
                        purely_horizontal = False
                    if abs(dx) > 1e-4:
                        purely_vertical = False

        ux, uy, rz = node.support

        # Roller restraining only Y: reaction is vertical, ground is beneath
        if not ux and uy:
            return 0.0, -1.0

        # Roller restraining only X: reaction is horizontal
        if ux and not uy:
            return (-1.0, 0.0) if vx >= 0 else (1.0, 0.0)

        # For fixed or pinned supports on purely horizontal members (e.g. continuous beams)
        if purely_horizontal and count > 0:
            # Fixed cantilever attached to wall at left/right end
            if ux and uy and rz and count == 1:
                return (-1.0, 0.0) if vx > 0 else (1.0, 0.0)
            # Continuous beam internal or pinned support: ground underneath
            return 0.0, -1.0

        # For members with vertical extent: ground is opposite to members
        if count > 0:
            fx = -vx
            fy = -vy
            if abs(fy) >= abs(fx):
                return (0.0, -1.0) if fy <= 0 else (0.0, 1.0)
            else:
                return (-1.0, 0.0) if fx <= 0 else (1.0, 0.0)

        return 0.0, -1.0

    def _draw_hatch_lines(self, ax, p_start, p_end, normal, tangent, hatch_len=0.015, n_hatches=6):
        """Draw 45-degree slanted ground hatch lines along a baseline."""
        p_start = np.array(p_start, dtype=float)
        p_end = np.array(p_end, dtype=float)
        normal = np.array(normal, dtype=float)
        tangent = np.array(tangent, dtype=float)

        # Vector pointing at 45 degrees into ground
        hatch_dir = (normal + tangent) / np.sqrt(2.0) * hatch_len

        for t in np.linspace(0.08, 0.92, n_hatches):
            base_pt = p_start + t * (p_end - p_start)
            tip_pt = base_pt + hatch_dir
            ax.plot([base_pt[0], tip_pt[0]], [base_pt[1], tip_pt[1]], color="#555555", linewidth=0.9, zorder=2)

    def _draw_fixed_support(self, ax, x, y, normal, tangent, size):
        """Draw a realistic fixed/clamped support (base plate + ground hatching)."""
        nx, ny = normal
        tx, ty = tangent

        width = 1.8 * size
        thickness = 0.35 * size

        center = np.array([x, y])
        norm_v = np.array([nx, ny])
        tang_v = np.array([tx, ty])

        # Base plate corners
        p1 = center - 0.5 * width * tang_v
        p2 = center + 0.5 * width * tang_v
        p3 = p2 + thickness * norm_v
        p4 = p1 + thickness * norm_v

        plate = Polygon(
            [p1, p2, p3, p4],
            closed=True,
            facecolor="#4A5568",
            edgecolor="#1A202C",
            linewidth=1.2,
            zorder=3,
        )
        ax.add_patch(plate)

        # Ground hatch lines behind the plate
        self._draw_hatch_lines(ax, p4, p3, norm_v, tang_v, hatch_len=0.5 * size, n_hatches=6)

    def _draw_pin_support(self, ax, x, y, normal, tangent, size):
        """Draw a classic pinned support (triangle + pivot circle + hatched ground line)."""
        nx, ny = normal
        tx, ty = tangent

        height = 1.0 * size
        base_width = 1.2 * size

        center = np.array([x, y])
        norm_v = np.array([nx, ny])
        tang_v = np.array([tx, ty])

        b_center = center + height * norm_v
        b1 = b_center - 0.5 * base_width * tang_v
        b2 = b_center + 0.5 * base_width * tang_v

        # Pinned triangle
        triangle = Polygon(
            [center, b2, b1],
            closed=True,
            facecolor="#EDF2F7",
            edgecolor="#2D3748",
            linewidth=1.2,
            zorder=3,
        )
        ax.add_patch(triangle)

        # Ground line extending slightly beyond triangle base
        ground_ext = 0.2 * base_width
        g1 = b1 - ground_ext * tang_v
        g2 = b2 + ground_ext * tang_v
        ax.plot([g1[0], g2[0]], [g1[1], g2[1]], color="#2D3748", linewidth=1.2, zorder=3)

        # Ground hatch lines
        self._draw_hatch_lines(ax, g1, g2, norm_v, tang_v, hatch_len=0.45 * size, n_hatches=5)

        # Pin circle at apex
        pin_radius = 0.16 * size
        circle = Circle(
            (x, y),
            radius=pin_radius,
            facecolor="white",
            edgecolor="#2D3748",
            linewidth=1.2,
            zorder=5,
        )
        ax.add_patch(circle)

    def _draw_roller_support(self, ax, x, y, normal, tangent, size):
        """Draw a roller support (triangle + wheels + track line + ground hatching)."""
        nx, ny = normal
        tx, ty = tangent

        height = 0.75 * size
        base_width = 1.05 * size

        center = np.array([x, y])
        norm_v = np.array([nx, ny])
        tang_v = np.array([tx, ty])

        b_center = center + height * norm_v
        b1 = b_center - 0.5 * base_width * tang_v
        b2 = b_center + 0.5 * base_width * tang_v

        # Triangle body
        triangle = Polygon(
            [center, b2, b1],
            closed=True,
            facecolor="#EDF2F7",
            edgecolor="#2D3748",
            linewidth=1.2,
            zorder=3,
        )
        ax.add_patch(triangle)

        # Pin circle at apex
        pin_radius = 0.16 * size
        circle = Circle(
            (x, y),
            radius=pin_radius,
            facecolor="white",
            edgecolor="#2D3748",
            linewidth=1.2,
            zorder=5,
        )
        ax.add_patch(circle)

        # Dual roller wheels
        wheel_r = 0.14 * size
        w1_center = b_center - 0.28 * base_width * tang_v + wheel_r * norm_v
        w2_center = b_center + 0.28 * base_width * tang_v + wheel_r * norm_v

        wheel1 = Circle(w1_center, radius=wheel_r, facecolor="white", edgecolor="#2D3748", linewidth=1.0, zorder=3)
        wheel2 = Circle(w2_center, radius=wheel_r, facecolor="white", edgecolor="#2D3748", linewidth=1.0, zorder=3)
        ax.add_patch(wheel1)
        ax.add_patch(wheel2)

        # Track line beneath wheels
        track_center = b_center + 2.0 * wheel_r * norm_v
        track_ext = 0.65 * base_width
        t1 = track_center - track_ext * tang_v
        t2 = track_center + track_ext * tang_v
        ax.plot([t1[0], t2[0]], [t1[1], t2[1]], color="#2D3748", linewidth=1.2, zorder=3)

        # Hatch lines under track
        self._draw_hatch_lines(ax, t1, t2, norm_v, tang_v, hatch_len=0.40 * size, n_hatches=5)

    def _draw_box_support(self, ax, x, y, normal, tangent, size, is_pin=False):
        """Draw a sleek rectangular shaded pad support."""
        nx, ny = normal
        tx, ty = tangent

        width = 1.6 * size
        depth = 0.5 * size

        center = np.array([x, y])
        norm_v = np.array([nx, ny])
        tang_v = np.array([tx, ty])

        p1 = center - 0.5 * width * tang_v
        p2 = center + 0.5 * width * tang_v
        p3 = p2 + depth * norm_v
        p4 = p1 + depth * norm_v

        box = Polygon(
            [p1, p2, p3, p4],
            closed=True,
            facecolor="#CBD5E0",
            edgecolor="#2D3748",
            linewidth=1.2,
            zorder=3,
        )
        ax.add_patch(box)

        # Hatch lines
        self._draw_hatch_lines(ax, p4, p3, norm_v, tang_v, hatch_len=0.45 * size, n_hatches=5)

        if is_pin:
            circle = Circle((x, y), radius=0.15 * size, facecolor="white", edgecolor="#2D3748", linewidth=1.2, zorder=5)
            ax.add_patch(circle)

    def _draw_support(self, ax, node, span, support_style="detailed"):
        """
        Draw boundary condition symbol at a node based on fixity.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Target plot axes.
        node : Node
            The Node object where support is located.
        span : float
            Characteristic structure span for scaling the support symbol.
        support_style : str, optional
            'detailed' for distinct fixed/pin/roller symbols, or 'box' for shaded box pads.
        """
        size = 0.035 * span
        x, y = node.x, node.y
        ux, uy, rz = node.support

        if not any([ux, uy, rz]):
            return

        nx, ny = self._get_support_normal(node)
        tx, ty = -ny, nx

        if support_style == "box":
            is_pin = (ux and uy and not rz) or (not ux and uy and not rz)
            self._draw_box_support(ax, x, y, (nx, ny), (tx, ty), size, is_pin=is_pin)
            return

        # Fixed Support (all 3 DOFs restrained)
        if ux and uy and rz:
            self._draw_fixed_support(ax, x, y, (nx, ny), (tx, ty), size)

        # Pinned Support (translation restrained, rotation free)
        elif ux and uy and not rz:
            self._draw_pin_support(ax, x, y, (nx, ny), (tx, ty), size)

        # Roller or single restrained translation
        else:
            self._draw_roller_support(ax, x, y, (nx, ny), (tx, ty), size)

    def _draw_point_arrow(self, ax, x, y, fx, fy, max_force, span, label_val=None, color="#C53030"):
        """
        Draw a concentrated force arrow pointing directly towards its application point (x, y).
        Includes a clean magnitude annotation badge placed outside the tail.
        """
        magnitude = np.hypot(fx, fy)
        if magnitude < 1e-12:
            return

        arrow_base_len = self.arrow_scale * span
        length = arrow_base_len * (0.65 + 0.35 * (magnitude / max_force))
        length = max(length, 0.06 * span)

        dir_x = fx / magnitude
        dir_y = fy / magnitude

        # Arrow points TOWARDS the application point (tail -> tip at node)
        tip_x, tip_y = x, y
        tail_x = tip_x - dir_x * length
        tail_y = tip_y - dir_y * length

        arrow = FancyArrowPatch(
            (tail_x, tail_y),
            (tip_x, tip_y),
            arrowstyle="-|>",
            mutation_scale=16,
            color=color,
            linewidth=2.0,
            zorder=6,
        )
        ax.add_patch(arrow)

        # Value annotation near tail
        if label_val is not None:
            text_x = tail_x - dir_x * (0.015 * span)
            text_y = tail_y - dir_y * (0.015 * span)

            if abs(dir_y) >= abs(dir_x):
                va = "bottom" if -dir_y > 0 else "top"
                ha = "center"
            else:
                ha = "left" if -dir_x > 0 else "right"
                va = "center"

            val_str = f"{abs(label_val):g}"
            ax.text(
                text_x,
                text_y,
                val_str,
                fontsize=8.5,
                fontweight="bold",
                color=color,
                ha=ha,
                va=va,
                bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="#CBD5E0", lw=0.6, alpha=0.92),
                zorder=7,
            )

    def _draw_moment_arc(self, ax, x, y, mz, span, color="#7B1FA2"):
        """Draw a curved circular arc with arrowhead for concentrated moments."""
        if abs(mz) < 1e-12:
            return

        r = 0.045 * span
        # Counterclockwise (positive) or clockwise (negative)
        if mz > 0:
            angles = np.linspace(-np.pi * 0.2, np.pi * 1.15, 30)
            tip_angle = angles[-1]
            # Tangent direction pointing CCW
            tangent_angle = tip_angle + np.pi / 2
        else:
            angles = np.linspace(np.pi * 1.2, -np.pi * 0.15, 30)
            tip_angle = angles[-1]
            # Tangent direction pointing CW
            tangent_angle = tip_angle - np.pi / 2

        arc_x = x + r * np.cos(angles)
        arc_y = y + r * np.sin(angles)
        ax.plot(arc_x, arc_y, color=color, linewidth=1.8, zorder=6)

        # Arrowhead at tip
        tip_x = x + r * np.cos(tip_angle)
        tip_y = y + r * np.sin(tip_angle)
        head_len = 0.015 * span
        head_p1 = (
            tip_x - head_len * np.cos(tangent_angle - 0.5),
            tip_y - head_len * np.sin(tangent_angle - 0.5),
        )
        head_p2 = (
            tip_x - head_len * np.cos(tangent_angle + 0.5),
            tip_y - head_len * np.sin(tangent_angle + 0.5),
        )
        head = Polygon([ (tip_x, tip_y), head_p1, head_p2 ], closed=True, facecolor=color, edgecolor=color, zorder=7)
        ax.add_patch(head)

        # Label
        label_str = f"M = {abs(mz):g}"
        offset_dir = 1.0 if mz > 0 else -1.0
        ax.text(
            x,
            y + offset_dir * 1.35 * r,
            label_str,
            fontsize=8.5,
            fontweight="bold",
            color=color,
            ha="center",
            va="bottom" if offset_dir > 0 else "top",
            bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="#CBD5E0", lw=0.6, alpha=0.92),
            zorder=7,
        )

    def _draw_distributed_load(self, ax, el, wx1, wx2, wy1, wy2, max_force, span, label_str=None, color="#2E7D32"):
        """
        Draw uniform or linearly varying distributed loads with translucent fill,
        neat boundary line, evenly spaced arrows, and clear load notations.
        """
        max_w = max(abs(wx1), abs(wx2), abs(wy1), abs(wy2))
        if max_w < 1e-12:
            return

        c, s = el.cos, el.sin
        # Element local axes: local_x = (c, s), local_y (transverse) = (-s, c)
        ex = np.array([c, s])
        ey = np.array([-s, c])

        base_height = self.arrow_scale * span * 0.65

        # Sample points along element
        n_arrows = 7
        t_vals = np.linspace(0.0, 1.0, n_arrows)

        tails_x = []
        tails_y = []
        bases_x = []
        bases_y = []

        for t in t_vals:
            # Local load intensity at t
            wx_t = (1.0 - t) * wx1 + t * wx2
            wy_t = (1.0 - t) * wy1 + t * wy2

            # Global load vector: wx*ex + wy*ey
            w_glob = wx_t * ex + wy_t * ey
            w_mag = np.hypot(w_glob[0], w_glob[1])

            # Point on beam
            px = el.node_i.x + t * (el.node_j.x - el.node_i.x)
            py = el.node_i.y + t * (el.node_j.y - el.node_i.y)
            bases_x.append(px)
            bases_y.append(py)

            if w_mag > 1e-9:
                u_x = w_glob[0] / w_mag
                u_y = w_glob[1] / w_mag
                h = base_height * (w_mag / max_w)
                h = max(h, 0.03 * span)

                # Tail is offset away; arrow points towards beam
                tx = px - u_x * h
                ty = py - u_y * h
                tails_x.append(tx)
                tails_y.append(ty)

                # Arrow pointing from tail to base on beam
                arrow = FancyArrowPatch(
                    (tx, ty),
                    (px, py),
                    arrowstyle="-|>",
                    mutation_scale=10,
                    color=color,
                    linewidth=1.1,
                    zorder=4,
                )
                ax.add_patch(arrow)
            else:
                tails_x.append(px)
                tails_y.append(py)

        # Translucent load area fill
        poly_x = bases_x + tails_x[::-1]
        poly_y = bases_y + tails_y[::-1]
        ax.fill(poly_x, poly_y, color=color, alpha=0.14, zorder=2)

        # Connecting boundary line across arrow tails
        ax.plot(tails_x, tails_y, color=color, linestyle="-", linewidth=1.1, zorder=3)

        # Load notation label centered above load boundary line
        if label_str:
            mid_idx = len(tails_x) // 2
            mid_tail_x = tails_x[mid_idx]
            mid_tail_y = tails_y[mid_idx]
            mid_base_x = bases_x[mid_idx]
            mid_base_y = bases_y[mid_idx]

            dir_norm_x = mid_tail_x - mid_base_x
            dir_norm_y = mid_tail_y - mid_base_y
            dn_len = np.hypot(dir_norm_x, dir_norm_y)
            if dn_len > 1e-9:
                lbl_x = mid_tail_x + (dir_norm_x / dn_len) * (0.015 * span)
                lbl_y = mid_tail_y + (dir_norm_y / dn_len) * (0.015 * span)
            else:
                lbl_x = mid_tail_x
                lbl_y = mid_tail_y + 0.02 * span

            ax.text(
                lbl_x,
                lbl_y,
                label_str,
                fontsize=8.5,
                fontweight="bold",
                color=color,
                ha="center",
                va="center",
                bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="#A5D6A7", lw=0.6, alpha=0.92),
                zorder=7,
            )

    def _draw_loads(self, ax, span):
        """
        Draw point loads, moments, and distributed loads acting on the structure.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Target plot axes.
        span : float
            Characteristic structure span for scaling the load arrows and symbols.
        """
        from fem2d.loads import PointLoad, ElementPointLoad, DistributedLoad, TriangularLoad

        # Determine max force for scaling arrows
        max_force = 0.0

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
                max_w = max(abs(load.wx), abs(load.wy))
                max_force = max(max_force, max_w * el.length)
            elif isinstance(load, TriangularLoad):
                el = load.element
                max_w = max(abs(load.wx1), abs(load.wx2), abs(load.wy1), abs(load.wy2))
                max_force = max(max_force, max_w * el.length)

        if max_force == 0.0:
            max_force = 1.0

        # 1. Distributed Loads & Triangular Loads
        for load in self.structure.loads:
            if isinstance(load, DistributedLoad):
                el = load.element
                val = load.wy if abs(load.wy) >= abs(load.wx) else load.wx
                label = f"w = {abs(val):g}"
                self._draw_distributed_load(
                    ax, el, load.wx, load.wx, load.wy, load.wy, max_force, span, label_str=label
                )
            elif isinstance(load, TriangularLoad):
                el = load.element
                w1 = load.wy1 if abs(load.wy1) >= abs(load.wx1) else load.wx1
                w2 = load.wy2 if abs(load.wy2) >= abs(load.wx2) else load.wx2
                if w1 == 0:
                    label = f"w = {abs(w2):g}"
                elif w2 == 0:
                    label = f"w = {abs(w1):g}"
                else:
                    label = f"w₁={abs(w1):g}, w₂={abs(w2):g}"
                self._draw_distributed_load(
                    ax, el, load.wx1, load.wx2, load.wy1, load.wy2, max_force, span, label_str=label
                )

        # 2. Point Loads from Structure.loads (Global coordinates)
        for load in self.structure.loads:
            if isinstance(load, PointLoad):
                node = load.node
                if load.fx != 0 or load.fy != 0:
                    val = np.hypot(load.fx, load.fy)
                    self._draw_point_arrow(
                        ax, node.x, node.y, load.fx, load.fy, max_force, span, label_val=val, color="#C53030"
                    )
                if load.mz != 0:
                    self._draw_moment_arc(ax, node.x, node.y, load.mz, span, color="#7B1FA2")

        # 3. Element Point Loads (Local to Global transformation)
        for load in self.structure.loads:
            if isinstance(load, ElementPointLoad):
                el = load.element
                c, s = el.cos, el.sin
                t = load.x / el.length if el.length > 0 else 0.0
                px = el.node_i.x + t * (el.node_j.x - el.node_i.x)
                py = el.node_i.y + t * (el.node_j.y - el.node_i.y)

                # Local to Global
                fx_global = load.px * c - load.py * s
                fy_global = load.px * s + load.py * c

                if fx_global != 0 or fy_global != 0:
                    val = np.hypot(fx_global, fy_global)
                    self._draw_point_arrow(
                        ax, px, py, fx_global, fy_global, max_force, span, label_val=val, color="#D97706"
                    )
                if load.mz != 0:
                    self._draw_moment_arc(ax, px, py, load.mz, span, color="#7B1FA2")

        # 4. Node Loads from node.load attribute (Global coordinates)
        for node in self.structure.nodes.values():
            fx, fy, mz = node.load
            if fx != 0 or fy != 0:
                val = np.hypot(fx, fy)
                self._draw_point_arrow(
                    ax, node.x, node.y, fx, fy, max_force, span, label_val=val, color="#C53030"
                )
            if mz != 0:
                self._draw_moment_arc(ax, node.x, node.y, mz, span, color="#7B1FA2")

    def _draw_member_hinges(self, ax, span):
        """Draw open circular hinge symbols at member ends with moment releases."""
        hinge_radius = 0.012 * span
        offset_dist = 0.035 * span

        for el in self.structure.elements.values():
            hinge_i = getattr(el, "hinge_i", False)
            hinge_j = getattr(el, "hinge_j", False)
            if hinge_i or hinge_j:
                L = el.length
                if L < 1e-9:
                    continue
                dx = (el.node_j.x - el.node_i.x) / L
                dy = (el.node_j.y - el.node_i.y) / L

                if hinge_i:
                    hx = el.node_i.x + dx * min(offset_dist, 0.2 * L)
                    hy = el.node_i.y + dy * min(offset_dist, 0.2 * L)
                    h_circle = Circle(
                        (hx, hy),
                        radius=hinge_radius,
                        facecolor="white",
                        edgecolor="#2D3748",
                        linewidth=1.2,
                        zorder=6,
                    )
                    ax.add_patch(h_circle)

                if hinge_j:
                    hx = el.node_j.x - dx * min(offset_dist, 0.2 * L)
                    hy = el.node_j.y - dy * min(offset_dist, 0.2 * L)
                    h_circle = Circle(
                        (hx, hy),
                        radius=hinge_radius,
                        facecolor="white",
                        edgecolor="#2D3748",
                        linewidth=1.2,
                        zorder=6,
                    )
                    ax.add_patch(h_circle)

    def draw(
        self,
        show_undeformed=True,
        show_deformed=True,
        show_loads=True,
        show_supports=True,
        show_node_labels=True,
        show_element_labels=False,
        show_grid=False,
        title=None,
        n_points=40,
        support_style="detailed",
        figsize=(10, 7),
        save_path=None,
        show=True,
        ax=None,
    ):
        """
        Draw the structural model with professional styling, loads, boundary conditions,
        and deformed shape.

        Parameters
        ----------
        show_undeformed : bool, optional
            Whether to draw the undeformed structure. Defaults to True.
        show_deformed : bool, optional
            Whether to draw the deformed structure. Defaults to True.
        show_loads : bool, optional
            Whether to draw applied loads. Defaults to True.
        show_supports : bool, optional
            Whether to draw boundary supports. Defaults to True.
        show_node_labels : bool, optional
            Whether to show node ID labels. Defaults to True.
        show_element_labels : bool, optional
            Whether to show element ID labels at member midpoints. Defaults to False.
        show_grid : bool, optional
            Whether to show background grid lines. Defaults to False.
        title : str, optional
            Custom figure title. If None, a professional dynamic title is generated.
        n_points : int, optional
            Number of points for plotting deformed shape curves. Defaults to 40.
        support_style : str, optional
            'detailed' for realistic fixed, pin, and roller symbols, or 'box' for shaded box pads.
        figsize : tuple of float, optional
            Figure dimensions in inches. Defaults to (10, 7).
        save_path : str, optional
            Path to save figure (e.g. 'structure.png' or 'structure.pdf'). Defaults to None.
        show : bool, optional
            Whether to call plt.show(). Defaults to True.
        ax : matplotlib.axes.Axes, optional
            Existing axes to draw on. If None, a new figure and axes are created.

        Returns
        -------
        matplotlib.axes.Axes
            The plot axes.
        """
        if plt is None:
            raise ImportError("matplotlib is required to draw structures")

        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
        else:
            fig = ax.figure

        span = self._get_span()

        all_x = [node.x for node in self.structure.nodes.values()]
        all_y = [node.y for node in self.structure.nodes.values()]

        # 1. Undeformed Geometry (Crisp Slate)
        if show_undeformed:
            first = True
            for el in self.structure.elements.values():
                x_vals = [el.node_i.x, el.node_j.x]
                y_vals = [el.node_i.y, el.node_j.y]
                ax.plot(
                    x_vals,
                    y_vals,
                    color="#2C3E50",
                    linewidth=2.0,
                    label="Undeformed" if first else "",
                    zorder=3,
                )
                first = False

        # 2. Deformed Geometry (Vibrant Red)
        has_deformed = show_deformed and self.structure.disp is not None
        if has_deformed:
            first = True
            for el in self.structure.elements.values():
                points = el.deformed_shape_points(self.structure.disp, n_points, self.scale)
                x_def = [p[0] for p in points]
                y_def = [p[1] for p in points]
                all_x.extend(x_def)
                all_y.extend(y_def)
                label_txt = (
                    f"Deformed (x{self.scale:g})"
                    if self.scale != 1.0
                    else "Deformed"
                )
                ax.plot(
                    x_def,
                    y_def,
                    color="#E53E3E",
                    linewidth=2.0,
                    linestyle="-",
                    label=label_txt if first else "",
                    zorder=4,
                )
                first = False

        # 3. Element Labels (Optional)
        if show_element_labels:
            for el in self.structure.elements.values():
                mid_x = (el.node_i.x + el.node_j.x) / 2.0
                mid_y = (el.node_i.y + el.node_j.y) / 2.0
                ax.text(
                    mid_x,
                    mid_y,
                    f"E{el.id}",
                    fontsize=7.5,
                    fontweight="bold",
                    color="#4A5568",
                    ha="center",
                    va="center",
                    bbox=dict(boxstyle="circle,pad=0.2", fc="#EDF2F7", ec="#CBD5E0", lw=0.6, alpha=0.85),
                    zorder=5,
                )

        # 4. Member Hinges
        self._draw_member_hinges(ax, span)

        # 5. Nodes
        for node in self.structure.nodes.values():
            ax.plot(
                node.x,
                node.y,
                marker="o",
                markersize=5,
                color="#1A202C",
                markeredgecolor="white",
                markeredgewidth=1.2,
                zorder=5,
            )
            if show_node_labels:
                offset = 0.015 * span
                ax.text(
                    node.x + offset,
                    node.y + offset,
                    f"{node.id}",
                    fontsize=8,
                    fontweight="bold",
                    color="#1A202C",
                    ha="left",
                    va="bottom",
                    bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.75),
                    zorder=6,
                )

        # 6. Boundary Supports
        if show_supports:
            for node in self.structure.nodes.values():
                self._draw_support(ax, node, span, support_style=support_style)

        # 7. Applied Loads
        if show_loads:
            self._draw_loads(ax, span)

        # 8. Legend
        handles, labels = ax.get_legend_handles_labels()
        if handles:
            ax.legend(
                handles,
                labels,
                loc="upper right",
                frameon=True,
                framealpha=0.92,
                facecolor="white",
                edgecolor="#CBD5E0",
                fontsize=8.5,
            )

        # 9. Title & Labels
        if title is None:
            if has_deformed:
                title = f"Structural Model & Deformed Shape (Scale: {self.scale:g}x)"
            else:
                title = "Structural Model"

        if title:
            ax.set_title(title, fontsize=12, fontweight="bold", color="#1A202C", pad=14)

        ax.set_xlabel("X", fontsize=9, color="#4A5568", labelpad=6)
        ax.set_ylabel("Y", fontsize=9, color="#4A5568", labelpad=6)

        # 10. Styling, Spines & Grid
        if show_grid:
            ax.grid(True, linestyle="--", linewidth=0.6, alpha=0.4, color="#CBD5E0")
        else:
            ax.grid(False)

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_color("#CBD5E0")
        ax.spines["bottom"].set_color("#CBD5E0")
        ax.tick_params(colors="#718096", labelsize=8)

        # Auto-compute bounding margins to prevent clipping of load arrows and supports
        if all_x and all_y:
            min_x, max_x = min(all_x), max(all_x)
            min_y, max_y = min(all_y), max(all_y)
            dx = max_x - min_x
            dy = max_y - min_y
            margin = max(dx, dy) * 0.18
            margin = max(margin, 0.15 * span)
            ax.set_xlim(min_x - margin, max_x + margin)
            ax.set_ylim(min_y - margin, max_y + margin)

        ax.set_aspect("equal", adjustable="box")
        plt.tight_layout()

        if save_path:
            fig.savefig(save_path, dpi=300, bbox_inches="tight")

        if show:
            plt.show()

        return ax

    def draw_mode_shape(
        self,
        mode=1,
        phi=None,
        period=None,
        omega=None,
        scale=None,
        auto_scale=True,
        show_undeformed=True,
        title=None,
        **kwargs,
    ):
        """
        Draw a vibration mode shape of the structure.

        Parameters
        ----------
        mode : int, optional
            1-based mode number to draw (e.g. 1 for fundamental mode, 2 for second mode, etc.). Defaults to 1.
        phi : numpy.ndarray, optional
            Matrix of mode shapes (columns are modes) or 1D array for the single mode shape.
            Can be in reduced free DOFs or full structure DOFs.
            If None, modal analysis is retrieved from `structure.modal_results` or computed automatically.
        period : float, optional
            Natural period T (in seconds) to display in the plot title.
        omega : float, optional
            Natural frequency (rad/s) to display if period is not provided.
        scale : float, optional
            Displacement magnification factor. If None and auto_scale is True, the mode shape
            is automatically scaled to produce an aesthetically pleasing deflection (~10% of span).
        auto_scale : bool, optional
            Whether to automatically scale the mode shape if scale is not specified. Defaults to True.
        show_undeformed : bool, optional
            Whether to show the undeformed geometry in the background. Defaults to True.
        title : str, optional
            Custom figure title. If None, a title with the mode number and period is generated.
        **kwargs : dict
            Additional keyword arguments forwarded to `draw()` (e.g. `show_node_labels`, `save_path`, `show`, `ax`).

        Returns
        -------
        matplotlib.axes.Axes
            The plot axes.
        """
        if mode < 1:
            raise ValueError("Mode number must be >= 1.")
        mode_idx = mode - 1

        # Ensure mode shape data is available
        if phi is None:
            if (
                hasattr(self.structure, "modal_results")
                and self.structure.modal_results is not None
            ):
                phi = self.structure.modal_results["phi"]
                # Also pull period/omega from stored results if caller didn't supply them
                stored_omega = self.structure.modal_results["omega"]
                if period is None and omega is None and mode_idx < len(stored_omega):
                    omega = stored_omega[mode_idx]
                    import math
                    period = 2 * math.pi / omega
            else:
                raise ValueError(
                    "No mode shape data found. Either pass 'phi' directly or first call "
                    "structure.set_modal_results(omega, phi) to register modal results."
                )


        phi = np.asarray(phi)
        if phi.ndim == 2:
            if mode_idx >= phi.shape[1]:
                raise IndexError(f"Mode {mode} requested, but phi only contains {phi.shape[1]} modes.")
            mode_vec = phi[:, mode_idx]
        elif phi.ndim == 1:
            mode_vec = phi
        else:
            raise ValueError("phi must be a 1D vector or 2D matrix.")

        # Ensure boundary conditions are applied
        if not hasattr(self.structure, "neq") or self.structure.free_dofs is None:
            self.structure.number_dofs()
            self.structure.apply_boundary_conditions()

        # Construct full displacement vector
        if len(mode_vec) == len(self.structure.free_dofs):
            disp_full = np.zeros(self.structure.neq)
            disp_full[self.structure.free_dofs] = mode_vec
        elif len(mode_vec) == self.structure.neq:
            disp_full = np.array(mode_vec, dtype=float)
        else:
            raise ValueError(
                f"Length of mode vector ({len(mode_vec)}) does not match either free DOFs ({len(self.structure.free_dofs)}) or full DOFs ({self.structure.neq})."
            )

        # Scale determination
        span = self._get_span()
        if scale is not None:
            effective_scale = scale
        elif auto_scale:
            max_trans = 0.0
            for node in self.structure.nodes.values():
                ux = disp_full[node.dofs[0]]
                uy = disp_full[node.dofs[1]]
                max_trans = max(max_trans, np.hypot(ux, uy))
            if max_trans > 1e-12:
                effective_scale = (0.10 * span) / max_trans
            else:
                effective_scale = 1.0
        else:
            effective_scale = self.scale

        # Title formatting
        if title is None:
            if period is not None and np.isfinite(period):
                title = f"Mode {mode} Shape (T = {period:.4f} s)"
            elif omega is not None and np.isfinite(omega):
                title = f"Mode {mode} Shape (ω = {omega:.2f} rad/s)"
            else:
                title = f"Mode {mode} Shape"

        # Temporary displacement and scale swap
        old_disp = self.structure.disp
        old_scale = self.scale
        self.structure.disp = disp_full
        self.scale = effective_scale
        show_loads = kwargs.pop("show_loads", False)

        try:
            ax = self.draw(
                show_undeformed=show_undeformed,
                show_deformed=True,
                show_loads=show_loads,
                title=title,
                **kwargs,
            )
        finally:
            self.structure.disp = old_disp
            self.scale = old_scale

        return ax


def plot_mode_shape(structure_or_frame, mode=1, phi=None, period=None, omega=None, scale=None, auto_scale=True, **kwargs):
    """
    Convenience function to plot a structural vibration mode shape.

    Parameters
    ----------
    structure_or_frame : Structure or SimpleFrame
        The structural model.
    mode : int, optional
        1-based mode number to draw. Defaults to 1.
    phi : numpy.ndarray, optional
        Mode shapes matrix or vector.
    period : float, optional
        Natural period for title annotation.
    omega : float, optional
        Natural circular frequency for title annotation.
    scale : float, optional
        Displacement magnification factor.
    auto_scale : bool, optional
        Whether to automatically scale the mode shape. Defaults to True.
    **kwargs : dict
        Additional options passed to DrawStructure.draw_mode_shape.

    Returns
    -------
    matplotlib.axes.Axes
        The plot axes.
    """
    struct = (
        structure_or_frame.structure
        if hasattr(structure_or_frame, "structure")
        else structure_or_frame
    )
    drawer = DrawStructure(struct)
    return drawer.draw_mode_shape(
        mode=mode,
        phi=phi,
        period=period,
        omega=omega,
        scale=scale,
        auto_scale=auto_scale,
        **kwargs,
    )
