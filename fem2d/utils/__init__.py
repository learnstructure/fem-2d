"""
Utilities package for model generation, frame creation, and visualization.
"""

from .simple_frame import SimpleFrame

try:
    from .draw_structure import DrawStructure
except ImportError:  # pragma: no cover - optional plotting dependency
    DrawStructure = None

__all__ = ["SimpleFrame", "DrawStructure"]
