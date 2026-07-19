"""
Deprecated. Redirects imports to the new fem2d/sections package.
"""

from fem2d.sections.section import Section
from fem2d.sections.fiber import FiberSection
from fem2d.sections.moment_curvature import MomentCurvatureSection

__all__ = [
    "Section",
    "FiberSection",
    "MomentCurvatureSection",
]
