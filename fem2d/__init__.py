"""
fem2d: A simple 2D finite element analysis library in Python.

This library provides classes and utilities for performing 2D finite element analysis (FEA)
of structural frames, including truss, beam, and spring elements with linear and 
geometrically non-linear solver options.

Author: Abinash Mandal
"""

from .structure import Structure
from .nodes import Node
from .materials import ElasticMaterial
from .sections import Section
from .results import Results
from .elements.beam import BeamElement
from .elements.beam_hinges import BeamWithHingesElement
from .elements.truss import TrussElement
from .elements.trussNL import TrussElementNL
from .elements.beamNL import BeamElementNL
from .utils.simple_frame import SimpleFrame
from .solver import NewtonRaphsonSolver

try:
    from .utils.draw_structure import DrawStructure
except ImportError:  # pragma: no cover - optional plotting dependency
    DrawStructure = None

__all__ = [
    "Structure",
    "Node",
    "ElasticMaterial",
    "Section",
    "Results",
    "BeamElement",
    "BeamWithHingesElement",
    "TrussElement",
    "TrussElementNL",
    "BeamElementNL",
    "SimpleFrame",
    "DrawStructure",
    "NewtonRaphsonSolver",
]
