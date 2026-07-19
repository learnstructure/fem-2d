"""
fem2d: A simple 2D finite element analysis library in Python.

This library provides classes and utilities for performing 2D finite element analysis (FEA)
of structural frames, including truss, beam, and spring elements with linear and 
geometrically non-linear solver options.

Author: Abinash Mandal
"""

from .structure import Structure
from .nodes import Node
from .materials import ElasticMaterial, UniaxialMaterial, BilinearMaterial, CSVMaterial
from .sections import Section, FiberSection, MomentCurvatureSection
from .results import Results
from .elements.beam import BeamElement
from .elements.beam_hinges import BeamWithHingesElement
from .elements.truss import TrussElement
from .elements.trussNL import TrussElementNL
from .elements.beamNL import BeamElementNL
from .elements.beam_materialNL import BeamElementMaterialNL
from .utils.simple_frame import SimpleFrame
from .solver import NewtonRaphsonSolver
from .buckling_analysis import buckling_analysis

try:
    from .utils.draw_structure import DrawStructure
except ImportError:  # pragma: no cover - optional plotting dependency
    DrawStructure = None

__all__ = [
    "Structure",
    "Node",
    "ElasticMaterial",
    "UniaxialMaterial",
    "BilinearMaterial",
    "CSVMaterial",
    "Section",
    "FiberSection",
    "MomentCurvatureSection",
    "Results",
    "BeamElement",
    "BeamWithHingesElement",
    "TrussElement",
    "TrussElementNL",
    "BeamElementNL",
    "BeamElementMaterialNL",
    "SimpleFrame",
    "DrawStructure",
    "NewtonRaphsonSolver",
    "buckling_analysis",
]
