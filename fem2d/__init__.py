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
    from .adapters import model_from_core, model_to_core, result_to_core, solve_combinations_and_envelope
except ImportError:  # struct_core is an optional dependency
    def _missing_adapter(*args, **kwargs):
        raise ImportError(
            "The 'struct_core' package is required for this feature.\n"
            "Install it with:  pip install struct_core\n"
            "or install the optional extra:  pip install fem2d[ecosystem]"
        )
    model_from_core = model_to_core = result_to_core = solve_combinations_and_envelope = _missing_adapter

try:
    from .utils.draw_structure import DrawStructure, plot_mode_shape
except ImportError:  # pragma: no cover - optional plotting dependency
    DrawStructure = None
    plot_mode_shape = None

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
    "plot_mode_shape",
    "NewtonRaphsonSolver",
    "buckling_analysis",
    "model_from_core",
    "model_to_core",
    "result_to_core",
]
