"""
Adapters bridging fem2d with external ecosystems.

Currently exposes:

- :func:`fem2d.adapters.struct_core_adapter.model_from_core` — convert a
  ``struct_core`` model into a ``fem2d.Structure``.
- :func:`fem2d.adapters.struct_core_adapter.model_to_core` — convert a
  ``fem2d.Structure`` into a ``struct_core.StructuralModel``.
- :func:`fem2d.adapters.struct_core_adapter.result_to_core` — convert fem2d
  analysis output into a ``struct_core.AnalysisResult``.

Re-exported at the top level as ``fem2d.model_from_core``,
``fem2d.model_to_core``, and ``fem2d.result_to_core`` for convenience.
"""

try:
    from .struct_core_adapter import model_from_core, model_to_core, result_to_core, solve_combinations_and_envelope
except ImportError:  # struct_core is an optional dependency
    def _missing(*args, **kwargs):
        raise ImportError(
            "The 'struct_core' package is required for this feature.\n"
            "Install it with:  pip install struct_core\n"
            "or install the optional extra:  pip install fem2d[ecosystem]"
        )
    model_from_core = model_to_core = result_to_core = solve_combinations_and_envelope = _missing

__all__ = ["model_from_core", "model_to_core", "result_to_core", "solve_combinations_and_envelope"]