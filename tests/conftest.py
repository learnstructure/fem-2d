"""
Pytest configuration for the fem2d test suite.

The ``structeng`` conda environment ships with a placeholder ``struct_core``
namespace package (v0.1.0). To exercise the real development version of
``struct_core`` that lives alongside ``fem-2d`` at the project root, we
prepend its source directory to ``sys.path`` before any test module imports
``struct_core``. This is opt-out via the environment variable
``FEM2D_KEEP_INSTALLED_STRUCT_CORE`` (e.g. for CI runners that install
``struct_core`` from PyPI).
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

_KEEP_ENV = "FEM2D_KEEP_INSTALLED_STRUCT_CORE"


def _maybe_prepend_dev_struct_core() -> None:
    if os.environ.get(_KEEP_ENV):
        return

    here = Path(__file__).resolve().parent
    repo_root = here.parent.parent  # tests/ -> fem-2d/ -> structural_engineering/
    dev_struct_core = repo_root / "struct_core"

    if not (dev_struct_core / "struct_core" / "__init__.py").is_file():
        # Source tree not present (e.g. isolated CI); leave sys.path untouched.
        return

    dev_path = str(dev_struct_core)
    if dev_path not in sys.path:
        sys.path.insert(0, dev_path)


_maybe_prepend_dev_struct_core()