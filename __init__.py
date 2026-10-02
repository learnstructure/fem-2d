"""
Dynamic forwarding to inner fem2d package.
This file automatically executes the inner fem2d/__init__.py.
You never need to edit or duplicate exports in this file.
"""

import os

_inner_dir = os.path.join(os.path.dirname(__file__), "fem2d")
if _inner_dir not in __path__:
    __path__.insert(0, _inner_dir)

_inner_init = os.path.join(_inner_dir, "__init__.py")
if os.path.exists(_inner_init):
    with open(_inner_init, "r", encoding="utf-8") as _f:
        exec(_f.read(), globals())
