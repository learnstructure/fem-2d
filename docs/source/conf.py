import os
import sys
sys.path.insert(0, os.path.abspath("../../"))

project = "fem2d"
copyright = "2026, Abinash Mandal"
author = "Abinash Mandal"
release = "0.2.0"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
    "sphinx.ext.todo",
    "myst_parser"   
]

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]