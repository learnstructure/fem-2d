"""
Deprecated. Redirects imports to the new fem2d/materials package.
"""

from fem2d.materials.material import UniaxialMaterial
from fem2d.materials.elastic import ElasticMaterial
from fem2d.materials.bilinear import BilinearMaterial
from fem2d.materials.csv_material import CSVMaterial

__all__ = [
    "UniaxialMaterial",
    "ElasticMaterial",
    "BilinearMaterial",
    "CSVMaterial",
]
