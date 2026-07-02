"""Buckling analysis utilities for FEM2D.

This module supports linear buckling analysis for truss and beam structures
using the standard eigenproblem:

    K_ff * phi = lambda * Kg_ff * phi

where Kg is the geometric stiffness matrix assembled with compression positive.
"""

import numpy as np
from scipy.linalg import eig
from .elements.truss import TrussElement
from .elements.beam import BeamElement


def _assemble_geometric_stiffness(structure):
    """Assemble the global geometric stiffness matrix Kg.

    The matrix is assembled using the element geometric stiffness matrices
    evaluated at the current axial forces.  Axial forces are taken as
    positive in compression (as returned by el.axial_force()).

    Parameters
    ----------
    structure : Structure
        A solved or unsolved FEM2D structure.

    Returns
    -------
    ndarray
        Global geometric stiffness matrix (neq x neq).
    """
    Kg = np.zeros((structure.neq, structure.neq))
    for el in structure.elements.values():
        if not isinstance(el, (TrussElement, BeamElement)):
            continue

        # axial force: compression positive (according to our convention)
        N = el.axial_force()
        Kg_el = el.geometric_stiffness(N)   # now returns positive definite for compression

        dofs = el.node_i.dofs + el.node_j.dofs
        for i, ii in enumerate(dofs):
            for j, jj in enumerate(dofs):
                Kg[ii, jj] += Kg_el[i, j]
    return Kg


def buckling_analysis(structure, num_modes=1):
    """Perform linear buckling analysis for truss/beam structures.

    The analysis solves the generalized eigenvalue problem:

        K_ff * phi = lambda * Kg_ff * phi

    where:
        K_ff  : elastic stiffness matrix (free dofs)
        Kg_ff : geometric stiffness matrix (free dofs) assembled with
                compression-positive axial forces.
        lambda : buckling load multiplier (positive for compression buckling)

    Parameters
    ----------
    structure : Structure
        A solved or unsolved FEM2D structure containing truss/beam elements.
        The structure must have boundary conditions applied and a valid load
        case (the load that defines the stress state for Kg).
    num_modes : int, optional
        Number of buckling modes to return. Defaults to 1.

    Returns
    -------
    tuple
        (buckling_factors, buckling_modes)
        - buckling_factors : 1D array of positive eigenvalues (sorted ascending)
        - buckling_modes   : 2D array (neq x num_modes) where each column is
                             a mode shape in full DOF space.
    """
    # Ensure the structure is solved and boundary conditions are applied
    if structure.disp is None:
        structure.solve()
    elif getattr(structure, "free_dofs", None) is None:
        structure.apply_boundary_conditions()

    if structure.K is None:
        structure.assemble_stiffness()

    if structure.free_dofs is None:
        structure.apply_boundary_conditions()

    # Extract free-dof submatrices
    free = structure.free_dofs
    K_ff = structure.K[np.ix_(free, free)]

    # Assemble Kg and extract free-dof part
    Kg = _assemble_geometric_stiffness(structure)
    Kg_ff = Kg[np.ix_(free, free)]

    if K_ff.size == 0 or Kg_ff.size == 0:
        raise ValueError("Insufficient DOFs for buckling analysis.")

    # Solve the generalized eigenvalue problem:
    #   K_ff * phi = mu * Kg_ff * phi
    # With our sign convention, mu = lambda (positive for compression buckling)
    eigenvals, eigenvecs = eig(K_ff, Kg_ff)

    # Keep only finite real eigenvalues (scipy returns complex)
    real_mask = np.isfinite(eigenvals) & (np.abs(np.imag(eigenvals)) < 1e-8)
    eigenvals = np.real(eigenvals[real_mask])
    eigenvecs = np.real(eigenvecs[:, real_mask])


    if eigenvals.size == 0:
        raise ValueError(
            "No eigenvalues found. Check boundary conditions and load case."
        )

    # Sort by absolute ascending buckling factor (critical load first)
    order = np.argsort(np.abs(eigenvals))
    eigenvals = eigenvals[order]
    eigenvecs = eigenvecs[:, order]

    # Keep only the requested number of modes
    n_modes = min(num_modes, eigenvals.size)
    buckling_factors = eigenvals[:n_modes]

    # Expand mode shapes from free DOFs to full DOF space
    buckling_modes = np.zeros((structure.neq, n_modes)) 
    for i in range(n_modes):
        mode = np.zeros(structure.neq)
        mode[free] = eigenvecs[:, i]
        buckling_modes[:, i] = mode

    return buckling_factors, buckling_modes