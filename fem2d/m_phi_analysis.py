"""
Moment-curvature driver for fiber sections under constant axial force.

Given a :class:`FiberSection`, this utility imposes a monotonic curvature
history and recovers the axial strain `eps_0` (i.e. the position of the
neutral axis) at each step so that the integrated axial force matches a
target `N_target` (defaults to 0.0, i.e. pure bending).
"""

import numpy as np


def run_moment_curvature(section, n_steps=400, kappa_max=None,
                         N_target=0.0, tol=1e-10, max_iter=50,
                         verbose=False):
    """
    Drive a fiber section through a curvature history under a constant
    axial force `N_target`.

    At each imposed curvature `kappa`, the axial strain `eps_0` is solved
    so that the resulting axial force `N(eps_0, kappa) = N_target`. The
    neutral axis is therefore free to shift as the section yields.

    Parameters
    ----------
    section : FiberSection
        Fiber-based cross-section. Each fiber must have a material with a
        ``get_response(strain, state)`` method (e.g. ``BilinearMaterial``).
    n_steps : int, optional
        Number of curvature increments between 0 and ``kappa_max``.
    kappa_max : float, optional
        Target maximum curvature. If ``None``, a multiple of the
        estimated yield curvature is used (requires the fiber material
        to expose ``E`` and ``sigma_y``).
    N_target : float, optional
        Constant axial force to enforce at each step. Defaults to 0.0
        (pure bending). Positive = tension, negative = compression.
    tol : float, optional
        Convergence tolerance on the residual axial force.
    max_iter : int, optional
        Maximum Newton iterations per step to recover ``eps_0``.
    verbose : bool, optional
        Print per-step information.

    Returns
    -------
    kappas : numpy.ndarray
        Imposed curvature history, shape ``(n_steps + 1,)``.
    moments : numpy.ndarray
        Resulting bending moment at each curvature, shape ``(n_steps + 1,)``.
    axial : numpy.ndarray
        Resulting axial force at each curvature, shape ``(n_steps + 1,)``.
        Should match ``N_target`` up to ``tol``.
    eps0 : numpy.ndarray
        Recovered axial strain at each curvature, shape ``(n_steps + 1,)``.
    state : list
        Final committed fiber states.
    """
    # Estimate a default kappa_max if not provided.
    if kappa_max is None:
        fiber_material = section.fibers[0]["material"]
        E0 = fiber_material.get_initial_tangent()
        sigma_y = getattr(fiber_material, "sigma_y", None)
        if sigma_y is None:
            kappa_max = 1.0
        else:
            # First yield curvature of the section ~ 2 * sigma_y / (E * h)
            depth = 2.0 * max(abs(f["y"]) for f in section.fibers)
            phi_y = 2.0 * sigma_y / (E0 * depth)
            kappa_max = 25.0 * phi_y

    kappas = np.linspace(0.0, kappa_max, n_steps + 1)
    moments = np.zeros(n_steps + 1)
    axial = np.zeros(n_steps + 1)
    eps0 = np.zeros(n_steps + 1)

    state = [None] * len(section.fibers)
    eps_0 = 0.0

    for i, kappa in enumerate(kappas):
        if verbose:
            print(f"step {i + 1}/{n_steps}: kappa = {kappa:.6e}")

        # Newton solve for eps_0 such that N(eps_0, kappa) = N_target
        eps_0_iter = eps_0
        for _ in range(max_iter):
            forces, tangent, new_state = section.get_response(
                np.array([eps_0_iter, kappa]), state
            )
            N = forces[0] - N_target
            if abs(N) < tol:
                break
            eps_0_iter -= N / tangent[0, 0]

        forces, _tangent, new_state = section.get_response(
            np.array([eps_0_iter, kappa]), state
        )
        state = new_state
        eps_0 = eps_0_iter
        eps0[i] = eps_0
        moments[i] = forces[1]
        axial[i] = forces[0]

    return kappas, moments, axial, eps0, state

