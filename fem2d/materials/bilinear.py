"""
Bilinear elasto-plastic material model with kinematic and/or isotropic hardening.
"""

import numpy as np
from fem2d.materials.material import UniaxialMaterial


class BilinearMaterial(UniaxialMaterial):
    """
    Uniaxial elasto-plastic material with bilinear hardening.
    Supports kinematic, isotropic, or mixed hardening.

    Attributes
    ----------
    E : float
        Initial elastic modulus.
    sigma_y : float
        Yield stress/strength.
    E_t : float
        Tangent modulus in the plastic range.
    beta : float, optional
        Hardening mixture parameter:
        - 0.0: Pure kinematic hardening (default)
        - 1.0: Pure isotropic hardening
        - intermediate values: Mixed hardening
    rho : float, optional
        Mass density of the material. Defaults to 0.0.
    """

    def __init__(self, E, sigma_y, E_t, beta=0.0, rho=0.0):
        """Initialize the BilinearMaterial."""
        self.E = E
        self.sigma_y = sigma_y
        self.E_t = E_t
        self.beta = beta
        self.rho = rho

        # Compute hardening modulus H
        # H is defined such that plastic tangent is E_t.
        # 1/E_t = 1/E + 1/H => H = E*E_t / (E - E_t)
        if self.E == self.E_t:
            self.H = 0.0
        else:
            self.H = (self.E * self.E_t) / (self.E - self.E_t)

        self.H_kin = self.H * (1.0 - self.beta)
        self.H_iso = self.H * self.beta

        # Initial state dictionary
        self.init_state = {
            "ep": 0.0,      # plastic strain
            "alpha": 0.0,   # back stress (kinematic hardening)
            "s_hard": 0.0,  # isotropic hardening stress increment
        }

    def get_response(self, deformation, state=None):
        """
        Compute stress and tangent modulus at a given strain using radial return mapping.

        Parameters
        ----------
        deformation : float
            Current total strain.
        state : dict, optional
            Committed history variables. If None, uses initial state.

        Returns
        -------
        stress : float
            Current stress.
        tangent : float
            Current tangent modulus.
        new_state : dict
            Updated/trial history variables.
        """
        if state is None:
            state = self.init_state.copy()

        ep_c = state["ep"]
        alpha_c = state["alpha"]
        s_hard_c = state["s_hard"]

        # 1. Trial elastic stress
        sigma_tr = self.E * (deformation - ep_c)

        # 2. Yield function
        f_tr = abs(sigma_tr - alpha_c) - (self.sigma_y + s_hard_c)

        if f_tr <= 1e-12:
            # Elastic response
            stress = sigma_tr
            tangent = self.E
            new_state = {
                "ep": ep_c,
                "alpha": alpha_c,
                "s_hard": s_hard_c,
            }
        else:
            # Plastic response (radial return mapping)
            delta_gamma = f_tr / (self.E + self.H_iso + self.H_kin)
            sign = np.sign(sigma_tr - alpha_c)

            stress = sigma_tr - sign * self.E * delta_gamma
            tangent = self.E_t

            new_state = {
                "ep": ep_c + delta_gamma * sign,
                "alpha": alpha_c + self.H_kin * delta_gamma * sign,
                "s_hard": s_hard_c + self.H_iso * delta_gamma,
            }

        return stress, tangent, new_state

    def get_initial_tangent(self):
        """Return initial elastic tangent E."""
        return self.E
