"""
Custom material model defined by an envelope curve from a CSV file.
"""

import os
import numpy as np
import pandas as pd
from fem2d.materials.material import UniaxialMaterial


class CSVMaterial(UniaxialMaterial):
    """
    Uniaxial material defined by a custom curve loaded from a CSV file.

    Supports:
    - Path-independent (nonlinear elastic) behavior.
    - Path-dependent (elastoplastic with elastic unloading/reloading) behavior.
    - Automatic symmetry mirroring for positive-only data curves.
    """

    def __init__(
        self,
        csv_path,
        x_col=0,
        y_col=1,
        path_dependent=False,
        extrapolate_slope="flat",
        rho=0.0,
    ):
        """
        Initialize the CSVMaterial.

        Parameters
        ----------
        csv_path : str
            Absolute or relative path to the CSV file.
        x_col : int or str, optional
            Column index or name for deformation/strain values. Defaults to 0.
        y_col : int or str, optional
            Column index or name for force/stress values. Defaults to 1.
        path_dependent : bool, optional
            If True, uses elastoplastic unloading/reloading. If False, follows
            the envelope curve for both loading and unloading. Defaults to False.
        extrapolate_slope : str, optional
            Extrapolation style beyond the bounds of the CSV points:
            - 'flat': Tangent slope becomes 0; force/stress remains at the boundary value.
            - 'linear': Extrapolates using the slope of the final segment.
            Defaults to 'flat'.
        rho : float, optional
            Mass density of the material. Defaults to 0.0.
        """
        self.csv_path = csv_path
        self.path_dependent = path_dependent
        self.extrapolate_slope = extrapolate_slope
        self.rho = rho

        # Load CSV
        if not os.path.exists(csv_path):
            raise FileNotFoundError(f"CSV material file not found: {csv_path}")

        df = pd.read_csv(csv_path)

        # Extract columns
        if isinstance(x_col, int):
            self.x_data = df.iloc[:, x_col].to_numpy(dtype=float)
        else:
            self.x_data = df[x_col].to_numpy(dtype=float)

        if isinstance(y_col, int):
            self.y_data = df.iloc[:, y_col].to_numpy(dtype=float)
        else:
            self.y_data = df[y_col].to_numpy(dtype=float)

        # Ensure sorted by x
        sort_idx = np.argsort(self.x_data)
        self.x_data = self.x_data[sort_idx]
        self.y_data = self.y_data[sort_idx]

        # Check if the curve needs (0,0) prepended/inserted
        if not np.any(np.isclose(self.x_data, 0.0)):
            # Insert (0,0)
            idx = np.searchsorted(self.x_data, 0.0)
            self.x_data = np.insert(self.x_data, idx, 0.0)
            self.y_data = np.insert(self.y_data, idx, 0.0)

        # Check for symmetry (if all values are non-negative, assume symmetric)
        self.is_symmetric = np.all(self.x_data >= 0.0)

        if self.is_symmetric:
            # Build full symmetric curve for internal interpolation
            # Negative coordinates: reflect the positive coordinates (excluding zero)
            pos_mask = self.x_data > 0.0
            x_neg = -self.x_data[pos_mask][::-1]
            y_neg = -self.y_data[pos_mask][::-1]

            self.x_full = np.concatenate([x_neg, self.x_data])
            self.y_full = np.concatenate([y_neg, self.y_data])
        else:
            self.x_full = self.x_data
            self.y_full = self.y_data

        # Determine initial elastic modulus E (from the first positive segment)
        zero_idx = np.where(np.isclose(self.x_full, 0.0))[0][0]
        if zero_idx < len(self.x_full) - 1:
            dx = self.x_full[zero_idx + 1] - self.x_full[zero_idx]
            dy = self.y_full[zero_idx + 1] - self.y_full[zero_idx]
            self.E_init = dy / dx if dx != 0.0 else 0.0
        else:
            self.E_init = 1.0  # fallback

        # Initial state for path-dependent unloading/reloading
        self.init_state = {
            "eps_max": 0.0,
            "eps_min": 0.0,
            "eps_c": 0.0,
            "sigma_c": 0.0,
        }

    def _interpolate_envelope(self, x):
        """Interpolate y and tangent slope dy/dx from full envelope curve."""
        x_pts = self.x_full
        y_pts = self.y_full

        # Handle boundaries / extrapolation
        if x < x_pts[0]:
            if self.extrapolate_slope == "flat":
                return y_pts[0], 0.0
            else:
                slope = (y_pts[1] - y_pts[0]) / (x_pts[1] - x_pts[0])
                return y_pts[0] + slope * (x - x_pts[0]), slope

        elif x > x_pts[-1]:
            if self.extrapolate_slope == "flat":
                return y_pts[-1], 0.0
            else:
                slope = (y_pts[-1] - y_pts[-2]) / (x_pts[-1] - x_pts[-2])
                return y_pts[-1] + slope * (x - x_pts[-1]), slope

        # Inside bounds: find segment
        idx = np.searchsorted(x_pts, x)
        if idx == 0:
            return y_pts[0], 0.0

        x0, x1 = x_pts[idx - 1], x_pts[idx]
        y0, y1 = y_pts[idx - 1], y_pts[idx]

        dx = x1 - x0
        slope = (y1 - y0) / dx if dx != 0.0 else 0.0
        y = y0 + slope * (x - x0)
        return y, slope

    def get_response(self, deformation, state=None):
        """Compute current force/stress and tangent slope."""
        if not self.path_dependent:
            # Path-independent: directly query envelope
            y, slope = self._interpolate_envelope(deformation)
            return y, slope, state

        # Path-dependent logic
        if state is None:
            state = self.init_state.copy()

        eps_max_c = state["eps_max"]
        eps_min_c = state["eps_min"]
        eps_c = state["eps_c"]
        sigma_c = state["sigma_c"]

        # 1. Check if we exceed historic envelope bounds
        if deformation > eps_max_c:
            # Loading on positive envelope
            y, slope = self._interpolate_envelope(deformation)
            new_state = {
                "eps_max": deformation,
                "eps_min": eps_min_c,
                "eps_c": deformation,
                "sigma_c": y,
            }
            return y, slope, new_state

        elif deformation < eps_min_c:
            # Loading on negative envelope
            y, slope = self._interpolate_envelope(deformation)
            new_state = {
                "eps_max": eps_max_c,
                "eps_min": deformation,
                "eps_c": deformation,
                "sigma_c": y,
            }
            return y, slope, new_state

        else:
            # Unloading/reloading elastically with slope E_init
            trial_sigma = sigma_c + self.E_init * (deformation - eps_c)

            # Cap trial stress to envelope values at current deformation
            y_max, _ = self._interpolate_envelope(eps_max_c)
            y_min, _ = self._interpolate_envelope(eps_min_c)

            # Capping: cannot exceed envelope stress at the boundaries
            if trial_sigma > y_max:
                trial_sigma = y_max
                slope = 0.0
            elif trial_sigma < y_min:
                trial_sigma = y_min
                slope = 0.0
            else:
                slope = self.E_init

            new_state = {
                "eps_max": eps_max_c,
                "eps_min": eps_min_c,
                "eps_c": deformation,
                "sigma_c": trial_sigma,
            }
            return trial_sigma, slope, new_state

    def get_initial_tangent(self):
        """Return the initial elastic slope E_init."""
        return self.E_init
