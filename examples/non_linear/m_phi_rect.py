"""
Moment-curvature analysis of a rectangular section using fiber discretization.

Units: SI (N, m). Stress in Pa, lengths in m, moments in N-m.
"""

import numpy as np
import matplotlib.pyplot as plt

from fem2d.materials import BilinearMaterial
from fem2d.sections import FiberSection, run_moment_curvature


# Section geometry
b = 0.10   # width  (10 cm)
h = 0.20   # depth  (20 cm)

# Material: bilinear elasto-plastic
E = 200.0e9
sigma_y = 250.0e6
E_t = 0.0e9
material = BilinearMaterial(E, sigma_y, E_t)

# Fiber section
section = FiberSection.from_rectangular_fibers(
    width=b, depth=h, material=material, n_fibers=40
)

# Closed-form reference values for a rectangle (elastic-perfectly plastic)
# I = b * h ** 3 / 12.0
# M_y = sigma_y * b * h ** 2 / 6.0
# M_p = sigma_y * b * h ** 2 / 4.0
# phi_y = M_y / (E * I)

# print(f"I     = {I:.6e} m^4")
# print(f"M_y   = {M_y:.3f} N-m ({M_y/1e3:.3f} kN-m)")
# print(f"M_p   = {M_p:.3f} N-m ({M_p/1e3:.3f} kN-m)")
# print(f"phi_y = {phi_y:.6e} 1/m")

# Moment-curvature analysis
kappas, moments, axial, eps0, _ = run_moment_curvature(section, n_steps=400, N_target=-2500e3)
# print(eps0)
print(kappas[-1], moments[-1], axial[-1])
# print(axial)
# Plot
fig, ax = plt.subplots(figsize=(8, 5.5))
ax.plot(kappas, moments, color="steelblue", lw=2.0, label="M-phi")
# ax.axhline(M_y, color="gray", ls="--", lw=1.0, label=r"$M_y$ (first yield)")
# ax.axhline(M_p, color="firebrick", ls="--", lw=1.0, label=r"$M_p$ (fully plastic)")
# ax.axvline(phi_y, color="gray", ls=":", lw=1.0, label=r"$\phi_y$")
ax.set_xlabel(r"Curvature  $\phi$  [1/m]")
ax.set_ylabel(r"Bending moment  $M$  [N-m]")
ax.set_title("Moment-Curvature of a Rectangular Section")
ax.grid(True, ls=":", alpha=0.6)
ax.legend(loc="lower right")
fig.tight_layout()
fig.savefig("rectangular_moment_curvature.png", dpi=150)
print("Saved figure to: rectangular_moment_curvature.png")
