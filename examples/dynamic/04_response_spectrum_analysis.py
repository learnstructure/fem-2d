"""
Example: Response Spectrum Analysis (RSA) of a 4-Story 2D Building Frame.

Demonstrates:
1. Modal analysis of a multi-story building frame with lumped floor masses.
2. Response Spectrum Analysis using continuous ASCE 7 design spectrum (SDS, SD1).
3. Response Spectrum Analysis using El Centro ground motion record.
4. Calculation of modal participation factors, modal effective masses, and cumulative MMPR.
5. Comparison between SRSS and CQC modal combination rules.
6. Extraction of peak member forces (axial, shear, bending moment).
"""

import numpy as np
from fem2d.utils.simple_frame import SimpleFrame
from fem2d import Section, ElasticMaterial
from structdyn.ground_motions.ground_motion import GroundMotion
from structdyn.utils.helpers import elcentro_chopra
from struct_load import SeismicParameters

# -----------------------------------------------------------------------------
# 1. Geometry and Material Properties (Units: kips, inches, seconds)
# -----------------------------------------------------------------------------
g = 386.09  # in/s²
bay_width = 16.0 * 12.0  # 192 in (16 ft)
story_height = 10.0 * 12.0  # 120 in (10 ft)
n_stories = 4

E = 6000.0  # ksi
material = ElasticMaterial(E)

# Concrete Cross Sections
# Beam: 20 in x 35 in
beam_sec = Section.from_rectangle(width=20.0, depth=35.0)
# Column: 20 in x 20 in
col_sec = Section.from_rectangle(width=20.0, depth=20.0)

# -----------------------------------------------------------------------------
# 2. Build Structural Model
# -----------------------------------------------------------------------------
frame = SimpleFrame()
nodes = []

# Base nodes (fixed)
node_lb = frame.add_node(1, 0.0, 0.0)
node_rb = frame.add_node(2, bay_width, 0.0)
nodes.extend([node_lb, node_rb])

# Story nodes
for story in range(1, n_stories + 1):
    y = story * story_height
    nl = frame.add_node(2 * story + 1, 0.0, y)
    nr = frame.add_node(2 * story + 2, bay_width, y)
    nodes.extend([nl, nr])

# Beams
for story in range(1, n_stories + 1):
    nl = nodes[2 * story]
    nr = nodes[2 * story + 1]
    frame.add_frame(story, nl.id, nr.id, E, beam_sec.A, beam_sec.I)

# Columns (Left and Right)
for story in range(n_stories):
    # Left column
    bl = nodes[2 * story]
    tl = nodes[2 * (story + 1)]
    frame.add_frame(n_stories + story + 1, bl.id, tl.id, E, col_sec.A, col_sec.I)

    # Right column
    br = nodes[2 * story + 1]
    tr = nodes[2 * (story + 1) + 1]
    frame.add_frame(2 * n_stories + story + 1, br.id, tr.id, E, col_sec.A, col_sec.I)

# Fixed base supports
frame.add_support(1, [1, 1, 1])
frame.add_support(2, [1, 1, 1])

# Lumped story masses: 450 kips per floor (225 kips per node)
weight_per_node = 450.0 / 2.0  # kips
mass_per_node = weight_per_node / g  # kips*s²/in (slug-equivalent)
small_inertia = 1e-6

for story in range(1, n_stories + 1):
    nl_id = 2 * story + 1
    nr_id = 2 * story + 2
    frame.structure.nodes[nl_id].set_mass(mass=mass_per_node, inertia=small_inertia)
    frame.structure.nodes[nr_id].set_mass(mass=mass_per_node, inertia=small_inertia)

print("=" * 70)
print("FEM2D: RESPONSE SPECTRUM ANALYSIS OF 4-STORY CONCRETE FRAME")
print("=" * 70)

# -----------------------------------------------------------------------------
# 3. RSA Method 1: Continuous ASCE 7 Design Spectrum (SDS=1.0g, SD1=0.6g)
# -----------------------------------------------------------------------------
print("\n--- [Analysis 1] ASCE 7 Design Spectrum (SDS = 1.0g, SD1 = 0.6g) ---")
asce_params = SeismicParameters(sds=1.0, sd1=0.6, tl=8.0, site_class="D", damping_ratio=0.05)
rsa_asce = frame.response_spectrum_analysis(
    params=asce_params,
    direction="x",
    scale_factor=g,  # Convert g to in/s²
    combination="CQC",
)

periods = rsa_asce["periods"]
gammas = rsa_asce["participation_factors"]
mmpr = rsa_asce["mass_participation_ratio"]

print("Modal Periods (s):")
for i, t in enumerate(periods[:4]):
    print(f"  Mode {i+1}: T = {t:.4f} s | Gamma = {gammas[i]:.4f} | Mass Ratio = {mmpr[i]*100:.2f}%")

print(f"Cumulative Effective Mass Participation: {np.sum(mmpr[:4])*100:.2f}%")
print(f"Combined Base Shear (CQC): {float(rsa_asce['combined_base_shear'][0]):.2f} kips")

# Compare CQC vs SRSS
srss_vb = rsa_asce["rsa"].combine_modal_responses(rsa_asce["modal_base_shear"], method="SRSS")
print(f"Combined Base Shear (SRSS): {float(srss_vb[0]):.2f} kips")

# Member Forces for Base Column (Element 5)
col_forces = rsa_asce["combined_element_forces"][5]
print("\nBase Column (Element 5) Combined Internal Forces:")
print(f"  Axial P:       {abs(col_forces[0]):.2f} kips")
print(f"  Shear V:       {abs(col_forces[1]):.2f} kips")
print(f"  Base Moment M: {abs(col_forces[2]):.2f} kip-in ({abs(col_forces[2])/12:.2f} kip-ft)")

# -----------------------------------------------------------------------------
# 4. RSA Method 2: Ground Motion Spectrum (El Centro Record)
# -----------------------------------------------------------------------------
print("\n--- [Analysis 2] RSA from El Centro Ground Motion ---")
elc = elcentro_chopra()
gm = GroundMotion.from_arrays(elc["acc (g)"], 0.02)

rsa_gm = frame.response_spectrum_analysis(
    gm=gm,
    ji=0.05,
    direction="x",
    scale_factor=g,  # Convert g to in/s²
    combination="CQC",
)

print(f"El Centro Base Shear (CQC):  {float(rsa_gm['combined_base_shear'][0]):.2f} kips")
srss_gm = rsa_gm["rsa"].combine_modal_responses(rsa_gm["modal_base_shear"], method="SRSS")
print(f"El Centro Base Shear (SRSS): {float(srss_gm[0]):.2f} kips")

# Roof lateral displacement (Node 9 or 10, X-translation)
node_roof = frame.structure.nodes[9]
roof_ux = frame.structure.disp[node_roof.dofs[0]]
print(f"Peak Roof Lateral Displacement: {roof_ux:.4f} in")
print("=" * 70)
