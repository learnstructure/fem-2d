# FEM2D

[![PyPI version](https://img.shields.io/pypi/v/fem2d)](https://pypi.org/project/fem2d)
[![Downloads](https://static.pepy.tech/badge/fem2d)](https://pepy.tech/project/fem2d)
[![GitHub repo](https://img.shields.io/badge/GitHub-fem--2d-blue?logo=github)](https://github.com/learnstructure/fem-2d)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

An open-source Python library for structural finite element analysis of 2D structures.

**FEM2D** is a Python package for performing 2D finite element analysis (FEA) of structural frames, including truss, beam, and spring elements with support for linear static analysis, geometrically non-linear analysis, and global mass matrix assembly for dynamic analysis.

Documentation: [![Documentation Status](https://readthedocs.org/projects/fem2d/badge/?version=latest)](https://fem2d.readthedocs.io/en/latest/)


## Features

- **Element Library:**
  - **Beam Element:** Elastic 2D beam elements supporting both classical Euler-Bernoulli and Timoshenko beam formulations with shear deformation ($\beta_s = \frac{12EI}{GA_s L^2}$). Supports uniform and varying member loads, moment releases (hinges at either end or general 3-DOF releases), and rotational/translational mass.
  - **BeamWithHingesElement:** Beam with optional moment releases at either end (retained for backward compatibility).
  - **Truss Element:** Pin-jointed bar elements with axial stiffness only. The same `TrussElement` class is used for linear and corotational geometrically-nonlinear analyses.
  - **Spring Element:** 2D elastic spring elements with customizable axial stiffness.
- **Analysis Types:**
  - **Linear Static Analysis:** Standard matrix analysis under nodal loads, distributed loads, and concentrated member loads.
  - **Geometrically Non-Linear Analysis:** Iterative solver using the Newton-Raphson scheme combined with corotational formulations for large displacement/rotation problems. Enabled by passing `geometric_nonlinear=True` to `Structure.solve(...)`.
  - **Mass Matrix Assembly:** Assembles global mass matrices (including rotational inertia and extra non-structural mass) to support modal and eigenvalue analysis.
- **Post-Processing & Visualization:**
  - **Pandas Integration:** Convert displacements, reactions, and element forces directly into pandas DataFrames for easy analysis and post-processing.
  - **Graphical Plots:** Publication-ready visualization of undeformed and deformed shapes, realistic boundary supports (fixed clamp, pin, roller), load notation badges, and shaded distributed load bands using Matplotlib.

## Project Structure & Architecture

The repository is modularly architected into distinct layers: core finite element abstractions, element formulations, constitutive material models, cross-sections, numerical solvers, external schema adapters, and visualization utilities:

```text
fem-2d/
├── fem2d/                       # Core finite element analysis package
│   ├── structure.py             # Global model container, DOF indexing, matrix assembly, and solvers
│   ├── nodes.py                 # Node definitions, boundary conditions [ux, uy, rz], mass/inertia
│   ├── loads.py                 # Concentrated forces/moments, UDLs, and triangular/varying UVLs
│   ├── solver.py                # Newton-Raphson nonlinear solver for large displacements
│   ├── results.py               # Post-processing, pandas DataFrame extraction, and text/HTML reports
│   ├── buckling_analysis.py     # Elastic buckling analysis via geometric stiffness (Ke + λ Kg)
│   ├── m_phi_analysis.py        # Cross-section moment-curvature (M-φ) analysis
│   ├── elements/                # Structural element library
│   │   ├── element.py           # Abstract base element class (ElementBase)
│   │   ├── beam.py              # 2D beam (axial + bending + Timoshenko shear deformation, releases, mass)
│   │   ├── beam_hinges.py       # Beam element with moment releases (hinge_i, hinge_j)
│   │   ├── truss.py             # Linear and corotational geometrically non-linear truss element
│   │   ├── spring.py            # 2D translational elastic spring element
│   │   ├── beamNL.py            # Geometrically non-linear beam formulations
│   │   ├── beam_materialNL.py   # Material-nonlinear beam element
│   │   └── trussNL.py           # Specialized nonlinear truss formulation
│   ├── materials/               # Constitutive material relationships
│   │   ├── material.py          # Abstract material base class
│   │   ├── elastic.py           # Linear elastic material (ElasticMaterial)
│   │   ├── bilinear.py          # Elastoplastic bilinear material (BilinearMaterial)
│   │   └── csv_material.py      # Custom stress-strain curve loaded from CSV (CSVMaterial)
│   ├── sections/                # Geometric cross-section definitions
│   │   ├── section.py           # Standard cross-section (A, Iz)
│   │   ├── fiber.py             # Discretized fiber section for nonlinear stress integration
│   │   └── moment_curvature.py  # User-defined moment-curvature section
│   ├── adapters/                # Interoperability with the structural engineering ecosystem
│   │   └── struct_core_adapter.py # Two-way converter for the struct_core JSON schema
│   └── utils/                   # High-level utilities and visualization
│       ├── simple_frame.py      # Streamlined high-level builder API (SimpleFrame)
│       └── draw_structure.py    # Publication-quality structure visualizer (DrawStructure)
├── examples/                    # Executable scripts and benchmarks
│   ├── linear/                  # Classical textbook frames, continuous beams, and trusses
│   ├── non_linear/              # Corotational trusses, buckling factors, and nonlinear materials
│   ├── dynamic/                 # Mass matrix assembly, modal analysis, and time-history examples
│   └── adapters/                # Model conversion and round-tripping with struct_core
├── tests/                       # Comprehensive automated test suite (pytest)
└── docs/                        # Sphinx documentation configuration and guides
```

### Module Breakdown

| Module / Directory | Primary Role | Key Classes & Functions |
| --- | --- | --- |
| `fem2d.structure` | Central coordinator for the FE model | `Structure` |
| `fem2d.nodes` | Spatial points and boundary constraints | `Node` |
| `fem2d.loads` | Concentrated loads, UDLs, and UVLs | `PointLoad`, `DistributedLoad`, `TriangularLoad`, `ElementPointLoad` |
| `fem2d.elements` | 2D element stiffness and mass formulations | `BeamElement`, `BeamWithHingesElement`, `TrussElement`, `SpringElement` |
| `fem2d.materials` | Material constitutive laws | `ElasticMaterial`, `BilinearMaterial`, `CSVMaterial` |
| `fem2d.sections` | Cross-section geometry & fiber discretization | `Section`, `FiberSection`, `MomentCurvatureSection` |
| `fem2d.solver` | Nonlinear iterative solvers | `NewtonRaphsonSolver` |
| `fem2d.buckling_analysis` | Critical elastic buckling factors | `buckling_analysis` |
| `fem2d.results` | Tabular output and pandas integration | `Results` |
| `fem2d.adapters` | Ecosystem data interchange (`struct_core`) | `model_from_core`, `model_to_core`, `result_to_core` |
| `fem2d.utils.simple_frame` | High-level quick-modeling interface | `SimpleFrame` |
| `fem2d.utils.draw_structure` | Engineering visualization & figure export | `DrawStructure` |

## Installation

### From Source (Developer Install)

1. Clone the repository:
   ```bash
   git clone https://github.com/learnstructure/fem-2d.git
   cd fem-2d
   ```

2. Install in editable mode along with development dependencies:
   ```bash
   pip install -e .[dev]
   ```

## Quick Start Examples

### 1. Linear Static Frame Analysis (High-Level API)

The `SimpleFrame` class provides a simplified API for building and solving structures.

```python
from fem2d import SimpleFrame
from fem2d.results import Results

# Initialize simple frame with engineering units
frame = SimpleFrame(units="kip, in")

# Define nodes (id, x, y)
frame.add_node(1, 0.0, 0.0)
frame.add_node(2, 0.0, 120.0)
frame.add_node(3, 120.0, 120.0)
frame.add_node(4, 120.0, 0.0)

# Properties
E = 30000.0  # ksi
A = 10.0     # sq. in.
I = 200.0    # in^4

# Add frame elements (id, node_i, node_j, E, A, I)
frame.add_frame(1, 1, 2, E, A, I)
frame.add_frame(2, 2, 3, E, A, I / 2)
frame.add_frame(3, 3, 4, E, A, I)

# Apply fixed supports at base nodes (node_id, [ux, uy, rz])
frame.add_support(1, [True, True, True])
frame.add_support(4, [True, True, True])

# Apply nodal loads (node_id, [Fx, Fy, Mz])
frame.add_node_load(2, [10.0, 0.0, 0.0])
frame.add_node_load(3, [0.0, 0.0, 5.0])

# Solve the structure
frame.solve()

# Retrieve results
results = Results(frame)
print("Node Displacements:\n", results.node_displacements())
print("Reactions:\n", results.reactions())
print("Element End Forces:\n", results.element_forces())
```

### 2. Geometrically Non-Linear Truss Analysis & Self-Weight

For advanced analyses, use the core `Structure` class. The same
`TrussElement` is used for both linear and corotational geometrically
non-linear analyses — pass `geometric_nonlinear=True` to
`Structure.solve(...)` to switch on the corotational formulation.

```python
from fem2d import Structure, Node, ElasticMaterial, TrussElement
from fem2d.results import Results
from fem2d.sections import Section

# Create structure with canonical engineering units
structure = Structure(units="kN, m")
node1 = Node(1, 0.0, 0.0)
node2 = Node(2, 4.0, 3.0)
node3 = Node(3, 8.0, 0.0)

structure.add_node(node1)
structure.add_node(node2)
structure.add_node(node3)

# Material and section (G is automatically derived from E and nu; specify unit_weight)
E = 200e6              # Material modulus (kN/m^2)
EA = 45155.0           # axial stiffness (kN)
A = EA / E             # cross-sectional area (m^2)
material = ElasticMaterial(E=E, nu=0.30, unit_weight=77.0)  # steel: 77 kN/m^3
section = Section(A)

# Add linear truss elements (corotational path is enabled via solve below)
structure.add_element(TrussElement(1, node1, node2, material, section))
structure.add_element(TrussElement(2, node2, node3, material, section))
structure.add_element(TrussElement(3, node3, node1, material, section))

# Support boundaries — compact (ux, uy) form
node1.set_support(1, 1)  # pinned
node3.set_support(0, 1)  # roller

# External vertical point force at Node 2
node2.set_load(fx=0.0, fy=-2000.0, mz=0.0)

# Apply element self-weight directly (auto-derives gravity g = 9.80665 m/s^2 from units)
structure.apply_self_weight(direction="-y")

# Run Newton-Raphson analysis with corotational formulation
structure.solve(geometric_nonlinear=True, tolerance=1e-8, max_iter=30)

# Print displacements and forces
results = Results(structure)
print(results.node_displacements())
print(results.element_forces())
```

### 3. Visualizing Structures

`DrawStructure` renders publication-ready Matplotlib plots with realistic support symbols (fixed clamps, pins, rollers), clear load notation badges, shaded distributed load bands, internal hinge releases, and deformed shapes:

```python
from fem2d import DrawStructure

# Initialize plotter with analyzed structure (set displacement scale factor)
plotter = DrawStructure(structure, scale=50.0)

# Render structure in Matplotlib window (or save directly to file)
plotter.draw(
    show_deformed=True,
    show_loads=True,
    show_supports=True,
    show_node_labels=True,
    show_grid=False,             # Clean background without distracting gridlines
    support_style="detailed",    # "detailed" (fixed clamp, pin, roller) or "box"
    save_path="deformed_shape.png"  # Optional: export 300 DPI high-res figure
)
```

### 4. Member End Releases (Internal Hinges)

Internal hinges and releases can be applied directly to `BeamElement` at either end via boolean flags or 3-DOF release vectors `[u1, u2, theta]`:

```python
from fem2d import BeamElement, BeamWithHingesElement

# Moment release at end j (classical internal hinge)
elem1 = BeamElement(1, node1, node2, material, section, hinge_j=True)

# Moment release at start end i
elem2 = BeamElement(2, node2, node3, material, section, hinge_i=True)

# General 3-DOF release vector [u1, u2, theta] (consistent with fem3d)
elem3 = BeamElement(3, node3, node4, material, section, releases_j=[False, False, True])

# Backward compatibility: BeamWithHingesElement is retained and fully supported
elem4 = BeamWithHingesElement(4, node4, node5, material, section, hinge_j=True)
```

### 5. Timoshenko Beam Formulation (Shear Deformation)

`BeamElement` supports both classical Euler-Bernoulli beam theory and Timoshenko beam theory with shear deformation effects. When `include_shear_deformation=True` (or when shear area $A_s$ or cross-section depth $h$ and Poisson's ratio $\nu$ are provided), `fem2d` includes the transverse shear flexibility factor:

$$\beta_s = \frac{12 E I}{G A_s L^2}$$

where $G = \frac{E}{2(1+\nu)}$ is the shear modulus, $A_s = \kappa A$ is the effective shear area ($\kappa \approx 5/6$ for rectangular sections), $I$ is the moment of inertia, and $L$ is the member length.

The local transverse stiffness terms are modified accordingly:

$$k_{12} = \frac{12 E I}{L^3 (1 + \beta_s)}, \quad k_6 = \frac{6 E I}{L^2 (1 + \beta_s)}$$

$$k_4 = \frac{(4 + \beta_s) E I}{L (1 + \beta_s)}, \quad k_2 = \frac{(2 - \beta_s) E I}{L (1 + \beta_s)}$$

```python
from fem2d import Structure, Node, Section, ElasticMaterial, BeamElement

structure = Structure(units="kip, in")
node1 = Node(1, 0.0, 0.0)
node2 = Node(2, 60.0, 0.0)  # Short/deep beam (L = 60 in)
structure.add_node(node1)
structure.add_node(node2)

# Material with shear modulus G derived from E and nu
material = ElasticMaterial(E=29000.0, nu=0.3)

# Rectangular section with shear area As automatically computed: As = (5/6) * b * h
section = Section.from_rectangle(width=10.0, depth=25.0, include_shear_deformation=True)

# Define beam element with shear deformation enabled
elem = BeamElement(
    1, node1, node2,
    material=material,
    section=section,
    include_shear_deformation=True,
)
structure.add_element(elem)
node1.set_support(1, 1, 1)  # Fixed support
node2.set_load(fy=-50.0)    # 50 kip tip load
structure.solve()
```

For deep beams or short spans ($L/h < 5$), shear deformation can significantly increase deflections and redistribute internal forces compared to classical slender beam theory.

### 6. Dynamic & Modal Analysis Bridge with `structdyn`

`fem2d` seamlessly bridges into `structdyn` for eigenvalue modal analysis and transient earthquake time-history integration:

```python
from fem2d import SimpleFrame
from structdyn import MDF, GroundMotion, elcentro_chopra

# Build 2D frame with lumped story mass
frame = SimpleFrame()
# ... add nodes, frames, supports, and lumped masses ...

# Create dynamic MDF system with 5% modal damping on the first 3 modes
mdf = MDF.from_fem2d(frame.structure, zeta=0.05, n_modes=3)

# 1. Eigenvalue Modal Analysis
omega, phi = mdf.modal.modal_analysis(n_modes=3)
freqs = omega / (2.0 * 3.14159)
periods = 1.0 / freqs
print("Natural Periods (s):", periods)

# 2. Transient Earthquake Response History
elc = elcentro_chopra()
gm = GroundMotion.from_arrays(elc["acc (g)"], dt=0.02, scale_factor=386.09)
response = mdf.find_response_ground_motion(gm, direction="x", method="newmark_beta")
print("Peak Roof Displacement (in):", response["u5"].abs().max())

# 3. Direct Response Spectrum Analysis (ASCE 7 CQC)
params = SeismicParameters(sds=1.0, sd1=0.6, tl=8.0, site_class="D", damping_ratio=0.05)
rsa_res = frame.response_spectrum_analysis(
    params=params, direction="x", scale_factor=386.09, combination="CQC"
)
print("Combined Base Shear (kips):", rsa_res["combined_base_shear"][0])
```

### 7. Ecosystem Integration with `struct_core`

`fem2d` ships with a three-function adapter to the
[struct_core](https://github.com/learnstructure/struct_core) schema so that
models and analysis results can be shared with other packages in the
ecosystem (visualizers, design tools, code-checkers, …).

Install the optional dependency:

```bash
pip install fem2d[ecosystem]
```

The adapter exposes three functions — all are pure conversions, no side
effects, no global state. Naming follows `<artefact>_<direction>_<target>`:

| Function | Direction | Purpose |
| --- | --- | --- |
| `model_from_core` | `struct_core → fem2d` | Build a `fem2d.Structure` from a `struct_core.StructuralModel` (or `Project`), supporting individual load cases, factored load combinations, and automatic member self-weight calculation via `struct_load`. |
| `model_to_core`   | `fem2d → struct_core` | Build a `struct_core.StructuralModel` from a `fem2d.Structure` / `SimpleFrame` / `Results`. |
| `result_to_core`  | `fem2d → struct_core` | Build a `struct_core.AnalysisResult` from an analyzed `fem2d` model. |
| `solve_combinations_and_envelope` | `struct_core → fem2d → struct_core` | Solve all factored load combinations and compile the governing max/min force envelope for design. |

There is intentionally **no** `result_from_core` — round-tripping an
`AnalysisResult` back into a `fem2d.Results` is out of scope.

#### 4.0 Solve Factored Load Combinations & Force Envelopes

```python
from fem2d import solve_combinations_and_envelope

# Solves all combinations in project.model.load_combinations, stores individual
# combination results, and compiles an enveloped AnalysisResult with governing Pu, Mu, Vu
comb_results, envelope_result = solve_combinations_and_envelope(project, envelope_case_id="Envelope")
print("Governing envelope beam moment:", envelope_result.element_results[2].max_moment)
```

#### 4.1 Convert analysis results for downstream consumers

```python
from fem2d import SimpleFrame, result_to_core
from struct_core import save_json

frame = SimpleFrame()
frame.add_node(1, 0.0, 0.0)
frame.add_node(2, 120.0, 0.0)
frame.add_frame(1, 1, 2, 30000.0, 10.0, 200.0)
frame.add_support(1, [True, True, True])
frame.add_node_load(2, [10.0, 0.0, 0.0])
frame.solve()

# fem2d → struct_core.AnalysisResult
result = result_to_core(frame, analysis_case_id="Static")

# Serialize to JSON for downstream consumers
save_json(result, "static_result.json")
```

#### 4.2 Build a struct_core model from a fem2d `Structure`

`model_to_core` automatically maps the structure's `units` to the `StructuralModel.units` schema and places applied loads under `default_load_case_id="LIVE"` by default, keeping the `DEAD` case reserved for self-weight materialization:

```python
from fem2d import SimpleFrame, model_to_core, result_to_core
from struct_core import Project, save_json

frame = SimpleFrame(units="kip, in")
frame.add_node(1, 0.0, 0.0)
frame.add_node(2, 100.0, 0.0)
frame.add_frame(1, 1, 2, 30000.0, 10.0, 200.0)
frame.add_support(1, [True, True, True])
frame.add_node_load(2, [0.0, -1000.0, 0.0])
frame.solve()

# fem2d → struct_core.StructuralModel (preserves units and places loads in LIVE case)
sc_model = model_to_core(frame, default_load_case_id="LIVE")

# Wrap in a Project and serialize
project = Project(model=sc_model)
project.metadata.title = "Cantilever beam"
save_json(project, "model.json")
```

#### 4.3 Build a `fem2d.Structure` from a `struct_core` model

`model_from_core` instantiates a `fem2d.Structure` with the project's canonical `Units`, maps materials and section geometries, and optionally filters by `load_case_id` (or aggregates all active load cases by default):

```python
from fem2d import model_from_core

# `my_project` is a struct_core.Project (or StructuralModel)
# Loads can be aggregated from all cases, or selected specifically:
structure = model_from_core(my_project, load_case_id=None)

# Direct self-weight calculation (auto-computes gravity g from units)
structure.apply_self_weight(direction="-y")

structure.solve()
```

#### 4.4 Round-trip a model through `struct_core`

The conversion functions maintain strict numerical and referential equivalence:

```python
from fem2d import SimpleFrame, model_from_core, model_to_core, result_to_core

frame = SimpleFrame(units="kip, in")
# ... build & solve ...

# fem2d → core → fem2d
sc_model = model_to_core(frame)
round_tripped = model_from_core(sc_model)
round_tripped.solve()
assert round_tripped.disp == pytest.approx(frame.structure.disp)
```

## Running Tests

Verify that your installation is working correctly by running the tests:

```bash
pytest
```

## Citing FEM2D

If you use `fem2d` in your academic research or professional work, please cite it as follows:

```
Mandal, A. (2026). FEM2D: An open-source Python library for structural analysis of 2D structures (v0.5.0). Zenodo. https://doi.org/10.5281/zenodo.20990850
```

Refer to [CITATION.cff](CITATION.cff) for the BibTeX format details.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.