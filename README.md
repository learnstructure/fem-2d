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
  - **Beam Element:** Elastic 2D Euler-Bernoulli beam elements including axial and bending stiffness (shear deformation is neglected). Supports uniform and varying member loads, moment releases (hinges), and rotational/translational mass.
  - **Truss Element:** Pin-jointed bar elements with axial stiffness only. The same `TrussElement` class is used for linear and corotational geometrically-nonlinear analyses.
  - **Spring Element:** 2D elastic spring elements with customizable axial stiffness.
- **Analysis Types:**
  - **Linear Static Analysis:** Standard matrix analysis under nodal loads, distributed loads, and concentrated member loads.
  - **Geometrically Non-Linear Analysis:** Iterative solver using the Newton-Raphson scheme combined with corotational formulations for large displacement/rotation problems. Enabled by passing `geometric_nonlinear=True` to `Structure.solve(...)`.
  - **Mass Matrix Assembly:** Assembles global mass matrices (including rotational inertia and extra non-structural mass) to support modal and eigenvalue analysis.
- **Post-Processing & Visualization:**
  - **Pandas Integration:** Convert displacements, reactions, and element forces directly into pandas DataFrames for easy analysis and post-processing.
  - **Graphical Plots:** Plot the undeformed/deformed configurations, support conditions, and applied loads using Matplotlib.

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

# Initialize simple frame
frame = SimpleFrame()

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

### 2. Geometrically Non-Linear Truss Analysis

For advanced analyses, use the core `Structure` class. The same
`TrussElement` is used for both linear and corotational geometrically
non-linear analyses — pass `geometric_nonlinear=True` to
`Structure.solve(...)` to switch on the corotational formulation.

```python
from fem2d import Structure, Node, ElasticMaterial, TrussElement
from fem2d.results import Results
from fem2d.sections import Section

# Create structure and nodes
structure = Structure()
node1 = Node(1, 0.0, 0.0)
node2 = Node(2, 4.0, 3.0)
node3 = Node(3, 8.0, 0.0)

structure.add_node(node1)
structure.add_node(node2)
structure.add_node(node3)

# Material and section
E = 200e6              # Material modulus (kN/m^2)
EA = 45155.0           # axial stiffness (kN)
A = EA / E             # cross-sectional area (m^2)
material = ElasticMaterial(E)
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

# Run Newton-Raphson analysis with corotational formulation
structure.solve(geometric_nonlinear=True, tolerance=1e-8, max_iter=30)

# Print displacements and forces
results = Results(structure)
print(results.node_displacements())
print(results.element_forces())
```

### 3. Visualizing Structures

You can easily generate visualization plots of undeformed and deformed structural shapes:

```python
from fem2d import DrawStructure

# Initialize plotter with analyzed structure
plotter = DrawStructure(structure, scale=0.05)

# Render structure in Matplotlib window
plotter.draw()
```

### 4. Ecosystem Integration with `struct_core`

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
| `model_from_core` | `struct_core → fem2d` | Build a `fem2d.Structure` from a `struct_core.StructuralModel` (or `Project`). |
| `model_to_core`   | `fem2d → struct_core` | Build a `struct_core.StructuralModel` from a `fem2d.Structure` / `SimpleFrame` / `Results`. |
| `result_to_core`  | `fem2d → struct_core` | Build a `struct_core.AnalysisResult` from an analyzed `fem2d` model. |

There is intentionally **no** `result_from_core` — round-tripping an
`AnalysisResult` back into a `fem2d.Results` is out of scope.

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

```python
from fem2d import SimpleFrame, model_to_core, result_to_core
from struct_core import Project, save_json

frame = SimpleFrame()
frame.add_node(1, 0.0, 0.0)
frame.add_node(2, 100.0, 0.0)
frame.add_frame(1, 1, 2, 30000.0, 10.0, 200.0)
frame.add_support(1, [True, True, True])
frame.add_node_load(2, [0.0, -1000.0, 0.0])
frame.solve()

# fem2d → struct_core.StructuralModel
sc_model = model_to_core(frame)

# Wrap in a Project and serialize
project = Project()
project.metadata.title = "Cantilever beam"
project.model = sc_model
save_json(project, "model.json")
```

#### 4.3 Build a `fem2d.Structure` from a `struct_core` model

```python
from fem2d import model_from_core

# `my_project` is a struct_core.Project (or StructuralModel)
structure = model_from_core(my_project)
structure.solve()
```

#### 4.4 Round-trip a model through `struct_core`

The three functions are designed to be inverse:

```python
from fem2d import SimpleFrame, model_from_core, model_to_core, result_to_core

frame = SimpleFrame()
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
Mandal, A. (2026). FEM2D: An open-source Python library for structural analysis of 2D structures (v0.4.0). Zenodo. https://doi.org/10.5281/zenodo.20990850
```

Refer to [CITATION.cff](CITATION.cff) for the BibTeX format details.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.