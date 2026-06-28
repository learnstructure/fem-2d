# FEM2D
An open-source Python library for structural finite element analysis of 2D structures.

**FEM2D** is a Python package for performing 2D finite element analysis (FEA) of structural frames, including truss, beam, and spring elements with support for linear static analysis, geometrically non-linear analysis, and global mass matrix assembly for dynamic analysis.

Documentation: [![Documentation Status](https://readthedocs.org/projects/fem2d/badge/?version=latest)](https://fem2d.readthedocs.io/en/latest/)


## Features

- **Element Library:**
  - **Beam Element:** Elastic 2D Euler-Bernoulli beam elements including axial, shear, and bending stiffness. Supports uniform and varying member loads, moment releases (hinges), and rotational/translational mass.
  - **Truss Element:** Pin-jointed bar elements with axial stiffness only.
  - **Spring Element:** 2D elastic spring elements with customizable axial stiffness.
- **Analysis Types:**
  - **Linear Static Analysis:** Standard matrix analysis under nodal loads, distributed loads, and concentrated member loads.
  - **Geometrically Non-Linear Analysis:** Iterative solver using the Newton-Raphson scheme combined with corotational formulations for large displacement/rotation problems.
  - **Mass Matrix Assembly:** Assemblies global mass matrices (including rotational inertia and extra non-structural mass) to support modal and eigenvalue analysis.
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

For advanced analyses, use the core `Structure` class with non-linear element formulations.

```python
from fem2d import Structure, Node, ElasticMaterial
from fem2d.elements.trussNL import TrussElementNL
from fem2d.results import Results

# Create structure and nodes
structure = Structure()
node1 = Node(1, 0.0, 0.0)
node2 = Node(2, 4.0, 3.0)
node3 = Node(3, 8.0, 0.0)

structure.add_node(node1)
structure.add_node(node2)
structure.add_node(node3)

# Materials and sections
E = 200e6     # Material Modulus (kN/m^2)
Area = 0.00022577  # Cross-sectional area (m^2)
material = ElasticMaterial(E)

# Create non-linear truss elements and add to structure
structure.add_element(TrussElementNL(1, node1, node2, material, Area))
structure.add_element(TrussElementNL(2, node2, node3, material, Area))
structure.add_element(TrussElementNL(3, node3, node1, material, Area))

# Support boundaries
node1.set_support(ux_fixed=True, uy_fixed=True)   # Pinned support
node3.set_support(ux_fixed=False, uy_fixed=True)  # Roller support

# External vertical point force at Node 2
node2.set_load(fx=0.0, fy=-2000.0, mz=0.0)

# Run non-linear Newton-Raphson analysis
structure.solve_nonlinear(tolerance=1e-8, max_iter=30)

# Print displacements
results = Results(structure)
print(results.node_displacements())
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

## Running Tests

Verify that your installation is working correctly by running the tests:

```bash
pytest
```

## Citing FEM2D

If you use `fem2d` in your academic research or professional work, please cite it as follows:

```
Mandal, A. (2026). FEM2D: An open-source Python library for structural analysis of 2D structures (Version 0.2.1) [Computer software]. https://github.com/learnstructure/fem-2d.git
```

Refer to [CITATION.cff](CITATION.cff) for the BibTeX format details.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.