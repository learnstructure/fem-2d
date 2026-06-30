FEM2D Documentation
===================

FEM2D – An open‑source Python library for 2‑D finite element analysis of structural frames.

It supports linear static analysis, non‑linear analysis, and dynamic analysis of 2‑D structures.


Features
--------

- **Element Library:**

  - **Beam Element:** Elastic 2D Euler-Bernoulli beam elements including axial, shear, and bending stiffness. Supports uniform and varying member loads, moment releases (hinges), and rotational/translational mass.
  - **Truss Element:** Pin-jointed bar elements with axial stiffness only.
  - **Spring Element:** 2D elastic spring elements with customizable axial stiffness.

- **Analysis Types:**

  - **Linear Static Analysis:** Standard matrix analysis under nodal loads, distributed loads, and concentrated member loads.
  - **Geometrically Non-Linear Analysis:** Iterative solver using the Newton-Raphson scheme combined with corotational formulations for large displacement/rotation problems.
  - **Mass Matrix Assembly:** Assembles global mass matrices (including rotational inertia and extra non-structural mass) to support modal and eigenvalue analysis.

- **Post-Processing & Visualization:**

  - **Pandas Integration:** Convert displacements, reactions, and element forces directly into pandas DataFrames for easy analysis and post-processing.
  - **Graphical Plots:** Plot the undeformed/deformed configurations, support conditions, and applied loads using Matplotlib.

Global and Local Axes
---------------------

**Global Axes:**
  X (+ve: right), Y (+ve: up), rotation (+ve: counter-clockwise)

**Local Axes:**
  x (+ve: from node i → node j), y (+ve: upwards when viewing element with left hand at node i and right hand at node j), rotation (+ve: counter-clockwise)

**Note:** Joint loads use global axes; member loads use local axes.

API Reference
-------------

.. toctree::
   :maxdepth: 2

   api
