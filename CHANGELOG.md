# Changelog

All notable changes to `fem2d` are recorded in this file. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project uses
[semantic versioning](https://semver.org/). The code of the paper is the git tag `v2.1-paper`.

## [Unreleased]

### Added
- `model_from_core` refuses models with shell elements (`tests/test_shell_refusal.py`): fem2d is a plane-frame solver.
- `fem2d.loads.bending_shape_functions`: the exact interpolation (deflection and cross-section rotation) of a
  two-node Timoshenko beam, which reduces to the Hermite cubics without shear deformation.
- `BeamElement.shear_parameter()`: the shear deformation constant beta_s, used by both the stiffness matrix and the
  equivalent nodal loads.
- `tests/test_timoshenko_member_loads.py`: the equivalent nodal loads of `DistributedLoad`, `TriangularLoad`, and
  `ElementPointLoad` equal the fixed-end forces of the flexibility (unit-load) method, with and without shear
  deformation, and a member point load gives the same solution as a mesh with a node at the load, with and without an
  end release.

### Changed
- `model_from_core` applies the struct_core stiffness modifier keys only (`A`, `I3`, `As2`; the out-of-plane keys
  do not apply in 2D). The non-canonical keys `I` and `As` are no longer read; struct_core now refuses them.

### Fixed
- Member loads on Timoshenko beams (`include_shear_deformation=True`): `DistributedLoad`, `TriangularLoad`, and
  `ElementPointLoad` computed their equivalent nodal loads with Euler-Bernoulli shape functions, exact only for
  full-span uniform loads. They now use the exact Timoshenko shapes. Kassimali (2021), Matrix Analysis of Structures,
  Section 9.7, notes that the usual fixed-end force expressions exclude shear deformation. Results of Euler-Bernoulli
  beams and of full-span uniform loads do not change.
- Member point moments (`ElementPointLoad(mz=...)`) were split between the end rotations by the lever rule (M b / L
  and M a / L) with no end shears. They now get the correct fixed-end forces, for example -M/4 at each end rotation and
  -/+ 1.5 M / L at the end translations for a moment at midspan.
- `model_from_core` applies the element `mass` modifier (it was ignored) and no longer lets the `A` modifier scale
  the mass. As in SAP2000 27.1, measured in the struct_link tests, the A modifier scales the axial stiffness only and
  the mass modifier the whole member mass.
- `TriangularLoad` now gives its equivalent nodal loads with the same quadrature as `DistributedLoad` (its closed form
  was Euler-Bernoulli only).
