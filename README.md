# Flexible-Manipulator

Recursive modeling of a flexible manipulator undergoing small deformations.

## Derivative comparisons

With CasADi on the MATLAB path, run either script. Each compares the position
and velocity Jacobians with automatic differentiation at one numeric state,
prints absolute/relative errors, and checks the relative errors against 1e-9.
There are no timing benchmarks.

- `src/tests/id_d_test.m`: `ID_flex_d` versus CasADi AD of `ID_flex`.
- `src/tests/id_FO_test.m`: `ID_flex_FO` versus CasADi AD of `ID_flex`.

## Simulation examples

- `src/tests/dynamics/test_FD_flex_7dof_grasp_cb_beam.m`: seven-DOF arm grasping a flexible beam.
- `src/tests/dynamics/flexible_base_manipulator.m`: floating rigid base, flexible beam, and manipulator.
- `src/tests/dynamics/test_7dof_all_flexible_exudyn.m`: seven-DOF arm with all flexible links; see the adjacent `README_7dof_all_flexible.md`.

The matching reference scripts are in `refs/`. Generate reference data before
running examples that require it. All generated data, figures, temporary files,
and other local tests are excluded from Git. Obsolete local tests and generated results have been removed. Shared test helpers remain in
`src/tests/dynamics/support/`.
