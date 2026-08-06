# Replication Materials

This folder contains public-facing replication materials for the
archive branch.

## Rendered Reports

- [Empirical FCI example](pdf/fci_example.pdf)
- [Simulation 1](pdf/sim1_replication.pdf)
- [Simulation 2](pdf/sim2_replication.pdf)

The empirical raw data cannot be redistributed, so the empirical QMD is
not included. The rendered empirical report is included to document the
analysis that was run locally. Self-contained HTML versions are also
included under `html/` for local viewing.

## Simulation QMDs

The simulation QMDs are in `qmd/`. They summarize saved condition-level
simulation results expected under `results/`. The large generated result
files are not included in the package build.

The installed package also provides `rdif_power_sim()` for users who want
to run smaller RDIF power simulations with the BH-start calibration used
in this archive branch.
