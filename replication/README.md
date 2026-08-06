# Replication Materials

This folder contains public-facing replication materials for the
archive branch.

## Rendered Reports

- [Empirical FCI example](html/fci_example.html)
- [Simulation 1](html/sim1_replication.html)
- [Simulation 2](html/sim2_replication.html)

The empirical raw data cannot be redistributed, so the empirical QMD is
not included. The rendered empirical report is included to document the
analysis that was run locally.

## Simulation QMDs

The simulation QMDs are in `qmd/`. They summarize saved condition-level
simulation results expected under `results/`. The large generated result
files are not included in the package build.

The installed package also provides `rdif_power_sim()` for users who want
to run smaller RDIF power simulations with the BH-start calibration used
in this archive branch.
