# robustDIF

Rendered online supplementary materials are available here as PDFs:

- [Empirical FCI example](replication/pdf/fci_example.pdf)
- [Simulation 1](replication/pdf/sim1_replication.pdf)
- [Simulation 2](replication/pdf/sim2_replication.pdf)

The empirical raw data are not included in this repository, so the
empirical QMD is not distributed. Simulation QMDs are included under
`replication/qmd/`. Self-contained HTML versions are also included under
`replication/html/` for local viewing.

# Installation

```{r, eval = FALSE}
## Install this archive branch from the review repository.
## remotes::install_github("<owner>/<repo>@<branch>")
library(robustDIF)
```

# Simulation Power

This archive branch includes a small user-facing simulation helper for
planning RDIF power studies:

```{r, eval = FALSE}
out <- rdif_power_sim(
  n.reps = 100,
  n.persons = 500,
  n.items = 20,
  n.dif = 4,
  item.delta = 0.5,
  slope.multiplier = 1.5,
  alpha.adjust = "BH",
  seed = 20260724)

out$summary
```
