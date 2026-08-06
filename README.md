# robustDIF

Rendered online supplementary materials are available here:

- [Empirical FCI example](https://htmlpreview.github.io/?https://github.com/peterhalpin/robustDIF/blob/Halpin2024-archive/replication/html/fci_example.html)
- [Simulation 1](https://htmlpreview.github.io/?https://github.com/peterhalpin/robustDIF/blob/Halpin2024-archive/replication/html/sim1_replication.html)
- [Simulation 2](https://htmlpreview.github.io/?https://github.com/peterhalpin/robustDIF/blob/Halpin2024-archive/replication/html/sim2_replication.html)

The empirical raw data are not included in this repository, so the
empirical QMD is not distributed. Simulation QMDs are included under
`replication/qmd/`.

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
