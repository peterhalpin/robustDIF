# Example data set for R-DIF functions.

A named list containing the maximum likelihood estimates and their
estimated covariance matrix, for the 2PL IRT model fitted to 5 items in
two independent groups. The first item has additive bias of .5 applied
to the intercept only. The groups have a mean difference of .5 standard
deviations on the latent trait. The variances of the latent trait are
equal in each group.

## Usage

``` r
rdif.eg
```

## Format

A named list with 4 components:

- par0:

  A `data.frame` or named `list` with `par0$a` containing the item
  slopes and `par0$d` containing the item intercepts, for the reference
  group.

- par1:

  The item parameter estimates of the comparison groups. See `par0` for
  formatting.

- vcov0:

  The covariance matrix of `par0`, formatted as either a `data.frame` or
  `matrix`. The parameters should be organized by item, with the slope
  parameter coming first and the intercept parameter coming second
  (e.g., `a.item1, d.item1, a.item2, d.item2, ...`).

- vcov1:

  The covariance matrix of `par1`. See `vcov0` for formatting.
