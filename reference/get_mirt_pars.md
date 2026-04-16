# Extract item parameter estimates and their covariance matrix from [`mirt`](https://philchalmers.github.io/mirt/reference/mirt.html).

Extract item parameter estimates and their covariance matrix from
[`mirt`](https://philchalmers.github.io/mirt/reference/mirt.html).

## Usage

``` r
get_mirt_pars(mirt.object, cluster = NULL)
```

## Arguments

- mirt.object:

  a [`mirt`](https://philchalmers.github.io/mirt/reference/mirt.html)
  object of class `SingleGroupClass` or `MultipleGroupClass`. Expected
  to be a 1-factor model with `SE = TRUE` and `itemtype` of any
  combination of `"2PL", "graded", or "gpcm"`.

- cluster:

  optional cluster-ID vector used to compute a cluster-robust covariance
  matrix using Oakes bread, empirical score outer products aggregated by
  cluster, and a CR1 finite-sample correction.

## Value

A three-element `list`:

- vector of parameter names taking the form "item.parameter"

- list (one element per group) of vectors of item parameter estimates

- list (one element per group) of covariance matrices of item parameter
  estimates
