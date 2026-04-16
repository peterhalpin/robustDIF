# Extract and format item parameter estimates and their covariance matrix

Takes a 1-factor model fit or list of 1-factor model fits from
[`mirt`](https://philchalmers.github.io/mirt/reference/mirt.html) or
[`cfa`](https://rdrr.io/pkg/lavaan/man/cfa.html) and formats the item
parameter estimates and their covariance matrix for use in other
`robustDIF` functions.

## Usage

``` r
get_model_parms(object, cluster = NULL)
```

## Arguments

- object:

  model fit from a multigroup analysis or list of model fits for each
  group for a 1-factor model. See Details.

- cluster:

  optional clustering IDs used to compute a cluster-robust covariance
  matrix for mirt fits. For `MultipleGroupClass`, provide an atomic
  vector of length equal to the number of observations. For a list of
  `SingleGroupClass` fits, provide a list of cluster-ID vectors (one per
  model).

## Value

A three-element `list`:

- `par.names`: list with `internal` and `original` parameter names.

- `est`: list (one element per group) of data frames containing item
  parameters by row (`a1`, `d1`, `d2`, ...).

- `var.cov`: list (one element per group) of covariance matrices for the
  corresponding parameter vectors.

## Details

The function takes a fitted 1-factor multigroup model or list of fitted
1-factor single group models. The factor must be standardized (i.e.,
variance = 1) and the covariance matrix be asymptotically correct.
Currently, the function accepts:

- a [`mirt`](https://philchalmers.github.io/mirt/reference/mirt.html)
  object of class `SingleGroupClass` or `MultipleGroupClass` with
  `SE = TRUE` (to return covariance matrix) and `itemtype` of any
  combination of `"2PL", "graded", or "gpcm"`.

- a `lavaan` object estimated from
  [`cfa`](https://rdrr.io/pkg/lavaan/man/cfa.html) with `std.lv = TRUE`.

When `cluster` is supplied for mirt fits, the covariance matrix is
computed using a cluster-robust sandwich estimator with Oakes bread and
cluster-summed empirical scores. A CR1 finite-sample correction is
applied: `(G/(G-1)) * ((N-1)/(N-p))`, where `G` is the number of
clusters, `N` is the number of observations, and `p` is the number of
free parameters. It is possible to use fits from other software with
`robustDIF` functions, but the parameter estimates and their covariance
matrices must be formatted identically to what is returned by
`get_model_parms`. For details, see the documentation for the example
dataset
[`rdif.eg`](https://peterhalpin.github.io/robustDIF/reference/rdif.eg.md).

## See also

[`rdif.eg`](https://peterhalpin.github.io/robustDIF/reference/rdif.eg.md)
