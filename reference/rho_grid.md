# Compute a grid of bi-square Rho values

Computes the objective function of the bi-square minimization problem in
a location parameter, theta. The theta values are obtained internally by
a grid search over the range of
[`y_fun`](https://peterhalpin.github.io/robustDIF/reference/y_fun.md).
Used for starting values and graphically diagnosing local solutions.

## Usage

``` r
rho_grid(mle, fun = "d_fun3", alpha = 0.05, grid.width = 0.01)
```

## Arguments

- mle:

  the output of
  [`get_model_parms`](https://peterhalpin.github.io/robustDIF/reference/get_model_parms.md)

- fun:

  one of `c("a_fun1", "a_fun2", "d_fun1", "d_fun2", "d_fun3")`. See
  description for details.

- alpha:

  the desired false positive rate for flagging items with DIF.

- grid.width:

  the width of grid points.

## Value

A named list with theta values and the corresponding Rho values.
