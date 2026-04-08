# Compute starting values for [`rdif`](https://peterhalpin.github.io/robustDIF/reference/rdif.md).

Compute starting values for
[`rdif`](https://peterhalpin.github.io/robustDIF/reference/rdif.md).

## Usage

``` r
get_starts(mle, fun = "d_fun3", alpha = 0.05)
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

## Value

A vector containing the median of
[`y_fun`](https://peterhalpin.github.io/robustDIF/reference/y_fun.md),
the least trimmed squares estimate of location for
[`y_fun`](https://peterhalpin.github.io/robustDIF/reference/y_fun.md)
with 50-percent trim rate, and the minimum of
[`rho_grid`](https://peterhalpin.github.io/robustDIF/reference/rho_grid.md).
