# Estimate IRT scale parameters using the robust DIF procedure.

Estimation can be performed using iteratively re-weighted least squares
(IRLS) or Newton-Raphson (NR). Currently, only IRLS is implemented. If
`starting.value = "all"`, three starting values are computed: the median
of
[`y_fun`](https://peterhalpin.github.io/robustDIF/reference/y_fun.md),
the least trimmed squares estimate of location for
[`y_fun`](https://peterhalpin.github.io/robustDIF/reference/y_fun.md)
with 50-percent trim rate, and the minimum of
[`rho_grid`](https://peterhalpin.github.io/robustDIF/reference/rho_grid.md).
The estimate is computed from each starting value, and the solution with
the lowest value of the bi-square objective function is returned. If
there are multiple solutions, they are stored `other.solutions`.

## Usage

``` r
rdif(
  mle,
  fun = "d_fun3",
  alpha = 0.05,
  starting.value = "all",
  tol = 1e-07,
  maxit = 100,
  method = "irls"
)
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

- starting.value:

  one of `c("med", "lts", "min_rho", "all")` or a numerical value to be
  used as the starting value. See description for details.

- tol:

  convergence criterion for comparing subsequent values of estimate

- maxit:

  maximum number of iterations

- method:

  one of `c("irls", "newton")`. Currently, only IRLS is implemented.

## Value

An `rdif` object.

## Details

Implements M-estimation of an IRT scale parameter using the bi-square
loss function. Also returns the bi-square weights for each item.

## Examples

``` r
# Item intercepts, using the built-in example dataset "rdif.eg"
if (FALSE) rdif(mle = rdif.eg, fun = "d_fun3") # \dontrun{}

# Item slopes
if (FALSE) rdif(mle = rdif.eg, fun = "a_fun1") # \dontrun{}
```
