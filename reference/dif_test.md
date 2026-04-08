# Wald tests of differential item functioning (DIF).

A Wald test of DIF on each item. Called internally by
[`rdif`](https://peterhalpin.github.io/robustDIF/reference/rdif.md)

## Usage

``` r
dif_test(object, theta = NULL, fun = "d_fun3")
```

## Arguments

- object:

  either the output of
  [`get_model_parms`](https://peterhalpin.github.io/robustDIF/reference/get_model_parms.md)
  or an `rdif` object from
  [`rdif`](https://peterhalpin.github.io/robustDIF/reference/rdif.md).

- theta:

  the estimated scaling parameter from
  [`rdif`](https://peterhalpin.github.io/robustDIF/reference/rdif.md)

- fun:

  one of `c("a_fun1", "a_fun2", "d_fun1", "d_fun2", "d_fun3")`. See
  description for details.

## Value

A data.frame whose rows containing the results of the test for each item
parameter.

## Examples

``` r
if (FALSE) { # \dontrun{
mod <- rdif(mle = rdif.eg)
dif_test(object = rdif.eg, theta = mod$est)
dif_test(mod)
} # }
```
