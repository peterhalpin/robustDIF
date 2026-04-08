# Wald test of differential test functioning.

A Wald test of the difference between the unweighted mean of the
[`y_fun`](https://peterhalpin.github.io/robustDIF/reference/y_fun.md)
and robust scaling parameter from
[`rdif`](https://peterhalpin.github.io/robustDIF/reference/rdif.md).
Called internally by
[`rdif`](https://peterhalpin.github.io/robustDIF/reference/rdif.md)

## Usage

``` r
delta_test(object, theta = NULL, k = NULL, fun = "d_fun3")
```

## Arguments

- object:

  either the output of
  [`get_model_parms`](https://peterhalpin.github.io/robustDIF/reference/get_model_parms.md)
  or an `rdif` object from
  [`rdif`](https://peterhalpin.github.io/robustDIF/reference/rdif.md).

- theta:

  the estimated scaling parameter from
  [`rdif`](https://peterhalpin.github.io/robustDIF/reference/rdif.md).
  Not needed when `object` is an `rdif` object.

- k:

  the tuning parameter from
  [`rdif`](https://peterhalpin.github.io/robustDIF/reference/rdif.md).
  Not needed when `object` is an `rdif` object.

- fun:

  one of `c("a_fun1", "a_fun2", "d_fun1", "d_fun2", "d_fun3")`. See
  description for details.

## Value

A data.frame that contains the output of the test.

## Examples

``` r
#
if (FALSE) { # \dontrun{
mod <- rdif(mle = rdif.eg)
delta_test(object = rdif.eg, theta = mod$est, k = mod$k)
delta_test(mod)
} # }
```
