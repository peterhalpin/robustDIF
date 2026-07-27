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
# \donttest{
mod <- rdif(mle = rdif.eg)
dif_test(object = rdif.eg, theta = mod$est)
#>                 delta         se      z.test        p.val
#> item1_d1  0.351129635 0.09867331  3.55850682 0.0003729691
#> item2_d1  0.006566087 0.06654194  0.09867593 0.9213955818
#> item3_d1  0.017345461 0.10887831  0.15931052 0.8734242308
#> item4_d1 -0.006742385 0.20492283 -0.03290207 0.9737526822
#> item5_d1 -0.357290253 0.24134142 -1.48043486 0.1387572332
dif_test(mod)
#>                 delta         se      z.test        p.val
#> item1_d1  0.351129635 0.09867331  3.55850682 0.0003729691
#> item2_d1  0.006566087 0.06654194  0.09867593 0.9213955818
#> item3_d1  0.017345461 0.10887831  0.15931052 0.8734242308
#> item4_d1 -0.006742385 0.20492283 -0.03290207 0.9737526822
#> item5_d1 -0.357290253 0.24134142 -1.48043486 0.1387572332
# }
```
