# Wald test of differential test functioning.

A Wald test of the difference between the unweighted mean of the
[`y_fun`](https://peterhalpin.github.io/robustDIF/reference/y_fun.md)
computed for all items, and the unweighted mean excluding `dif.items`.

## Usage

``` r
delta_test_from_dif(mle, dif.items, fun = "d_fun3")
```

## Arguments

- mle:

  the output of
  [`get_model_parms`](https://peterhalpin.github.io/robustDIF/reference/get_model_parms.md)

- dif.items:

  the indices of the items with DIF.

- fun:

  one of `c("a_fun1", "a_fun2", "d_fun1", "d_fun2", "d_fun3")`. See
  description for details.

## Value

A data.frame that contains the output of the test.
