# The covariance matrix of IRT scaling functions.

When evaluating the covariance matrix under the null hypothesis of no
DIF, the optional argument `theta` can be provided. It replaces the
item-specific scaling functions in the gradient computation. Type should
be the same as used in
[`y_fun`](https://peterhalpin.github.io/robustDIF/reference/y_fun.md).

## Usage

``` r
vcov_y(mle, theta = NULL, fun = "d_fun3")
```

## Arguments

- mle:

  the output of
  [`get_model_parms`](https://peterhalpin.github.io/robustDIF/reference/get_model_parms.md)

- theta:

  (optional) the scaling parameter. Replaces item-specific scaling
  functions if provided.

- fun:

  one of `c("a_fun1", "a_fun2", "d_fun1", "d_fun2", "d_fun3")`. See
  description for details.

## Value

The covariance matrix of `y_fun`.

## See also

[`y_fun`](https://peterhalpin.github.io/robustDIF/reference/y_fun.md)
