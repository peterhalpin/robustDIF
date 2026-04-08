# The gradient matrix of [`d_fun`](https://peterhalpin.github.io/robustDIF/reference/d_fun.md).

The gradient is taken with respect to the item parameters and organized
to be conformable with `Matrix::bdiag(mle$var.cov)`. When evaluating the
gradient under the null hypothesis of no DIF, the optional argument
`theta` can be provided. It replaces the item-specific values of d_fun
in the gradient computation.

## Usage

``` r
grad_d(mle, theta = NULL, type = 3)
```

## Arguments

- mle:

  the output of
  [`get_model_parms`](https://peterhalpin.github.io/robustDIF/reference/get_model_parms.md)

- theta:

  (optional) the scaling parameter. Replaces item-specific values of
  d_fun if provided.

- type:

  a number in `1:3` indicating which version of delta to compute. See
  description for details.

## Value

A matrix in which the columns are the gradient vectors of
[`d_fun`](https://peterhalpin.github.io/robustDIF/reference/d_fun.md),
for each item and threshold.

## See also

[`d_fun`](https://peterhalpin.github.io/robustDIF/reference/d_fun.md)
