# The gradient matrix of [`a_fun`](https://peterhalpin.github.io/robustDIF/reference/a_fun.md).

The gradient is taken with respect to the item parameters and organized
to be conformable with `Matrix::bdiag(mle$var.cov)`. When evaluating the
gradient under the null hypothesis of no DIF, the optional argument
`theta` can be provided. It replaces the item-specific values of a_fun
in the gradient computation.

## Usage

``` r
grad_a(mle, theta = NULL, log = F)
```

## Arguments

- mle:

  the output of
  [`get_model_parms`](https://peterhalpin.github.io/robustDIF/reference/get_model_parms.md)

- theta:

  (optional) the scaling parameter. Replaces item-specific values of
  alpha if provided.

- log:

  logical: return of log(a2/a1)?

## Value

A matrix in which the columns are the gradient vectors of
[`a_fun`](https://peterhalpin.github.io/robustDIF/reference/a_fun.md),
for each item.

## See also

[`a_fun`](https://peterhalpin.github.io/robustDIF/reference/a_fun.md)
