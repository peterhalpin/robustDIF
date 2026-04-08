# The R-DIF scaling function for item slopes.

Computes the scaling function `a2/a1` for item slopes (a) in groups g =
{1, 2}

## Usage

``` r
a_fun(mle, log = F)
```

## Arguments

- mle:

  the output of
  [`get_model_parms`](https://peterhalpin.github.io/robustDIF/reference/get_model_parms.md)

- log:

  logical: return of log(a2/a1)?

## Value

The vector of scaling function values.
