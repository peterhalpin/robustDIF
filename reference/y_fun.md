# R-DIF scaling functions for item intercepts/thresholds and slopes.

Computes the a scaling function using the item thresholds (d) and slopes
(a) in groups = {0, 1}. The options are:

- `"a_fun1"`: computes `a2/a1`

- `"a_fun2"`: computes `log(a2/a1)`

- `"d_fun1"`: computes `(d2 - d1)/a1` for each threshold

- `"d_fun2"`: computes `(d2 - d1)/a2` for each threshold

- `"d_fun3"`: computes `(d2 - d1)/sqrt((a1^2 + a2^2)/2)` for each
  threshold

## Usage

``` r
y_fun(mle, fun = "d_fun3")
```

## Arguments

- mle:

  the output of
  [`get_model_parms`](https://peterhalpin.github.io/robustDIF/reference/get_model_parms.md)

- fun:

  one of `c("a_fun1", "a_fun2", "d_fun1", "d_fun2", "d_fun3")`. See
  description for details.

## Value

A vector of scaling function values.

## Details

Computes the scaling function specified by `fun`, for each item.
