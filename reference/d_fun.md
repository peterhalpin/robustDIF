# The R-DIF scaling functions for item intercepts / thresholds.

Computes the scaling function to be used from the item thresholds (d) in
groups = {1, 2}. The options are:

- `type = 1`: computes `d.fun1 = (d2 - d1)/a1`

- `type = 2`: computes `d.fun2 = (d2 - d1)/a2`

- `type = 3`: computes `d.fun3 = (d2 - d1)/sqrt{(a1^2 + a2^2)/2}`

## Usage

``` r
d_fun(mle, type = 3)
```

## Arguments

- mle:

  the output of
  [`get_model_parms`](https://peterhalpin.github.io/robustDIF/reference/get_model_parms.md)

- type:

  a number in `1:3` indicating which version of delta to compute. See
  description for details.

## Value

The vector of scaling function values.

## Details

Computes `(d2 - d1)/a` for each threshold (d) of each item in groups g =
{1, 2}. The parameter 'a' depends on the value of `type`:
