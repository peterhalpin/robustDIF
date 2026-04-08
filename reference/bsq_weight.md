# The bi-square weight function.

If `abs(u) > k` , `psi(u) = 0`. Else, `psi(u) = (1 - (u/k)^2)^2`.

## Usage

``` r
bsq_weight(u, k = 1.96)
```

## Arguments

- u:

  Can be a single value, vector, or matrix.

- k:

  The tuning parameter. Can be a scalar or the same dimension as `u`.

## Value

The bi-square psi-prime function.
