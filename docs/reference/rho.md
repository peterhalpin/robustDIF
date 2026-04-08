# The bi-square rho function.

If `abs(u) > k` , `rho(u) = 1`. Else, `psi(u) = 1 - (1 - (u/k)^2)^3`.

## Usage

``` r
rho(u, k = 1.96)
```

## Arguments

- u:

  Can be a single value, vector, or matrix.

- k:

  The tuning parameter. Can be a scalar or the same dimension as `u`.

## Value

The bi-square rho function.
