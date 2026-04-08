# The derivative of the bi-square psi function.

If `abs(u) > k` , `psi_prime(u) = 0`. Else,
`psi_prime(u) = (1 - (u/k)^2)^2 - (2u/k)^2 (1 - (u/k)^2)`.

## Usage

``` r
psi_prime(u, k = 1.96)
```

## Arguments

- u:

  Can be a single value, vector, or matrix.

- k:

  The tuning parameter. Can be a scalar or the same dimension as `u`.

## Value

The bi-square psi-prime function.
