# Extract item parameter estimates and their covariance matrix from `lavaan`.

Extract item parameter estimates and their covariance matrix from
`lavaan`.

## Usage

``` r
get_lavaan_pars(lavaan.object)
```

## Arguments

- lavaan.object:

  an object of class `lavaan`. Expected to be a 1-factor model estimated
  with `std.lv = TRUE`.

## Value

A three-element `list`:

- vector of parameter names taking the form "item.parameter"

- list (one element per group) of vectors of item parameter estimates

- list (one element per group) of covariance matrices of item parameter
  estimates
