# Helper function used to format parameters estimates

Helper function used to format parameters estimates

## Usage

``` r
format_pars(pars, names.vec, type)
```

## Arguments

- pars:

  numeric vector of item parameter estimates

- names.vec:

  character vector item names

- type:

  character; are `pars` from `"lavaan"` or `"mplus"`?

## Value

data.frame of item parameter estimates

## See also

\[robustDIF::get_lavaan_params()\], \[robustDIF::get_mplus_params()\]
