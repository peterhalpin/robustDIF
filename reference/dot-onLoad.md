# Registers S3 methods at load time: - print for class "rdif" - plot for class "rdif" - summary for class "rdif"

Registers S3 methods at load time: - print for class "rdif" - plot for
class "rdif" - summary for class "rdif"

## Usage

``` r
.onLoad(libname, pkgname)
```

## Arguments

- libname:

  Character string with the path to the package library.

- pkgname:

  Character string with the package name.

## Value

No return value, called for side effects when the package loads.
