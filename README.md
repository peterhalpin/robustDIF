
# robustDIF

Functions for assessing differential item functioning (DIF) in item
response theory (IRT) models using methods from robust statistics. Based
on the paper: Halpin, P.F. (2022) Differential Item Functioning Via
Robust Scaling.

Currently `robustDIF` only supports the two-parameter logistc (2PL) IRT
model in two independent groups.

## Installation

``` r
# Development version
install.packages("devtools")
devtools::install_github("peterhalpin/robustDIF")
library(robustDIF)
```

## Example

The main user-facing functions are illustrated below using the built-in
example data set `rdif.eg`. Check out the documentation for `rdif.eg`
and `get_irt_pars` for how to format your own data for use with
`robustDIF`

This is an R Markdown document. Markdown is a simple formatting syntax
for authoring HTML, PDF, and MS Word documents. For more details on
using R Markdown see <http://rmarkdown.rstudio.com>.

When you click the **Knit** button a document will be generated that
includes both content as well as the output of any embedded R code
chunks within the document. You can embed an R code chunk like this:

``` r
summary(cars)
```

    ##      speed           dist       
    ##  Min.   : 4.0   Min.   :  2.00  
    ##  1st Qu.:12.0   1st Qu.: 26.00  
    ##  Median :15.0   Median : 36.00  
    ##  Mean   :15.4   Mean   : 42.98  
    ##  3rd Qu.:19.0   3rd Qu.: 56.00  
    ##  Max.   :25.0   Max.   :120.00

## Including Plots

You can also embed plots, for example:

![](README_files/figure-gfm/pressure-1.png)<!-- -->

Note that the `echo = FALSE` parameter was added to the code chunk to
prevent printing of the R code that generated the plot.
