# robustDIF
Peter F Halpin
2026-01-22

# robustDIF

An R package that provides functions for assessing differential item
functioning (DIF) in item response theory (IRT) models using methods
from robust statistics. Based on the paper: 

Halpin, P.F. (2022) Differential Item Functioning Via Robust Scaling. Arxiv Preprint. [https://arxiv.org/abs/2207.04598](https://arxiv.org/abs/2207.04598). Published in Psychometrika in 2024 under the same title.  

Please see the `README.html` file for a more detailed discussion of the methodology and recent updates. 

# Installation

```{r, eval = FALSE}
install.packages("remotes")
remotes::install_github("peterhalpin/robustDIF")
library(robustDIF)
```

```{r, echo = F}
#library(robustDIF)
```

# Example Dataset

The main user-facing functions are illustrated below using the built-in
example dataset `rdif.eg`. In the example dataset, there are a total of
five binary items and the first item has DIF on the intercept
(threshold) equal to 1/2 SD on the latent trait. The latent trait was
generated from *N*(0,1) in the reference group and *N*(.5, 1) in the
comparison group.

Check out the documentation for `rdif.eg` and `get_model_parms` for more
info about how to format data for use with `robustDIF`.

# Choice of scaling functions 

The RDIF procedure is a robust version of "mean/mean" IRT-based scaling
(see Kolen & Brennan, 2014, chap. 6). Robustness is achieved by
down-weighting items that exhibit DIF. Following estimation of the
scaling parameter, Wald tests of DIF can be implemented. The type of DIF
that is tested depends on the choice of scaling parameter.

Let `theta ~ N(mu, sigma)` in groups 1 and 2 and let `a_i` and `d_i` 
denote the item slope and intercept/thresholds respectively. The choices
of scaling parameters and corresponding item-level scaling functions are
as follows. 

| Name     | Scaling parameter                                   | Item-level scaling function                                              | Test for  DIF on      |
|----------|----------------------------------------------------|---------------------------------------------------------------------------|-----------------|
| `a_fun1` | `sigma2 / sigma1`                                  | `a2_i / a1_i`                                                             | slopes       |
| `a_fun2` | `log(sigma2 / sigma1)`                             | `log(a2_i / a1_i)`                                                       | slopes       |
| `d_fun1` | `(mu2 - mu1) / sigma1`                             | `(d2_i - d1_i) / a1_i`                                  | intercepts   |
| `d_fun2` | `(mu2 - mu1) / sigma2`                             | `(d2_i - d1_i) / a2_i`                                 | intercepts   |
| `d_fun3` | `(mu2 - mu1) / sqrt((sigma1^2 + sigma2^2) / 2)`   | `(d2_i - d1_i) / sqrt((a1_i^2 + a2_i^2) / 2)`            | slope  + intercept |



The illustration uses `d_fun3`. Equality between the scaling parameter and the item-level scaling function requires assuming the item does not exhibit DIF on either the slope or intercept, although the function is more sensitive to intercept DIF than slope DIF. 

## Robust Scaling

```{r}
rdif(mle = rdif.eg, fun = "d_fun3", alpha = .05)
```

The estimated scaling parameter is equal to approximately 0.38. We can see that the first item was downweighted to zero when estimating the scaling parameter. Note that the standard error of the scaling parameter is not provided by `rdif` but is available via `delta_test`. 

## Wald Tests of DIF

Inferences about DIF can be made via Wald tests of whether the item-level scaling function differs from the estimated scaling parameter. 

```{r}
dif_test(mle = rdif.eg, fun = "d_fun3")
```

In this output `delta` is the estimated scaling parameter subtracted from the item-level scaling function.  

## Testing Whether DIF Affects Conclusions about Impact

One may wish to know whether any DIF, if present, affects conclusions
about impact (i.e., how the groups differ on the latent trait). 
One way to approach this question is to compare a naive estimate of impact that ignores DIF to 
an estimator that is robust to DIF. In `robustDIF`, the naive estimator is chosen to be the
unweighted mean of the item-level functions (i.e., the usual IRT moment-based scaling procedure). The robust estimator is as above.  The null hypothesis that both estimators are
consistent for the "true" scaling parameter leads to the following test.

```{r}
delta_test(mle = rdif.eg, fun = "d_fun3")
```

In this output delta = `naive.est` - `rdif.est`

## The Rho Function

For data analyses, it is useful to check whether the robust estimator
has a clear global minimum before proceeding to make inferences about DIF / Impact.  
In `robustDIF`, the minimization problem is described by `rho_fun`.

```{r}
rho <- rho_grid(mle = rdif.eg, fun = "d_fun3", grid.width = .01)
plot(rho$theta, rho$rho, type = "l")
```

Note the minimizing values of theta in the plots corresponds to the
values reported by `rdif` above. 
