# Testing for impact

## Differential item functioning

In the Item Response Theory (IRT) literature, Differential Item
Functioning (DIF) is an approach to assessing situations where response
values to an assessment differ as a function of an external covariate;
for example, gender or treatment condition. In many contexts, the main
goal of DIF analysis is to evaluate whether the items on an assessment
are biased in regards to these external covariates. Traditional DIF
methods require analysts pre-specify a set of anchor items (items
assumed to not have DIF). In the robust DIF method, no such requirement
is made. (See technical notes for more details)

The following example demonstrates how to use `robustDIF` to investigate
DIF and potential impact across treatment condition on a set of
mathematical achievement items.

## Effect of an Adaptive Game-Based Math Learning App on Students’ Learning

In a recent block randomized control trial, Bang, Li, and Flynn (2022)
investigated whether using My Math Academy, a personalized adaptive
learning app, could improve math achievement and learning engagement in
K-1 students (n=505 treatment, n=481 control). My Math Academy
supplemented school curricula with game-based activities, performance
dashboards for teachers, and offline additional learning activities.
Between February and June 2019, 41 classrooms from 11 elementary schools
participated. Treatment teachers were tasked with using My Math Academy
for at least 60 minutes per week. Math achievement was measured using
items selected from the Certica assessment item bank at both pre-test
and post-test for both conditions. Findings indicated students using My
Math Academy had greater math achievement gains than students who did
not.

All data were collected from the Item Response Warehouse
(<https://itemresponsewarehouse.org/>).

## Calculating multiple group 2PL model

The following code utilizes `mirt` to build 2PL IRT models for future
testing of `robustDIF`:

``` r
# Subset data to just items
items <- data[,-39]

# Calculate 1-factor 2PL models, using treatment to split groups and specifying SE=TRUE for the covariance matrix.
mirt.dat <- multipleGroup(items,
                      model = 1,
                      group=data$treat,
                      itemtype = "2PL",
                      SE=TRUE)
# Plot the IRFs
plot(mirt.dat, type = "trace", facet = F)
```

![IRF](figures/mirt_plot.png)

IRF

A useful first step is to investigate the IRFs of the multiple group 2PL
model using [`plot()`](https://rdrr.io/r/graphics/plot.default.html).

## The Robust DIF procedure

The
[`get_model_parms()`](https://peterhalpin.github.io/robustDIF/reference/get_model_parms.md)
function from `robustDIF` can now be used to extract the estimates from
the `mirt` object. After, robust DIF can be investigated using the
[`rdif()`](https://peterhalpin.github.io/robustDIF/reference/rdif.md)
function. Users supply a significance level by setting `alpha` (here,
`.05`) and testing for DIF on slope (discrimination), intercept
(difficulty), or both with `fun`. Here, we choose `d_fun3` to test for
DIF on both.

``` r
# Save model parameters
parms <- get_model_parms(mirt.dat)

# Investigate DIF
mod <- rdif(mle = parms, fun = "d_fun3", alpha = .05)
# Print estimate
print(mod)
```

    ## Est: 0.1226934     SE: 0.07434124

``` r
# Print summary
summary(mod)
```

    ## Robust Scaling and Differential Item Functioning.
    ## 
    ## Data: 38 items 
    ## Estimation ended after  24  iterations
    ## Single solution found
    ## 
    ## Est: 0.123    SE: 0.0743 
    ## 
    ## Results from Wald Tests of DIF:
    ##                 delta         se      z.test        p.val
    ## item1_d1   0.09285407 0.14521912  0.63940664 0.5225584296
    ## item2_d1  -0.09048042 0.12665285 -0.71439704 0.4749816995
    ## item3_d1  -0.20016615 0.11147067 -1.79568451 0.0725447062
    ## item4_d1  -0.01502962 0.15825763 -0.09496933 0.9243391874
    ## item5_d1   0.03525616 0.15834030  0.22266070 0.8237995932
    ## item6_d1  -0.04327801 0.18176764 -0.23809522 0.8118072386
    ## item7_d1  -0.10431595 0.10434728 -0.99969973 0.3174558425
    ## item8_d1   0.08815631 0.10604524  0.83130849 0.4057993780
    ## item9_d1  -0.03237231 0.13040411 -0.24824610 0.8039440015
    ## item10_d1 -0.31758571 0.09942323 -3.19428090 0.0014017965
    ## item11_d1 -0.08136337 0.19205515 -0.42364588 0.6718240927
    ## item12_d1  0.11811708 0.20383831  0.57946460 0.5622757258
    ## item13_d1 -0.11100245 0.07504851 -1.47907592 0.1391200265
    ## item14_d1 -0.04689359 0.06946288 -0.67508851 0.4996195344
    ## item15_d1  0.05314942 0.08707391  0.61039430 0.5416006462
    ## item16_d1 -0.18112391 0.11786067 -1.53676299 0.1243513585
    ## item17_d1 -0.22699593 0.13921062 -1.63059346 0.1029761310
    ## item18_d1  0.30608446 0.47382700  0.64598357 0.5182899966
    ## item19_d1  0.02570377 0.12260780  0.20964218 0.8339469553
    ## item20_d1 -0.03059457 0.11151030 -0.27436539 0.7838038410
    ## item21_d1  0.14305968 0.15488568  0.92364695 0.3556701595
    ## item22_d1  0.09457635 0.16317216  0.57961081 0.5621770997
    ## item23_d1  0.12233002 0.10121235  1.20864717 0.2267984274
    ## item24_d1 -0.10386980 0.26185639 -0.39666706 0.6916129940
    ## item25_d1  0.41146895 0.16919289  2.43195174 0.0150177081
    ## item26_d1  1.53923012 0.89267914  1.72428150 0.0846570371
    ## item27_d1  0.92867227 0.33442653  2.77690966 0.0054878430
    ## item28_d1  0.17229595 0.12093534  1.42469476 0.1542454552
    ## item29_d1  0.33277965 0.13656601  2.43676774 0.0148191945
    ## item30_d1  0.86082824 0.37496578  2.29575146 0.0216900964
    ## item31_d1  0.63261467 0.17002927  3.72062225 0.0001987325
    ## item32_d1 -0.33722912 0.14013518 -2.40645589 0.0161081491
    ## item33_d1  0.02591741 0.17939249  0.14447322 0.8851268003
    ## item34_d1  0.17575687 0.21062048  0.83447188 0.4040151275
    ## item35_d1  0.02202072 0.10209143  0.21569604 0.8292246868
    ## item36_d1 -0.18297001 0.16305709 -1.12212241 0.2618103975
    ## item37_d1  0.69401393 0.35140583  1.97496419 0.0482722122
    ## item38_d1 -0.36777386 0.15407353 -2.38700227 0.0169863866

The [`print()`](https://rdrr.io/r/base/print.html) function provides the
scaling parameter (`0.12`) and standard error (`0.07`) estimated using
iteratively reweighted least squares with Tukey’s bisquare, and
[`summary()`](https://rdrr.io/r/base/summary.html) provides additional
information regarding Wald tests on each of the items. Significant
p-values indicate that, at the chosen `alpha`, the item was flagged for
DIF. Those items are downweighted to zero during estimation of the
scaling parameter. `delta` is the estimated scaling parameter subtracted
from the item-level scaling function value.

The items that indicate DIF on both intercepts and slopes are: item 10,
item 25, item 27, item 29, item 30, item 31, item 32, item 37, and item
38 - almost a quarter (24%) of test items. This is a nontrivial amount
of item-level bias, and can be considered a form of measurement
noninvariance.

## The Rho Function

It is useful to use the
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) function to
visually inspect the Rho Function for a clear global minimum before
proceeding with analyses and making inferences about DIF.

``` r
# Plot Rho Function
plot(mod)
```

![](Testing%20for%20impact_files/figure-html/unnamed-chunk-4-1.png)

Here, the Rho Function has a clear global minimum.

## Testing for impact

Researchers may wish to know if any DIF present affects conclusions
about impact, or how the groups differ on the latent trait. In
`robustDIF`, the function
[`delta_test()`](https://peterhalpin.github.io/robustDIF/reference/delta_test.md)
accomplishes this by comparing a naive estimate of impact (here, the
unweighted mean of the item-level scaling functions) to the robust
estimate above.

``` r
delta_test(mod)
```

    ##   naive.est   naive.se  rdif.est    rdif.se     delta   delta.se   z.test
    ## 1 0.2385313 0.07744853 0.1226934 0.07434124 0.1158379 0.04237391 2.733709
    ##         p.val
    ## 1 0.006262537

In this output, `delta = naive.est - rdif.est`. The null hypothesis
(that both estimators are consistent for the true scaling parameter) is
rejected at `p=0.006`. There is evidence that DIF affects conclusions
about impact: treatment condition impacts the distribution of the latent
trait. Comparisons of total test scores between the two groups may not
necessarily be accurate comparisons of mathematical achievement ability.
