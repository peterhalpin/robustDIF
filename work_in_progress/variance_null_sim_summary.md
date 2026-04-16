# Variance Estimator Comparison Under Null DIF (2PL)

## Design
- Groups: 2
- `N` per group: 500
- Items: 15 (2PL, binary)
- Null setup: no DIF on any item (same item parameters in both groups)
- Item parameters re-drawn each replication:
  - `a_j ~ lognormal(0, 0.25^2)`
  - `d_j ~ N(0, 1)`
- Latent trait: `theta_g ~ N(0, 1)` in each group
- Main run: `R = 1000` replications

## Methods compared for `Var(theta)`
- `precision`: inverse-variance weights from `diag(vcov_y(...))`
- `irls_fixed`: IRLS bisquare weights at `theta_hat`, treated as fixed
- `irls_theta0`: IRLS bisquare weights at true `theta = 0`, treated as fixed
- `psi_prime`: derivative weights (unclamped)
- `psi_prime_clamped`: derivative weights with clamping `pmax(psi_prime, 0)`

## Main Results (`R = 1000`)
Source: `work_in_progress/results_variance_null_r1000_core_v2/summary_methods.csv`

- Empirical `Var(theta)`: `0.007368`
- Mean estimated `Var(theta)`:
  - `precision`: `0.005578` (ratio `0.757`)
  - `irls_fixed`: `0.005885` (ratio `0.799`)
  - `irls_theta0`: `0.005980` (ratio `0.812`)
  - `psi_prime`: `0.012187` (ratio `1.654`)
  - `psi_prime_clamped`: `0.006644` (ratio `0.902`)

- Empirical mean `Var(y - theta)`: `0.031364`
- Mean estimated `Var(y - theta)`:
  - `precision`: ratio `0.957`
  - `irls_fixed`: ratio `0.968`
  - `irls_theta0`: ratio `0.969`
  - `psi_prime`: ratio `1.171`
  - `psi_prime_clamped`: ratio `0.993`

- Empirical `Var(ybar - theta)`: `0.002181`
- Mean estimated `Var(ybar - theta)`:
  - `precision`: ratio `0.216`
  - `irls_fixed`: ratio `0.372`
  - `irls_theta0`: ratio `0.393`
  - `psi_prime`: ratio `3.292`
  - `psi_prime_clamped`: ratio `0.735`

Source: `work_in_progress/results_variance_null_r1000_core_v2/summary_other.csv`
- Mean item-level `Var(y_i)`:
  - empirical: `0.035243`
  - model-based: `0.035670`
- `Var(ybar)`:
  - empirical: `0.006060`
  - model-based: `0.006132`

## Bootstrap Comparator (parametric, smaller benchmark)
- Dedicated run: `R = 100`, bootstrap on first 20 reps, `B = 20` each
- Source: `work_in_progress/results_variance_null_r100_boot/summary_bootstrap.csv`
- Bootstrap mean `Var(theta)` ratio: `1.549` (overestimation in this smaller benchmark)

## Recommendation
For practical use in this framework, `psi_prime_clamped` is the best overall compromise across targets.

Why:
- It is the closest to empirical `Var(theta)` among stable estimators in the `R=1000` run.
- It is also closest for `Var(y - theta)` (nearly unbiased).
- It still underestimates `Var(ybar - theta)`, but far less severely than precision/IRLS variants and avoids the instability/overestimation seen with raw `psi_prime`.

Additional notes:
- Precision and IRLS-fixed variants systematically underestimate `Var(theta)` and strongly underestimate `Var(ybar - theta)`.
- Unclamped `psi_prime` appears unstable and overestimates all key targets.
- Treating IRLS weights at true `theta=0` improves little over IRLS-at-`theta_hat`.
