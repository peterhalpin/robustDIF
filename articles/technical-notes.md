# robustDIF Technical Notes

These notes outline an updated version of the robust DIF procedure
described in Halpin (2022) that is now implemented by the `robustDIF`
package.

Assume two groups of respondents with sample sizes $n_{1}$ and $n_{2}$
and let $n = n_{1} + n_{2}$. Also let
$Y = \left( Y_{1},\ldots,Y_{m} \right)^{\top}$ denote item-level
statistics derived from the parameter estimates of items
$i = 1,\ldots m$. The asymptotic arguments presented below assume that
$n_{1}$ and $n_{2}$ go to infinity at the same rate but that the number
of items $m$ remains fixed and finite.

The $Y_{i}$ are chosen so that
$\sqrt{n}\left( \mathbf{Y} - \theta_{0}\mathbf{1} \right)\overset{d}{\rightarrow}N\left( 0,\Sigma_{0} \right)$
under the null hypothesis of no DIF on any item. In this set-up, DIF
means that $Y_{i}$ converges in probability to a fixed value other than
$\theta_{0}$. Some specific choices of $Y_{i}$ are detailed in Halpin
(2024, 2025) and Halpin & Gilbert (2025). The notation
$V_{0}(\mathbf{Y}) = \Sigma_{0}/n$ is used to denote the finite sample
variance of $\mathbf{Y}$, and similarly for other statistics.

The robust DIF procedure can be seen as solving two interrelated
problems. First, it provides an M-estimator of $\theta_{0}$ that is
highly robust to DIF. Second, it provides a procedure for flagging items
with DIF, which happens automatically as a by-product of estimating
$\theta_{0}$. Standard Wald tests of DIF are also available.

Halpin (2022) used the (unstated) assumption that $\Sigma_{0}$ is
diagonal and attempted to combine efficiency and robustness in a way
that did not clearly separate these two opposing considerations. These
notes address these shortcomings and implement a new, simpler and more
general, version of the robust DIF procedure. The analytical details are
stated so that they apply to any bounded, redescending loss function,
although for computation the focus is on Tukey’s bisquare. Differences
with Halpin (2022) are pointed out as they arise.

## The Updated Robust DIF Procedure

### Defining the Estimator

The robust estimator $\widetilde{\theta}$ can be defined in three
related ways. Let
$u_{i} = u_{i}(\theta) = \left( Y_{i} - \theta \right)/s_{i}$ where
$s_{i} > 0$ are item-specific scaling factors to be chosen subsequently.
The three definitions of $\widetilde{\theta}$ are as follows.

1.  The minimizing argument of the loss function:
    $$R(\theta) = \sum\limits_{i = 1}^{m}\rho\left( u_{i}(\theta) \right).$$
    This definition is useful for deriving results about the robustness
    of $\widetilde{\theta}$. For redescending loss functions there
    exists constants $c$ and $k$ such that
    $\rho\left( u_{i} \right) = c$ whenever
    $\left| u_{i} \right| \geq k$. It is usual to scale $\rho$ so that
    $c = 1$. The constant $k$ is treated as a tuning parameter that
    serves to identify outliers (i.e., items with DIF).

2.  The solution to the estimating equation:
    $$\Psi(\theta) = \sum\limits_{i = 1}^{m}\psi\left( u_{i}(\theta) \right)/s_{i} = 0,$$
    where $\psi(u) = \rho\prime(u)$. The influence function
    $\Psi(\theta)$ is important for obtaining the variance of
    $\widetilde{\theta}$.

3.  A weighted mean that is obtained by defining the weights
    $w_{i}^{*}(\theta) = \psi\left( u_{i}(\theta) \right)/u_{i}(\theta)$
    and substituting these into the estimating equation to get:
    $$\theta = \sum\limits_{i = 1}^{m}w_{i}(\theta)\, Y_{i}\qquad\text{with}\qquad w_{i}(\theta) = \frac{w_{i}^{*}(\theta)/s_{i}^{2}}{\sum\limits_{j = 1}^{m}w_{j}^{*}(\theta)/s_{j}^{2}}.$$
    By convention, $w_{i}^{*}(0) \equiv 1$ to avoid division by zero.
    Also note that, when $\left| u_{i} \right| \geq k$,
    $\psi\left( u_{i} \right) = w_{i}^{*}(\theta) = 0$, so that outliers
    (as defined by $k$) are “redescended” to zero. The weighted mean is
    useful for computation via iteratively re-weighted least squares
    (IRLS). Especially for redescending loss functions, IRLS can be much
    more stable than Newton-based methods that solve $\Psi(\theta) = 0$
    (because $\Psi\prime(\theta)$ can approach zero, leading the Newton
    steps to diverge).

### Choosing the Scaling Factors $s_{i}$

The scaling factors $s_{i}$ are required to ensure that
$\widetilde{\theta}$ is equivariant under re-scaling of the $Y_{i}$. In
conventional applications, the $Y_{i}$ are “raw” data points and the
scaling factors $s_{i} = s$ are chosen to be an ancillary estimate of
the scale of the $Y_{i}$ (e.g., the median absolute deviation or MAD).
In this situation, the scaling factors are constant over $i$, so they
factor out of the estimating equation and cancel out in the numerator
and denominator of the normalized weights $w_{i}$.

In the present application, item-specific scaling factors $s_{i}$ are
available because we can derive $\Sigma_{0}$ (the asymptotic covariance
matrix of the $Y_{i}$ under the null hypothesis of no DIF) by applying
the Delta method to the item parameter estimates. As shown above, this
somewhat complicates the relationship among the different definitions of
$\widetilde{\theta}$, because it is now important to keep track of how
the item-specific scaling factors appear in the different definitions.
However, the item-specific scaling factors are worth this additional
complication for the following three reasons, all of which were noted in
Halpin (2022).

- First, obtaining item-specific scaling factors $s_{i}$ analytically
  from $\Sigma_{0}$ means that we no longer require an ancillary
  estimate based on the scale of the realized values of $Y_{i}$. This is
  important because it leads the resulting estimator to be highly robust
  to DIF. This is also the main detail that separates the proposed
  approach from that considered (and dismissed) by Stocking and Lord
  (1984). Stocking and Lord used $s = MAD\left( Y_{i} \right)$. However,
  the MAD has a breakdown point of 1/4, so any estimator of $\theta$
  that uses $s = MAD\left( Y_{i} \right)$ will breakdown if $\geq 1/4$
  of the items exhibit DIF (see Huber and Roncetti, 2009, chapter 6). By
  contrast $\Sigma_{0}$ does not depend directly on the scale of the
  realized values of $Y_{i}$ – it can be computed directly from the item
  parameters. The overall result is that the robustness of the resulting
  M-estimator no longer depends on that of an ancillary scaling factor.

- The value of $Y_{i}$ appears in the expression for
  $V_{0}\left( Y_{i} \right)$. Thus,
  ${\widehat{V}}_{0}\left( Y_{i} \right)$ may yet be contaminated by
  DIF, leading to the potential for “masking”. Although this problem is
  not as severe as when using an ancillary estimate like the MAD, it is
  still a potential concern. This problem can be avoided as follows. The
  null hypothesis that the item does not exhibit DIF gives
  $Y_{i}\overset{p}{\rightarrow}\theta_{0}$. This motivates using the
  substitution $Y_{i} = \theta^{\star}$ when estimating $\Sigma_{0}$,
  where $\theta^{\star}$ is a consistent, high-breakdown estimate of
  $\theta$ (e.g., the median). The overall result is a plug-in estimator
  of $V_{0}\left( Y_{i} \right)$ that is robust to DIF.

- Third, using item-specific scaling factors implies that we can
  downweight items with DIF at the desired asymptotic false positive
  rate during IRLS-based estimation. For example, if we choose
  $s_{i} = \sqrt{V_{0}\left( Y_{i} \right)}$ and $k$ as the
  $1 - \alpha/2$ quantile of $N(0,1)$, then items are down-weighted to
  zero once $\left| Y_{i} \right|$ lies beyond $(1 - \alpha) \times 100$
  confidence interval (CI) centered at $\theta$. In this way, DIF
  detection arises as a by-product of robust scaling.

Halpin (2022) chose $s_{i} = V_{0}\left( Y_{i} \right)$ based on a
(flawed) argument about efficiency in the absence of DIF, which also
complicated the choice of the tuning parameter $k$. The argument was
flawed because (a) it was based on the unstated assumption that
$\Sigma_{0}$ is diagonal, which is not true for many IRT estimators, and
(b) it did not account for the scaling factors that appear outside the
influence function in the estimating equation (see point 2 in the
previous section).

To address these issues, the robust DIF estimator has been updated to
use $s_{i} = \sqrt{V_{0}\left( Y_{i} \right)}$ with tuning parameter $k$
based on the asymptotic CI rationale outlined above. This approach
ignores sampling variation in $\widetilde{\theta}$ when downweighting
items with DIF. A more accurate downweighting procedure could instead be
based on
$s_{i} = \sqrt{V_{0}\left( Y_{i} - \widetilde{\theta} \right)}$. Since
$Y_{i}$ and $\widetilde{\theta}$ are positively correlated,
$V_{0}\left( Y_{i} - \widetilde{\theta} \right) \leq V_{0}\left( Y_{i} \right)$.
Thus, the flagging procedure based on $V_{0}\left( Y_{i} \right)$ is
somewhat anti-conservative. To address this issue, it is recommended to
compute item-by-item tests of DIF using a standard Wald test following
estimation of $\widetilde{\theta}$:

$$z_{i} = \frac{Y_{i} - \widetilde{\theta}}{\sqrt{V_{0}\left( Y_{i} - \widetilde{\theta} \right)}}.$$

Note that modifying the estimator $\widetilde{\theta}$ to instead use
$s_{i} = \sqrt{V_{0}\left( Y_{i} - \widetilde{\theta} \right)}$ is
possible (see Halpin), but this leads to complications obtaining its
asymptotic distribution. In practice, these complications do not appear
to be worth the trouble as there is little change in finite sample
performance when using the simpler approach outlined above.

### The Asymptotic Distribution of $\widetilde{\theta}$

Halpin (2022) obtained the asymptotic distribution of
$\widetilde{\theta}$ using the Delta method and the implicit function
theorem. The derivation is recapitulated here for general (i.e.,
non-diagonal) $\Sigma_{0}$, and the results presented in Halpin (2022)
are seen to follow from the assumption that $\Sigma_{0}$ is diagonal.

The estimator $\widetilde{\theta} = g(\mathbf{Y})$ is implicitly defined
as the solution to the estimating equation

$$\Psi(\theta;\mathbf{Y}) = \sum\limits_{i = 1}^{m}\frac{\psi\left( {(Y_{i} - \theta})/s_{i} \right)}{s_{i}} = 0.$$

Let the asymptotic distribution of $\mathbf{Y}$ be denoted as

$$\sqrt{n}(\mathbf{Y} - {\mathbf{μ}})\overset{p}{\rightarrow}N(0,\Sigma)$$
with the null hypothesis of no DIF leading to
${\mathbf{μ}} = \theta_{0}\mathbf{1}$ and $\Sigma = \Sigma_{0}$. Also
let $\theta^{\star}$ be defined as any solution to the population
estimating equation
$E_{\mathbf{Y}}\left\lbrack \Psi(\theta;\mathbf{Y}) \right\rbrack = 0$.
There may be multiple local solutions when using a redescending loss
function, and the asymptotic results described here apply to any local
solution. In practice, local minima can be diagnosed by plotting
$R(\theta)$ over a grid of $\theta$ values.

The following assumptions are used:

- A1: $\psi(u)$ continuously differentiable.
- A2: $\psi(u)$ is odd (i.e. $\psi(-u) = -\psi(u)$).
- A3: $\psi\prime(0) \neq 0$.
- A4: $\Psi\prime\left( \theta^{\star};{\mathbf{μ}} \right) \neq 0$.

A1 allows the Delta method to be applied to $g(\mathbf{Y})$. A2 ensures
that $\theta_{0}$ is a solution to
$E_{\mathbf{Y}}\left\lbrack \Psi(\theta;\mathbf{Y}) \right\rbrack$ under
the null hypothesis. A3 and continuity (A1) imply that the population
estimating equation is monotone around $\theta_{0}$ and hence that
$\theta_{0}$ is a locally unique solution. A4 is required by the
implicit function theorem. Under the null hypothesis, A3 implies A4, but
in general the two assumptions are distinct.

Applying the Delta method gives (using A1)

$$\sqrt{n}\left( \widetilde{\theta} - \theta^{\star} \right)\overset{p}{\rightarrow}N\left( 0,\ \nabla g({\mathbf{μ}})^{\top}\Sigma\;\nabla g({\mathbf{μ}}) \right).$$
The gradient of $g(\mathbf{Y})$ is obtained from the implicit function
theorem (using A4):

$$\nabla g(\mathbf{Y}) = -\left( \frac{\partial\Psi(\theta;\mathbf{Y})}{\partial\theta} \right)^{-1}\frac{\partial\Psi(\theta;\mathbf{Y})}{\partial\mathbf{Y}}.$$

Evaluating the partial derivatives gives

$$\frac{\partial\Psi(\theta;\mathbf{Y})}{\partial\theta} = -\sum\limits_{i = 1}^{m}\frac{\psi\prime\left( \left( Y_{i} - \theta \right)/s_{i} \right)}{s_{i}^{2}}$$

and
$$\frac{\partial\Psi(\theta;\mathbf{Y})}{\partial Y_{i}} = \frac{\psi\prime\left( \left( Y_{i} - \theta \right)/s_{i} \right)}{s_{i}^{2}}.$$

Therefore the gradient $\nabla g(\mathbf{Y})$ has elements

$$\frac{\partial g(\mathbf{Y})}{\partial Y_{i}} = \frac{\psi\prime\left( \left( Y_{i} - \theta \right)/s_{i} \right)/s_{i}^{2}}{\sum\limits_{j = 1}^{m}\psi\prime\left( \left( Y_{j} - \theta \right)/s_{j} \right)/s_{j}^{2}}.$$

The foregoing results provide the general (i.e., non-null) distribution
of $\widetilde{\theta}$.

Next we consider the null distribution. First we show that
$\theta^{\star} = \theta_{0}$ and then derive
$V_{0}\left( \widetilde{\theta} \right)$. The null hypothesis is that
${\mathbf{μ}} = \theta_{0}\mathbf{1}$. This implies that the
standardized residuals
$U_{0i} = \left( Y_{i} - \theta_{0} \right)/s_{i}$ are symmetrically
distributed about zero. Combined with the assumption that $\psi(u)$ is
odd (A2), this gives $\theta^{\star} = \theta_{0}$ by the following
argument:

$$E_{\mathbf{Y}}\left\lbrack \Psi\left( \theta_{0};\mathbf{Y} \right) \right\rbrack = \sum\limits_{i}E_{\mathbf{Y}}\left\lbrack \psi\left( U_{0i} \right) \right\rbrack = \sum\limits_{i}E_{\mathbf{Y}}\left\lbrack \psi\left( -U_{0i} \right) \right\rbrack = \sum\limits_{i}E_{\mathbf{Y}}\left\lbrack -\psi\left( U_{0i} \right) \right\rbrack = -E_{\mathbf{Y}}\left\lbrack \Psi\left( \theta_{0};\mathbf{Y} \right) \right\rbrack.$$
The second equality follows from the symmetry of $U_{0i}$ about zero and
the third from A2. The chain of equalities shows that
$E_{\mathbf{Y}}\left\lbrack \Psi\left( \theta_{0};\mathbf{Y} \right) \right\rbrack = 0$.
Together A1 and A3 ensure that this is a locally unique solution.

To obtain the null variance, we evaluate the gradient at
$\mathbf{Y} = {\mathbf{μ}} = \theta_{0}\mathbf{1}$:

$$\left. \frac{\partial g(\mathbf{Y})}{\partial Y_{i}} \right|_{\mathbf{Y} = \theta_{0}\mathbf{1}} = \frac{\psi\prime\left( \left( \theta_{0} - \theta_{0} \right)/s_{i} \right)/s_{i}^{2}}{\sum\limits_{j = 1}^{m}\psi\prime\left( \left( \theta_{0} - \theta_{0} \right)/s_{j} \right)/s_{j}^{2}} = \frac{1/s_{i}^{2}}{\sum\limits_{j = 1}^{m}1/s_{j}^{2}}.$$

The second equality follows from A3, which implies that $\psi\prime(0)$
is equal to a non-zero constant that factors out of the numerator and
denominator.

Finally, using $s_{i} = \sqrt{V_{0}\left( Y_{i} \right)}$ the gradient
$\nabla g({\mathbf{μ}})$ becomes a vector of precision weights. Letting
$\mathbf{p} = \left( p_{1},\ldots,p_{n} \right)^{\top}$ denote the
vector of precision weights, we can write the asymptotic null
distribution of $\widetilde{\theta}$ as

$$\sqrt{n}\left( \widetilde{\theta} - \theta_{0} \right)\overset{p}{\rightarrow}N\left( 0,\mathbf{p}^{\top}\Sigma_{0}\;\mathbf{p} \right).$$

Under the additional assumption that $\Sigma_{0}$ is diagonal, the
resulting expression for the null variance of $\widetilde{\theta}$ is

$$V_{0}\left( \widetilde{\theta} \right) = \sum\limits_{i = 1}^{m}p_{i}^{2}\, V_{0}\left( Y_{i} \right) = \sum\limits_{i = 1}^{m}\left( \frac{1/V_{0}\left( Y_{i} \right)}{\sum\limits_{j = 1}^{m}1/V_{0}\left( Y_{i} \right)} \right)^{2}V_{0}\left( Y_{i} \right) = \left( \sum\limits_{j = 1}^{m}1/V_{0}\left( Y_{i} \right) \right)^{-1}.$$

This is the result given in part (a) of Theorem 1 in Halpin (2022).

A similar argument gives the asymptotic null distribution of
$Y_{i} - \widetilde{\theta}$ as

$$\sqrt{n}\left( Y_{i} - \widetilde{\theta} \right)\overset{p}{\rightarrow}N\left( 0,\ \left( \mathbf{e}_{i} - \mathbf{p} \right)^{\top}\Sigma_{0}\;\left( \mathbf{e}_{i} - \mathbf{p} \right) \right)$$
where $\mathbf{e}_{i}$ is the $i$-th column of the identity matrix.
Under the additional assumption that $\Sigma_{0}$ is diagonal, the
resulting expression for the variance is

$$V_{0}\left( Y_{i} - \widetilde{\theta} \right) = V_{0}\left( Y_{i} \right) + V_{0}\left( \widetilde{\theta} \right) - 2p_{i}V_{0}\left( Y_{i} \right) = V_{0}\left( Y_{i} \right) + V_{0}\left( \widetilde{\theta} \right) - 2V_{0}\left( \widetilde{\theta} \right) = V_{0}\left( Y_{i} \right) - V_{0}\left( \widetilde{\theta} \right)$$
This is the result given in part (b) of Theorem 1 in Halpin (2022).

Halpin (2022) used the same overall approach to compare two esimates of
$\theta$ – the unweighted mean of $Y_{i}$ and the robust estimator
outlined above.

## Implementation via IRLS

The robust DIF estimator $\widetilde{\theta}$ can be computed using IRLS
based on the weighted mean definition given above. The IRLS algorithm is
as follows:

1.  Initialize $\theta^{(0)}$ and set iteration counter $t = 0$.  
2.  Compute standardized residuals
    $u_{i}^{(t)} = \left( Y_{i} - \theta^{(t)} \right)/s_{i}$.
3.  Compute weights
    $w_{i}\left( \theta^{(t)} \right) = \psi\left( u_{i}^{(t)} \right)/u_{i}^{(t)}$
    with $w_{i}^{*}(0) \equiv 1$.
4.  Update
    $\theta^{(t + 1)} = \sum_{i = 1}^{m}w_{i}\left( \theta^{(t)} \right)\, Y_{i}$
    with

$$w_{i}\left( \theta^{(t)} \right) = \frac{w_{i}^{*}\left( \theta^{(t)} \right)/s_{i}^{2}}{\sum\limits_{j = 1}^{m}w_{j}^{*}\left( \theta^{(t)} \right)/s_{j}^{2}}$$.

5.  If $\left| \theta^{(t + 1)} - \theta^{(t)} \right| < \epsilon$ stop;
    else set $t = t + 1$ and return to step 2.

Once the algorithm has converged, set
$\widetilde{\theta} = \theta^{(t + 1)}$. Item-level DIF tests can be
conducted using the Wald test statistic given above or the
multi-parameter variant given in Halpin (2022).

## References

Halpin, P.F. (2022) Differential Item Functioning Via Robust Scaling.
Arxiv Preprint. <https://arxiv.org/abs/2207.04598>. Published in
Psychometrika in 2024 under the same title.

Halpin, P.F. (2024) Differential Test Functioning Via Robust Scaling.
Arxiv Preprint. <https://arxiv.org/abs/2409.03502>.

Halpin, P. F., & Gilbert, J. (2025). Testing Whether Reported Treatment
Effects Are Unduly Influenced by Item-Level Heterogeneity. PsyArxiv
Preprint. <https://doi.org/10.31234/osf.io/9ru45_v1>
