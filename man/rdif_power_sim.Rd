\name{rdif_power_sim}
\alias{rdif_power_sim}
\title{Simulate RDIF power for a 2PL DIF design.}
\usage{
rdif_power_sim(
  n.reps = 100,
  n.persons = 500,
  n.items = 20,
  n.dif = 4,
  item.delta = 0.5,
  slope.multiplier = 1,
  alpha = 0.05,
  alpha.adjust = c("BH", "none"),
  fun = "d_fun3",
  reference.mean = 0,
  reference.var = 1,
  comparison.mean.range = c(-0.25, 0.25),
  comparison.var.range = c(0.5, 1.5),
  a.range = c(0.5, 2),
  b.range = c(-1.5, 1.5),
  seed = NULL,
  n.cores = 1,
  keep.mle = FALSE
)
}
\arguments{
\item{n.reps}{number of Monte Carlo replications.}

\item{n.persons}{number of persons per group, or length-two vector for
reference and comparison group sample sizes.}

\item{n.items}{number of binary 2PL items.}

\item{n.dif}{number of items with DIF.}

\item{item.delta}{target item-level effect size for DIF items on the
\code{d_fun3} scale.}

\item{slope.multiplier}{multiplicative slope DIF applied to DIF items in the
comparison group.}

\item{alpha}{nominal Type I error rate.}

\item{alpha.adjust}{passed to \code{\link[robustDIF]{rdif}}.}

\item{fun}{item-level scaling function. Currently only \code{"d_fun3"} is
supported because \code{item.delta} is imposed on that scale.}

\item{reference.mean, reference.var}{latent mean and variance in the reference
group.}

\item{comparison.mean.range, comparison.var.range}{length-two ranges from
which the comparison-group latent mean and variance are sampled in each
replication.}

\item{a.range, b.range}{ranges for uniformly sampled reference-group slopes
and difficulties.}

\item{seed}{optional random seed.}

\item{n.cores}{number of cores for Unix-like parallel execution. On Windows,
the function uses sequential execution.}

\item{keep.mle}{logical; if \code{TRUE}, store formatted MLE objects for each
replication.}
}
\value{
A list with replicate-level DTF results, item-level diagnostics, a compact
summary table, optionally saved MLEs, and the matched call.
}
\description{
Runs a small simulation study for the RDIF DTF Wald test. The function is
intended for planning and sensitivity analyses like those reported in the
Halpin2024 archive branch. It uses two independently fitted 2PL models,
extracts item parameters with \code{\link[robustDIF]{get_model_parms}}, and
applies \code{\link[robustDIF]{rdif}} with optional BH-start calibration.
}
\examples{
\dontrun{
out <- rdif_power_sim(n.reps = 25, n.persons = 500, n.dif = 4,
                      item.delta = 0.5, alpha.adjust = "BH")
out$summary
}
}
