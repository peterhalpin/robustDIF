#' Simulate RDIF power for a 2PL DIF design.
#'
#' Runs a small simulation study for the RDIF DTF Wald test. The function is
#' intended for planning and sensitivity analyses like those reported in the
#' Halpin2024 archive branch. It uses two independently fitted 2PL models,
#' extracts item parameters with \code{\link[robustDIF]{get_model_parms}}, and
#' applies \code{\link[robustDIF]{rdif}} with optional BH-start calibration.
#'
#' @param n.reps number of Monte Carlo replications.
#' @param n.persons number of persons per group, or length-two vector for
#'   reference and comparison group sample sizes.
#' @param n.items number of binary 2PL items.
#' @param n.dif number of items with DIF.
#' @param item.delta target item-level effect size for DIF items on the
#'   \code{d_fun3} scale.
#' @param slope.multiplier multiplicative slope DIF applied to DIF items in the
#'   comparison group.
#' @param alpha nominal Type I error rate.
#' @param alpha.adjust passed to \code{\link[robustDIF]{rdif}}.
#' @param fun item-level scaling function. Currently only \code{"d_fun3"} is
#'   supported because \code{item.delta} is imposed on that scale.
#' @param reference.mean,reference.var latent mean and variance in the reference
#'   group.
#' @param comparison.mean.range,comparison.var.range length-two ranges from
#'   which the comparison-group latent mean and variance are sampled in each
#'   replication.
#' @param a.range,b.range ranges for uniformly sampled reference-group slopes
#'   and difficulties.
#' @param seed optional random seed.
#' @param n.cores number of cores for Unix-like parallel execution. On Windows,
#'   the function uses sequential execution.
#' @param keep.mle logical; if \code{TRUE}, store formatted MLE objects for each
#'   replication.
#'
#' @return A list with replicate-level DTF results, item-level diagnostics, a
#'   compact summary table, optionally saved MLEs, and the matched call.
#' @examples
#' \dontrun{
#' out <- rdif_power_sim(n.reps = 25, n.persons = 500, n.dif = 4,
#'                       item.delta = 0.5, alpha.adjust = "BH")
#' out$summary
#' }
#' @export
rdif_power_sim <- function(n.reps = 100,
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
                           keep.mle = FALSE) {
  alpha.adjust <- match.arg(alpha.adjust)
  rdif_power_check_args(
    n.reps, n.persons, n.items, n.dif, item.delta, slope.multiplier,
    alpha, fun, reference.var, comparison.mean.range, comparison.var.range,
    a.range, b.range)

  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (length(n.persons) == 1L) {
    n.persons <- rep(n.persons, 2)
  }

  seeds <- sample.int(.Machine$integer.max, n.reps)
  run_one <- function(rep.index) {
    set.seed(seeds[rep.index])
    rdif_power_rep(
      rep.index = rep.index,
      n.persons = n.persons,
      n.items = n.items,
      n.dif = n.dif,
      item.delta = item.delta,
      slope.multiplier = slope.multiplier,
      alpha = alpha,
      alpha.adjust = alpha.adjust,
      fun = fun,
      reference.mean = reference.mean,
      reference.var = reference.var,
      comparison.mean.range = comparison.mean.range,
      comparison.var.range = comparison.var.range,
      a.range = a.range,
      b.range = b.range,
      keep.mle = keep.mle)
  }

  if (.Platform$OS.type != "windows" && n.cores > 1) {
    reps <- parallel::mclapply(
      seq_len(n.reps),
      run_one,
      mc.cores = n.cores,
      mc.set.seed = FALSE)
  } else {
    reps <- lapply(seq_len(n.reps), run_one)
  }

  delta <- rdif_power_rbind(lapply(reps, `[[`, "delta"))
  items <- rdif_power_rbind(lapply(reps, `[[`, "items"))
  mles <- if (keep.mle) lapply(reps, `[[`, "mle") else NULL

  list(
    delta = delta,
    items = items,
    summary = rdif_power_summary(delta, items, alpha),
    mles = mles,
    call = match.call())
}

rdif_power_check_args <- function(n.reps, n.persons, n.items, n.dif,
                                  item.delta, slope.multiplier, alpha, fun,
                                  reference.var, comparison.mean.range,
                                  comparison.var.range, a.range, b.range) {
  if (!identical(fun, "d_fun3")) {
    stop("rdif_power_sim currently supports fun = \"d_fun3\" only.",
         call. = FALSE)
  }
  scalar_positive <- function(x) {
    is.numeric(x) && length(x) == 1L && is.finite(x) && x > 0
  }
  if (!scalar_positive(n.reps) || n.reps != as.integer(n.reps)) {
    stop("n.reps must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(n.persons) || !length(n.persons) %in% c(1L, 2L) ||
      any(!is.finite(n.persons)) || any(n.persons <= 0) ||
      any(n.persons != as.integer(n.persons))) {
    stop("n.persons must be a positive integer or a length-two integer vector.",
         call. = FALSE)
  }
  if (!scalar_positive(n.items) || n.items != as.integer(n.items)) {
    stop("n.items must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(n.dif) || length(n.dif) != 1L || !is.finite(n.dif) ||
      n.dif < 0 || n.dif > n.items || n.dif != as.integer(n.dif)) {
    stop("n.dif must be an integer between 0 and n.items.", call. = FALSE)
  }
  if (!is.numeric(item.delta) || length(item.delta) != 1L ||
      !is.finite(item.delta)) {
    stop("item.delta must be a single finite numeric value.", call. = FALSE)
  }
  if (!scalar_positive(slope.multiplier)) {
    stop("slope.multiplier must be a positive finite scalar.", call. = FALSE)
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("alpha must be between 0 and 1.", call. = FALSE)
  }
  if (!scalar_positive(reference.var)) {
    stop("reference.var must be positive.", call. = FALSE)
  }
  for (x in list(comparison.mean.range, comparison.var.range, a.range, b.range)) {
    if (!is.numeric(x) || length(x) != 2L || any(!is.finite(x)) ||
        x[1] > x[2]) {
      stop("range arguments must be finite length-two numeric vectors.",
           call. = FALSE)
    }
  }
  invisible(TRUE)
}

rdif_power_rep <- function(rep.index, n.persons, n.items, n.dif, item.delta,
                           slope.multiplier, alpha, alpha.adjust, fun,
                           reference.mean, reference.var,
                           comparison.mean.range, comparison.var.range,
                           a.range, b.range, keep.mle) {
  latent.mean.1 <- reference.mean
  latent.var.1 <- reference.var
  latent.mean.2 <- stats::runif(1, comparison.mean.range[1],
                                comparison.mean.range[2])
  latent.var.2 <- stats::runif(1, comparison.var.range[1],
                               comparison.var.range[2])
  latent.sd.1 <- sqrt(latent.var.1)
  latent.sd.2 <- sqrt(latent.var.2)
  pooled.sd <- sqrt((latent.var.1 + latent.var.2) / 2)
  delta.0 <- (latent.mean.2 - latent.mean.1) / pooled.sd

  a0 <- stats::runif(n.items, a.range[1], a.range[2])
  b0 <- sort(stats::runif(n.items, b.range[1], b.range[2]))
  dif.items <- if (n.dif > 0) sort(sample(seq_len(n.items), n.dif)) else integer(0)
  target.item.delta <- rep(0, n.items)
  target.item.delta[dif.items] <- item.delta

  a1 <- a0
  a1[dif.items] <- a0[dif.items] * slope.multiplier
  d0 <- -b0 * a0
  d1 <- rdif_power_d_from_target_delta(
    a0 = a0,
    d0 = d0,
    a1 = a1,
    target.item.delta = target.item.delta,
    latent.mean.1 = latent.mean.1,
    latent.mean.2 = latent.mean.2,
    latent.sd.1 = latent.sd.1,
    latent.sd.2 = latent.sd.2)

  true.effects <- rdif_power_true_effects(
    a0 = a0,
    d0 = d0,
    a1 = a1,
    d1 = d1,
    latent.mean.1 = latent.mean.1,
    latent.mean.2 = latent.mean.2,
    latent.sd.1 = latent.sd.1,
    latent.sd.2 = latent.sd.2)
  true.delta <- mean(true.effects$item.delta)

  theta0 <- stats::rnorm(n.persons[1], latent.mean.1, latent.sd.1)
  theta1 <- stats::rnorm(n.persons[2], latent.mean.2, latent.sd.2)
  dat0 <- mirt::simdata(a0, d0, n.persons[1], "2PL", Theta = matrix(theta0))
  dat1 <- mirt::simdata(a1, d1, n.persons[2], "2PL", Theta = matrix(theta1))

  out <- tryCatch({
    fit0 <- mirt::mirt(dat0, 1, SE = TRUE, verbose = FALSE)
    fit1 <- mirt::mirt(dat1, 1, SE = TRUE, verbose = FALSE)
    mle <- get_model_parms(list(fit0, fit1))
    rdif.out <- rdif(
      mle,
      fun = fun,
      alpha = alpha,
      alpha.adjust = alpha.adjust)
    oracle <- rdif_power_delta_from_dif(mle, fun = fun, dif.items = dif.items)
    y <- y_fun(mle, fun)
    p.adj <- stats::p.adjust(rdif.out$dif.test$p.val, method = "BH")

    delta <- rdif_power_delta_rows(
      rep.index, n.persons, n.items, n.dif, item.delta, slope.multiplier,
      latent.mean.1, latent.mean.2, latent.sd.1, latent.sd.2, delta.0,
      true.delta, rdif.out$delta.test, oracle, alpha, rdif.out)
    items <- rdif_power_item_rows(
      rep.index, n.persons, n.items, n.dif, item.delta, slope.multiplier,
      latent.mean.1, latent.mean.2, latent.sd.1, latent.sd.2, delta.0,
      true.delta, seq_along(y), dif.items, true.effects$item.delta,
      rdif.out$dif.test$p.val, p.adj, rdif.out$weights, rdif.out$k,
      rdif.out$alpha.i, alpha)
    list(delta = delta, items = items, mle = if (keep.mle) mle else NULL)
  }, error = function(e) {
    meta <- rdif_power_meta(
      rep.index, n.persons, n.items, n.dif, item.delta, slope.multiplier,
      latent.mean.1, latent.mean.2, latent.sd.1, latent.sd.2, delta.0,
      true.delta)
    delta <- data.frame(
      meta,
      method = "RDIF",
      naive.est = NA_real_,
      naive.se = NA_real_,
      rdif.est = NA_real_,
      rdif.se = NA_real_,
      delta = NA_real_,
      delta.se = NA_real_,
      z.test = NA_real_,
      p.val = NA_real_,
      rejected = NA,
      covered = NA,
      error = conditionMessage(e))
    items <- data.frame(
      meta,
      item = seq_len(n.items),
      true_dif = seq_len(n.items) %in% dif.items,
      true_Delta_i = true.effects$item.delta,
      p.val = NA_real_,
      p.adjusted = NA_real_,
      flagged = NA,
      rdif_weight = NA_real_,
      k = NA_real_,
      alpha.i = NA_real_,
      error = conditionMessage(e))
    list(delta = delta, items = items, mle = NULL)
  })

  out
}

rdif_power_meta <- function(rep.index, n.persons, n.items, n.dif, item.delta,
                            slope.multiplier, latent.mean.1, latent.mean.2,
                            latent.sd.1, latent.sd.2, delta.0, true.delta) {
  data.frame(
    rep = rep.index,
    n_persons_1 = n.persons[1],
    n_persons_2 = n.persons[2],
    n_items = n.items,
    n_dif = n.dif,
    p_dif = n.dif / n.items,
    target_item_delta = item.delta,
    slope_multiplier = slope.multiplier,
    latent_mean_1 = latent.mean.1,
    latent_mean_2 = latent.mean.2,
    latent_sd_1 = latent.sd.1,
    latent_sd_2 = latent.sd.2,
    delta_0 = delta.0,
    true_Delta = true.delta)
}

rdif_power_delta_rows <- function(rep.index, n.persons, n.items, n.dif,
                                  item.delta, slope.multiplier, latent.mean.1,
                                  latent.mean.2, latent.sd.1, latent.sd.2,
                                  delta.0, true.delta, rdif.test,
                                  oracle.test, alpha, rdif.out) {
  meta <- rdif_power_meta(
    rep.index, n.persons, n.items, n.dif, item.delta, slope.multiplier,
    latent.mean.1, latent.mean.2, latent.sd.1, latent.sd.2, delta.0,
    true.delta)
  rdif.row <- data.frame(
    meta,
    method = "RDIF",
    rdif.test,
    rejected = rdif.test$p.val < alpha,
    covered = abs(rdif.test$delta - true.delta) <=
      stats::qnorm(.975) * rdif.test$delta.se,
    n_flagged = sum(stats::p.adjust(rdif.out$dif.test$p.val, "BH") <= alpha,
                    na.rm = TRUE),
    mean_k = mean(rdif.out$k, na.rm = TRUE),
    min_k = min(rdif.out$k, na.rm = TRUE),
    max_k = max(rdif.out$k, na.rm = TRUE),
    multiple_solutions = rdif.out$multiple.solutions,
    error = NA_character_)
  true.row <- data.frame(
    meta,
    method = "True",
    oracle.test,
    rejected = oracle.test$p.val < alpha,
    covered = abs(oracle.test$delta - true.delta) <=
      stats::qnorm(.975) * oracle.test$delta.se,
    n_flagged = n.dif,
    mean_k = NA_real_,
    min_k = NA_real_,
    max_k = NA_real_,
    multiple_solutions = NA,
    error = NA_character_)
  rdif_power_rbind(list(rdif.row, true.row))
}

rdif_power_item_rows <- function(rep.index, n.persons, n.items, n.dif,
                                 item.delta, slope.multiplier, latent.mean.1,
                                 latent.mean.2, latent.sd.1, latent.sd.2,
                                 delta.0, true.delta, item, dif.items,
                                 true.item.delta, p.val, p.adjusted, weights,
                                 k, alpha.i, alpha) {
  meta <- rdif_power_meta(
    rep.index, n.persons, n.items, n.dif, item.delta, slope.multiplier,
    latent.mean.1, latent.mean.2, latent.sd.1, latent.sd.2, delta.0,
    true.delta)
  data.frame(
    meta,
    item = item,
    true_dif = item %in% dif.items,
    true_Delta_i = true.item.delta,
    p.val = p.val,
    p.adjusted = p.adjusted,
    flagged = p.adjusted <= alpha,
    rdif_weight = weights,
    k = k,
    alpha.i = alpha.i,
    error = NA_character_)
}

rdif_power_d_from_target_delta <- function(a0, d0, a1, target.item.delta,
                                           latent.mean.1, latent.mean.2,
                                           latent.sd.1, latent.sd.2) {
  a0.fit <- a0 * latent.sd.1
  d0.fit <- d0 + a0 * latent.mean.1
  a1.fit <- a1 * latent.sd.2
  pooled.sd <- sqrt((latent.sd.1^2 + latent.sd.2^2) / 2)
  delta.0 <- (latent.mean.2 - latent.mean.1) / pooled.sd
  target.y <- delta.0 + target.item.delta
  denom <- sqrt((a0.fit^2 + a1.fit^2) / 2)
  d1.fit <- d0.fit + target.y * denom
  d1.fit - a1 * latent.mean.2
}

rdif_power_true_effects <- function(a0, d0, a1, d1, latent.mean.1,
                                    latent.mean.2, latent.sd.1,
                                    latent.sd.2) {
  a0.fit <- a0 * latent.sd.1
  d0.fit <- d0 + a0 * latent.mean.1
  a1.fit <- a1 * latent.sd.2
  d1.fit <- d1 + a1 * latent.mean.2
  pooled.sd <- sqrt((latent.sd.1^2 + latent.sd.2^2) / 2)
  delta.0 <- (latent.mean.2 - latent.mean.1) / pooled.sd
  item.y <- (d1.fit - d0.fit) / sqrt((a0.fit^2 + a1.fit^2) / 2)
  list(y = item.y, item.delta = item.y - delta.0, delta.0 = delta.0)
}

rdif_power_delta_from_dif <- function(mle, fun = "d_fun3", dif.items) {
  y <- y_fun(mle, fun)
  vcov.y <- vcov_y(mle, fun = fun)
  n <- length(y)
  n.anchor <- n - length(dif.items)
  if (n.anchor <= 0) {
    out <- data.frame(
      naive.est = mean(y),
      naive.se = sqrt(t(rep(1 / n, n)) %*% vcov.y %*% rep(1 / n, n))[1, 1],
      rdif.est = NA_real_,
      rdif.se = NA_real_,
      delta = NA_real_,
      delta.se = NA_real_,
      z.test = NA_real_,
      p.val = NA_real_)
    return(out)
  }
  w.bar <- rep(1 / n, n)
  w.anchor <- rep(1 / n.anchor, n)
  w.anchor[dif.items] <- 0
  delta <- sum((w.bar - w.anchor) * y)
  se.delta <- sqrt(t(w.bar - w.anchor) %*% vcov.y %*% (w.bar - w.anchor))[1, 1]
  z <- delta / se.delta
  data.frame(
    naive.est = sum(w.bar * y),
    naive.se = sqrt(t(w.bar) %*% vcov.y %*% w.bar)[1, 1],
    rdif.est = sum(w.anchor * y),
    rdif.se = sqrt(t(w.anchor) %*% vcov.y %*% w.anchor)[1, 1],
    delta = delta,
    delta.se = se.delta,
    z.test = z,
    p.val = 2 * stats::pnorm(-abs(z)))
}

rdif_power_summary <- function(delta, items, alpha = 0.05) {
  if (!nrow(delta)) {
    return(data.frame())
  }
  delta.split <- split(delta, delta$method)
  delta.summary <- rdif_power_rbind(lapply(delta.split, function(x) {
    data.frame(
      method = x$method[1],
      n_reps = nrow(x),
      rejection_rate = mean(x$rejected, na.rm = TRUE),
      bias = mean(x$delta - x$true_Delta, na.rm = TRUE),
      rmse = sqrt(mean((x$delta - x$true_Delta)^2, na.rm = TRUE)),
      coverage = mean(x$covered, na.rm = TRUE),
      mean_delta_se = mean(x$delta.se, na.rm = TRUE))
  }))

  item.summary <- data.frame(
    method = "RDIF",
    item_tpr = if (any(items$true_dif, na.rm = TRUE)) {
      mean(items$flagged[items$true_dif], na.rm = TRUE)
    } else {
      NA_real_
    },
    item_fpr = if (any(!items$true_dif, na.rm = TRUE)) {
      mean(items$flagged[!items$true_dif], na.rm = TRUE)
    } else {
      NA_real_
    },
    mean_rdif_weight_dif = if (any(items$true_dif, na.rm = TRUE)) {
      mean(items$rdif_weight[items$true_dif], na.rm = TRUE)
    } else {
      NA_real_
    },
    mean_rdif_weight_anchor = if (any(!items$true_dif, na.rm = TRUE)) {
      mean(items$rdif_weight[!items$true_dif], na.rm = TRUE)
    } else {
      NA_real_
    })

  merge(delta.summary, item.summary, by = "method", all.x = TRUE)
}

rdif_power_rbind <- function(x) {
  x <- x[!vapply(x, is.null, logical(1))]
  if (!length(x)) {
    return(data.frame())
  }
  do.call(rbind, x)
}
