suppressPackageStartupMessages({
  library(mirt)
  library(Matrix)
  library(pkgload)
})

pkgload::load_all("/Users/peterhalpin/git/robustDIF", export_all = FALSE, quiet = TRUE)

inv_logit <- function(x) 1 / (1 + exp(-x))

simulate_group_data <- function(n, a, d, theta_mean = 0, theta_sd = 1) {
  theta <- rnorm(n, mean = theta_mean, sd = theta_sd)
  p <- sapply(seq_along(a), function(j) inv_logit(a[j] * theta + d[j]))
  y <- matrix(rbinom(n * length(a), size = 1, prob = as.vector(p)), nrow = n, ncol = length(a))
  colnames(y) <- paste0("item", seq_along(a))
  y
}

normalize_weights <- function(w_num) {
  s <- sum(w_num)
  if (!is.finite(s) || abs(s) < 1e-12) {
    return(rep(NA_real_, length(w_num)))
  }
  w_num / s
}

safe_scalar_var <- function(w, vc) {
  if (any(!is.finite(w))) {
    return(NA_real_)
  }
  as.numeric(t(w) %*% vc %*% w)
}

compute_method_weights <- function(y, vcov_y_theta, theta_hat, k) {
  var_y <- diag(vcov_y_theta)
  m <- length(y)

  w_prec <- normalize_weights(1 / var_y)

  u_hat <- (y - theta_hat) / sqrt(var_y)
  w_irls_hat <- normalize_weights(robustDIF:::bsq_weight(u_hat, k) / var_y)
  w_psi_hat <- normalize_weights(robustDIF:::psi_prime(u_hat, k) / var_y)
  w_psi_clamp_hat <- normalize_weights(pmax(robustDIF:::psi_prime(u_hat, k), 0) / var_y)

  u0 <- y / sqrt(var_y)
  w_irls_theta0 <- normalize_weights(robustDIF:::bsq_weight(u0, k) / var_y)

  list(
    precision = w_prec,
    irls_fixed = w_irls_hat,
    irls_theta0 = w_irls_theta0,
    psi_prime = w_psi_hat,
    psi_prime_clamped = w_psi_clamp_hat
  )
}

estimate_components <- function(mle, fun = "d_fun3", alpha = 0.05) {
  mod <- rdif(mle = mle, fun = fun, alpha = alpha)
  theta_hat <- mod$est
  k <- mod$k

  y <- y_fun(mle, fun = fun)
  m <- length(y)
  vbar <- rep(1 / m, m)

  # Under null, requested true generating value.
  vcov_y_theta0 <- vcov_y(mle, theta = 0, fun = fun)
  var_y_est <- diag(vcov_y_theta0)
  ybar <- mean(y)
  ybar_var_hat <- as.numeric(t(vbar) %*% vcov_y_theta0 %*% vbar)

  w_methods <- compute_method_weights(
    y = y,
    vcov_y_theta = vcov_y_theta0,
    theta_hat = theta_hat,
    k = k
  )

  method_out <- lapply(w_methods, function(w) {
    theta_var <- safe_scalar_var(w, vcov_y_theta0)
    A <- diag(m) - matrix(1, nrow = m, ncol = 1) %*% t(w)
    y_minus_theta_var <- diag(A %*% vcov_y_theta0 %*% t(A))
    ybar_minus_theta_var <- safe_scalar_var(vbar - w, vcov_y_theta0)

    list(
      theta_var_hat = theta_var,
      y_minus_theta_var_hat = y_minus_theta_var,
      ybar_minus_theta_var_hat = ybar_minus_theta_var
    )
  })

  list(
    mod = mod,
    y = y,
    ybar = ybar,
    ybar_var_hat = ybar_var_hat,
    var_y_hat = var_y_est,
    methods = method_out
  )
}

parametric_bootstrap_theta_var <- function(mle, B = 200, fun = "d_fun3", alpha = 0.05, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  g1 <- mle$est$group.1
  g2 <- mle$est$group.2
  # Null: groups have equal parameters; use pooled item parameters.
  a <- (g1$a1 + g2$a1) / 2
  d <- (g1[, -1, drop = FALSE] + g2[, -1, drop = FALSE]) / 2
  d <- as.numeric(d[, 1])

  n1 <- 500
  n2 <- 500
  thetas <- numeric(B)
  for (b in seq_len(B)) {
    y1 <- simulate_group_data(n = n1, a = a, d = d, theta_mean = 0, theta_sd = 1)
    y2 <- simulate_group_data(n = n2, a = a, d = d, theta_mean = 0, theta_sd = 1)
    dat <- rbind(y1, y2)
    grp <- rep(c("g1", "g2"), c(n1, n2))
    fit <- multipleGroup(
      data = dat,
      model = 1,
      group = grp,
      itemtype = "2PL",
      SE = TRUE,
      verbose = FALSE,
      technical = list(NCYCLES = 500)
    )
    mle_b <- get_model_parms(fit)
    thetas[b] <- rdif(mle = mle_b, fun = fun, alpha = alpha)$est
  }
  var(thetas)
}

run_sim <- function(
  R = 1000,
  n_per_group = 500,
  n_items = 15,
  fun = "d_fun3",
  alpha = 0.05,
  seed = 20260226,
  bootstrap_reps = 200,
  bootstrap_B = 200,
  checkpoint_every = 25,
  out_dir = "/Users/peterhalpin/git/robustDIF/work_in_progress/results_variance_null"
) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  set.seed(seed)

  methods <- c("precision", "irls_fixed", "irls_theta0", "psi_prime", "psi_prime_clamped")

  theta <- numeric(R)
  ybar <- numeric(R)
  ybar_var_hat <- numeric(R)
  y_vec <- matrix(NA_real_, nrow = R, ncol = n_items)
  y_minus_theta <- matrix(NA_real_, nrow = R, ncol = n_items)
  ybar_minus_theta <- numeric(R)

  theta_var_hat <- matrix(NA_real_, nrow = R, ncol = length(methods), dimnames = list(NULL, methods))
  ybar_minus_theta_var_hat <- matrix(NA_real_, nrow = R, ncol = length(methods), dimnames = list(NULL, methods))
  y_minus_theta_var_hat <- array(NA_real_, dim = c(R, n_items, length(methods)), dimnames = list(NULL, NULL, methods))
  y_var_hat <- matrix(NA_real_, nrow = R, ncol = n_items)

  pb_theta_var <- rep(NA_real_, R)

  for (r in seq_len(R)) {
    # Random 2PL parameters each replication.
    a <- exp(rnorm(n_items, mean = 0, sd = 0.25))
    d <- rnorm(n_items, mean = 0, sd = 1)

    y1 <- simulate_group_data(n = n_per_group, a = a, d = d)
    y2 <- simulate_group_data(n = n_per_group, a = a, d = d)
    dat <- rbind(y1, y2)
    grp <- rep(c("g1", "g2"), c(n_per_group, n_per_group))

    fit <- multipleGroup(
      data = dat,
      model = 1,
      group = grp,
      itemtype = "2PL",
      SE = TRUE,
      verbose = FALSE,
      technical = list(NCYCLES = 500)
    )

    mle <- get_model_parms(fit)
    comp <- estimate_components(mle, fun = fun, alpha = alpha)

    theta[r] <- comp$mod$est
    ybar[r] <- comp$ybar
    ybar_var_hat[r] <- comp$ybar_var_hat
    y_vec[r, ] <- comp$y
    y_minus_theta[r, ] <- comp$y - comp$mod$est
    ybar_minus_theta[r] <- comp$ybar - comp$mod$est
    y_var_hat[r, ] <- comp$var_y_hat

    for (m in methods) {
      theta_var_hat[r, m] <- comp$methods[[m]]$theta_var_hat
      ybar_minus_theta_var_hat[r, m] <- comp$methods[[m]]$ybar_minus_theta_var_hat
      y_minus_theta_var_hat[r, , m] <- comp$methods[[m]]$y_minus_theta_var_hat
    }

    if (r <= bootstrap_reps) {
      pb_theta_var[r] <- tryCatch(
        parametric_bootstrap_theta_var(mle, B = bootstrap_B, fun = fun, alpha = alpha, seed = seed + r),
        error = function(e) NA_real_
      )
    }

    if (r %% checkpoint_every == 0 || r == R) {
      saveRDS(
        list(
          meta = list(
            R = R,
            n_per_group = n_per_group,
            n_items = n_items,
            fun = fun,
            alpha = alpha,
            seed = seed,
            bootstrap_reps = bootstrap_reps,
            bootstrap_B = bootstrap_B,
            completed = r
          ),
          theta = theta[1:r],
          ybar = ybar[1:r],
          ybar_var_hat = ybar_var_hat[1:r],
          y_vec = y_vec[1:r, , drop = FALSE],
          y_minus_theta = y_minus_theta[1:r, , drop = FALSE],
          ybar_minus_theta = ybar_minus_theta[1:r],
          theta_var_hat = theta_var_hat[1:r, , drop = FALSE],
          ybar_minus_theta_var_hat = ybar_minus_theta_var_hat[1:r, , drop = FALSE],
          y_minus_theta_var_hat = y_minus_theta_var_hat[1:r, , , drop = FALSE],
          y_var_hat = y_var_hat[1:r, , drop = FALSE],
          pb_theta_var = pb_theta_var[1:r]
        ),
        file = file.path(out_dir, "sim_checkpoint.rds")
      )
      cat(sprintf("Completed %d/%d reps\n", r, R))
    }
  }

  # Empirical ("truth" in this simulation) variances.
  theta_var_emp <- var(theta)
  y_var_emp <- apply(y_vec, 2, var)
  ybar_var_emp <- var(ybar)
  y_minus_theta_var_emp <- apply(y_minus_theta, 2, var)
  ybar_minus_theta_var_emp <- var(ybar_minus_theta)

  summarize_method <- function(method) {
    data.frame(
      method = method,
      theta_var_emp = theta_var_emp,
      theta_var_hat_mean = mean(theta_var_hat[, method], na.rm = TRUE),
      theta_ratio = mean(theta_var_hat[, method], na.rm = TRUE) / theta_var_emp,
      y_minus_theta_var_emp_mean = mean(y_minus_theta_var_emp),
      y_minus_theta_var_hat_mean = mean(y_minus_theta_var_hat[, , method], na.rm = TRUE),
      y_minus_theta_ratio = mean(y_minus_theta_var_hat[, , method], na.rm = TRUE) / mean(y_minus_theta_var_emp),
      ybar_minus_theta_var_emp = ybar_minus_theta_var_emp,
      ybar_minus_theta_var_hat_mean = mean(ybar_minus_theta_var_hat[, method], na.rm = TRUE),
      ybar_minus_theta_ratio = mean(ybar_minus_theta_var_hat[, method], na.rm = TRUE) / ybar_minus_theta_var_emp
    )
  }

  summary_methods <- do.call(rbind, lapply(methods, summarize_method))

  summary_other <- data.frame(
    quantity = c("mean_var_y", "var_ybar"),
    empirical = c(mean(y_var_emp), ybar_var_emp),
    model_based_mean = c(mean(y_var_hat, na.rm = TRUE), mean(ybar_var_hat, na.rm = TRUE))
  )

  summary_boot <- data.frame(
    n_boot_reps = sum(is.finite(pb_theta_var)),
    theta_var_emp = theta_var_emp,
    pb_theta_var_mean = mean(pb_theta_var, na.rm = TRUE),
    pb_theta_ratio = mean(pb_theta_var, na.rm = TRUE) / theta_var_emp
  )

  write.csv(summary_methods, file.path(out_dir, "summary_methods.csv"), row.names = FALSE)
  write.csv(summary_other, file.path(out_dir, "summary_other.csv"), row.names = FALSE)
  write.csv(summary_boot, file.path(out_dir, "summary_bootstrap.csv"), row.names = FALSE)

  saveRDS(
    list(
      summary_methods = summary_methods,
      summary_other = summary_other,
      summary_boot = summary_boot,
      theta = theta,
      ybar = ybar,
      ybar_var_hat = ybar_var_hat,
      y_vec = y_vec,
      y_minus_theta = y_minus_theta,
      ybar_minus_theta = ybar_minus_theta,
      theta_var_hat = theta_var_hat,
      ybar_minus_theta_var_hat = ybar_minus_theta_var_hat,
      y_minus_theta_var_hat = y_minus_theta_var_hat,
      y_var_hat = y_var_hat,
      pb_theta_var = pb_theta_var
    ),
    file = file.path(out_dir, "sim_results_full.rds")
  )

  list(
    summary_methods = summary_methods,
    summary_other = summary_other,
    summary_boot = summary_boot,
    out_dir = out_dir
  )
}

# Intentionally no auto-run block.
# Call run_sim(...) explicitly from an interactive session or wrapper script.
