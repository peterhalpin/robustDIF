# Internal input validation helpers for public robustDIF functions.

check_mle <- function(mle) {
  if (!is.list(mle)) {
    stop("`mle` must be a list returned by `get_model_parms()`.", call. = FALSE)
  }
  if (!all(c("est", "var.cov") %in% names(mle))) {
    stop("`mle` must contain `est` and `var.cov` components.", call. = FALSE)
  }
  if (!is.list(mle$est) || !is.list(mle$var.cov) || length(mle$est) < 2) {
    stop("`mle$est` and `mle$var.cov` must be lists with one entry per group (>= 2 groups).", call. = FALSE)
  }
  if (length(mle$est) != length(mle$var.cov)) {
    stop("`mle$est` and `mle$var.cov` must have the same number of groups.", call. = FALSE)
  }
  invisible(TRUE)
}

check_fun <- function(fun) {
  allowed <- c("a_fun1", "a_fun2", "d_fun1", "d_fun2", "d_fun3")
  if (!is.character(fun) || length(fun) != 1 || !fun %in% allowed) {
    stop("`fun` must be one of: 'a_fun1', 'a_fun2', 'd_fun1', 'd_fun2', 'd_fun3'.", call. = FALSE)
  }
  invisible(TRUE)
}

check_alpha <- function(alpha) {
  if (!is.numeric(alpha) || length(alpha) != 1 || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a single numeric value in (0, 1).", call. = FALSE)
  }
  invisible(TRUE)
}

check_theta <- function(theta, allow_null = FALSE) {
  if (allow_null && is.null(theta)) {
    return(invisible(TRUE))
  }
  if (!is.numeric(theta) || length(theta) != 1 || !is.finite(theta)) {
    stop("`theta` must be a single finite numeric value.", call. = FALSE)
  }
  invisible(TRUE)
}

check_k <- function(k) {
  if (!is.numeric(k) || length(k) != 1 || !is.finite(k) || k <= 0) {
    stop("`k` must be a single positive finite numeric value.", call. = FALSE)
  }
  invisible(TRUE)
}

check_grid_width <- function(grid.width) {
  if (!is.numeric(grid.width) || length(grid.width) != 1 || !is.finite(grid.width) || grid.width <= 0) {
    stop("`grid.width` must be a single positive finite numeric value.", call. = FALSE)
  }
  invisible(TRUE)
}

check_starting_value <- function(starting.value) {
  allowed <- c("med", "lts", "min_rho", "all")
  if (is.character(starting.value) && length(starting.value) == 1 && starting.value %in% allowed) {
    return(invisible(TRUE))
  }
  if (is.numeric(starting.value) && length(starting.value) >= 1 && all(is.finite(starting.value))) {
    return(invisible(TRUE))
  }
  stop("`starting.value` must be one of 'med', 'lts', 'min_rho', 'all', or a finite numeric value.", call. = FALSE)
}

check_method <- function(method) {
  if (!is.character(method) || length(method) != 1 || !method %in% c("irls", "newton")) {
    stop("`method` must be either 'irls' or 'newton'.", call. = FALSE)
  }
  if (method != "irls") {
    stop("Only `method = 'irls'` is currently implemented.", call. = FALSE)
  }
  invisible(TRUE)
}
