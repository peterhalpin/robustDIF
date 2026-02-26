
# -------------------------------------------------------------------
#' S3 print method for objects of class "rdif"
#' Prints the estimated scaling parameter from rdif. (to be added: SE)
#'
#' @param x An object of class 'rdif', the saved output from the \code{rdif()} function.
#' @param ... Additional arguments passed through from the generic.
#'
#' @return The input object, returned invisibly.
#'
#' @examples
#' \dontrun{
#' mod <- rdif(rdif.eg)
#' print(mod)
#' }
print.rdif <- function(x, ...) {
  if (!inherits(x, "rdif")) stop("Object is not of class 'rdif'")
  est <- x$est
  se <- x[["delta.test"]][["rdif.se"]]
  if (is.null(est)) {
    cat("Estimated scaling parameter is null.\n")
  } else {
    cat("Est:", est, "    SE:", se, "\n")
  }
  invisible(x)
}

# -------------------------------------------------------------------
#' S3 summary method for objects of class "rdif"
#'
#' @param object An object of class 'rdif', a saved list of values from \code{rdif()}
#'
#' @return A printed summary of values
#'
#' @examples
#' \dontrun{
#' mod <- rdif(rdif.eg)
#' summary(mod)
#' }

summary.rdif <- function(object, ...) {

  df.name <- "unknown"
  if (!is.null(object[["data"]][["est"]][["group.1"]])) {
    df.name <- paste0(nrow(object[["data"]][["est"]][["group.1"]]), " items")
  }
  n.iter <- as.character(object$n.iter)
  n.sols <- object$multiple.solutions
  theta <- as.character(round(object[["est"]], 3))
  se <- as.character(round(object[["delta.test"]][["rdif.se"]], 4))
  walds <- as.data.frame(object[["dif.test"]])

 cat("Robust Scaling and Differential Item Functioning.\n\n")
 cat("Data:", df.name, "\n")
 cat("Estimation ended after ", n.iter, " iterations\n")
 if (n.sols) {
   cat("Multiple solutions found\n\n")
 }
 if (!n.sols) {
   cat("Single solution found\n\n")
 }
 cat("Est:", theta, "   SE:", se, "\n\n")

 cat("Results from Wald Tests of DIF:\n")
 print(walds)

 invisible(NULL)

}

# -------------------------------------------------------------------
#' S3 plot method for objects of class "rdif"
#' Plots the rho function of the output from rdif.
#'
#' @param x A saved output from \code{rdif()} that is to be plotted.
#' @param ... Additional arguments to be passed to \code{plot()}
#'
#' @return The input object, returned invisibly. Called for plotting side effects.
#'
#' @examples
#' \dontrun{
#' # Assuming "rdif.eg" is a list of parameter values
#' rho <- rdif(mle = rdif.eg, fun = "d_fun3")
#' plot(rho)
#' }
plot.rdif <- function(x, ...) {
  if (!inherits(x, "rdif")) stop("Object is not of class 'rdif'")

  # Validate output and plot
  if (!is.null(x[["rho.plot"]][["theta"]]) && !is.null(x[["rho.plot"]][["rho"]])) {
    plot(x[["rho.plot"]][["theta"]], x[["rho.plot"]][["rho"]], type = "l", xlab = "theta", ylab = "Rho", ...)
  } else {
    stop("rho_grid did not return 'theta' and 'rho' components suitable for plotting.")
  }

  invisible(x)
}

# -------------------------------------------------------------------
#' Registers S3 methods at load time:
#' - print for class "rdif"
#' - plot for class "rdif"
#' - summary for class "rdif"
#'
#' @param libname Character string with the path to the package library.
#' @param pkgname Character string with the package name.
#'
#' @return No return value, called for side effects when the package loads.
#'
.onLoad <- function(libname, pkgname) {
  envir <- parent.env(environment())
  ns <- asNamespace(pkgname)

  # Register print method for class "rdif" if available.
  if (exists("print.rdif", mode = "function", envir = envir)) {
    registerS3method("print", "rdif", get("print.rdif", envir = envir), envir = ns)
  }

  # Register plot method for class "rdif" if available.
  if (exists("plot.rdif", mode = "function", envir = envir)) {
    registerS3method("plot", "rdif", get("plot.rdif", envir = envir), envir = ns)
  }

  # Register summary method for class "rdif" if available.
  if (exists("summary.rdif", mode = "function", envir = envir)) {
    registerS3method("summary", "rdif", get("summary.rdif", envir = envir), envir = ns)
  }
}
