
# -------------------------------------------------------------------
#' S3 print method for objects of class "rdif"
#' Prints the estimated scaling parameter from rdif. (to be added: SE)
#'
#' @param object An object of class 'rdif', the saved output from the \code{rdif()} function.
#'
#' @return An object of class rdif, a string with rdif estimate value.
#' @export
#'
#' @examples
#' \dontrun {
#' mod <- rdif(rdif.eg)
#' print(mod)
#' }
print.rdif <- function(object) {
  if (!inherits(object, "rdif")) stop("Object is not of class 'rdif'")
  est <- object$est
  se <- object[["delta.test"]][["rdif.se"]]
  if (is.null(est)) {
    cat("Estimated scaling parameter is null.\n")
  } else {
    cat("Est:", est, "    SE:", se, "\n")
  }
  invisible(object)
}

# -------------------------------------------------------------------
#' S3 summary method for objects of class "rdif"
#'
#' @param object An object of class 'rdif', a saved list of values from \code{rdif()}
#'
#' @return A printed summary of values
#' @export
#'
#' @examples
#' \dontrun{
#' mod <- rdif(rdif.eg)
#' summary(mod)
#' }

summary.rdif <- function(object) {

  df.name <- as.character(object[["df"]])
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
#' @param object A saved output from \code{rdif()} that is to be plotted.
#' @param ... Additional arguments to be passed to \code{plot()}
#'
#' @return An object of class rdif, a plot with components:
#' \itemize {
#' \item{\code{xlab}}{theta}
#' \item{\code{ylab}}{rho}
#' \item{\code{type}}{line}
#' }
#' @export
#'
#' @examples
#' \dontrun{
#' # Assuming "rdif.eg" is a list of parameter values
#' rho <- rdif(mle = rdif.eg, fun = "d_fun3", grid.width = .01)
#' plot(rho)
#' }
plot.rdif <- function(object, ...) {
  if (!inherits(object, "rdif")) stop("Object is not of class 'rdif'")

  # Validate output and plot
  if (!is.null(object[["rho.plot"]][["theta"]]) && !is.null(object[["rho.plot"]][["rho"]])) {
    plot(object[["rho.plot"]][["theta"]], object[["rho.plot"]][["rho"]], type = "l", xlab = "theta", ylab = "Rho", ...)
  } else {
    stop("rho_grid did not return 'theta' and 'rho' components suitable for plotting.")
  }

  invisible(object)
}

# -------------------------------------------------------------------
#' Registers S3 methods at load time:
#' - print for class "rdif"
#' - plot for class "rdif"
#' - summary for class "rdif"
#'
#' @param libname
#' @param pkgname
#'
#' @return
#' @export
#'
.onLoad <- function(libname, pkgname) {
  envir <- parent.env(environment())
  ns <- asNamespace(pkgname)

  # Register print method for class "rdif" if available.
  if (exists("print.rdif", mode = "function", envir = envir)) {
    registerS3method("print", "rdif", get("print.rdif", envir = envir), envir = ns)
  }

  # Register plot method for class "rdif" if available.
  if (exists("plot.rdif.rho", mode = "function", envir = envir)) {
    registerS3method("plot", "rdif", get("plot.rdif.rho", envir = envir), envir = ns)
  }

  # Register summary method for class "rdif" if available.
  if (exists("summary.rdif", mode = "function", envir = envir)) {
    registerS3method("summary", "rdif", get("summary.rdif", envir = envir), envir = ns)
  }
}
