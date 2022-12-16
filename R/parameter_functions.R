# -------------------------------------------------------------------
#' Extract and format item parameter estimates and their covariance matrix
#'
#' @description
#' Takes a list of 1-factor model fits from \code{\link[mirt]{mirt}} (with \code{SE = TRUE}), \code{\link[lavaan]{cfa}} (with \code{std.lv = TRUE}), or \code{\link[MplusAutomation]{readModels}} (with TECH3) and formats the item parameter estimates and their covariance matrix. All \code{robustDIF} functions assume that the estimates were obtained by maximum likelihood and the covariance is asymptotically correct.
#'
#' Note that the only type of fit currently supported is the \code{SingleGroupClass} of the \code{\link[mirt]{mirt}} package.
#'
#' It is possible to use fits from other software with \code{robustDIF} functions, but the parameter estimates and their covariance matrices must be formatted in a particular way. For more details, see the documentation for the example dataset \code{\link[robustDIF]{rdif.eg}}.
#'
#'
#' @param object a list of model fits (current implementation only supports objects from \code{\link[mirt]{mirt}} with \code{itemtype} of any combination of \code{"2PL", "graded", or "gpcm"} and \code{SE = TRUE})
#' @return A list of item parameter estimates and covariance matrices for each group.
#'
#' @seealso \code{\link[robustDIF]{rdif.eg}}
#' @export
#'
# -------------------------------------------------------------------

get_params <- function(object) {

  if(inherits(object, "list")){
    if(inherits(object[[1]], "SingleGroupClass")){
      temp <- lapply(object, get_mirt_params)

      out <- list(param.names = temp[[1]]$param.names,
                  param.est = lapply(temp, "[[", 2),
                  var.cov = lapply(temp, "[[", 3))

    } else if(inherits(object[[1]], "mplus.model")){
      out <- lapply(object, get_mplus_params)

    } else if(inherits(object[[1]], "lavaan")){
      out <- lapply(object, get_lavaan_params)
    }

  } else if(inherits(object, "MultipleGroupClass")){
    out <- get_mirt_params(object)

  }
  return(out)
}

#-------------------------------------------------------------------
#' Extract item parameter estimates and their covariance matrix from \code{\link[mirt]{mirt}}.
#'
#' @param mirt.object a \code{\link[mirt]{mirt}} object (\code{SingleGroupClass}). Expected to be a 1-factor model with \code{itemtype} of any combination of \code{"2PL", "graded", or "gpcm"}.
#' @return A three-element \code{list}:
#' \itemize{
#' \item vector of parameter names taking the form "item.param"
#' \item vector of item parameter estimates
#' \item covariance matrix of item parameter estimates
#' }
#'
#' @importFrom mirt coef
#' @importFrom mirt vcov
#' @export
# -------------------------------------------------------------------

get_mirt_params <- function(mirt.object){

  # check for 1-factor model
  if(mirt.object@Model$nfact != 1){
    stop("mirt.object must be a 1-factor model.")
  }
  # need to check itemtype? mirt.object@Model$itemtype

  ## number of items
  n.items <- mirt.object@Data$nitems

  if(inherits(mirt.object, "MultipleGroupClass")) mo <- mirt.object@ParObjects$pars[[1]] else mo <- mirt.object

  ## extract slope and intercept parameter names for each item
  item.params <- lapply(1:n.items, function(x){
    parnames <- mo@ParObjects$pars[[x]]@parnames
    parnames[grepl("^[ad][1-9]|[ad]$", parnames)]
  })
  original.names <- colnames(mirt.object@Data$data) # original item names
  internal.names <- paste0("item", 1:n.items) # item names to be used for robustDIF functions

  ## combining item and parameter names into single vector - one for each type of item name
  param.names <- list(internal = unlist(lapply(1:n.items, function(x) paste(internal.names[[x]], item.params[[x]], sep = "."))),
                      original = unlist(lapply(1:n.items, function(x) paste(original.names[[x]], item.params[[x]], sep = "."))))

  ## extracting parameter estimates
  param.est <- mirt.object@Internals$shortpars # call made by mirt::extract.mirt(object, "parvec")

  if(inherits(mirt.object, "MultipleGroupClass")){
    param.est <- split(param.est, cut(seq_along(param.est), mirt.object@Data$ngroups, labels = FALSE))
  }


  ## extracting vcov matrix and removing dimnames
  v <- mirt::vcov(mirt.object)
  if(inherits(mirt.object, "MultipleGroupClass")){

    gs <- paste0("g", 1:mirt.object@Data$ngroups)
    ipg <- paste(rep(param.names[[1]], times = length(gs)), # create all item.prameter.group names
                 rep(gs, each = length(param.names[[1]])),
                 sep = ".")

    row.names(v) <- colnames(v) <- ipg

    ## subset vcov by group
    v <- lapply(gs, function(x) v[grepl(x, row.names(v)), grepl(x, colnames(v))])
    v <- lapply(v, function(x){
      row.names(x) <- colnames(x) <- NULL
      return(x)
    })

  } else{
    row.names(v) <- colnames(v) <- NULL
  }


  return(list(param.names = param.names,
              param.est = param.est,
              var.cov = v))
}

#-------------------------------------------------------------------
#' Helper function used within \code{get_mirt_params} to format parameters
#'
#' @param par.mat matrix of item parameter estimates extracted from a 1-factor \code{\link[mirt]{mirt}} object
#' @param item.names character vector of item names
#' @return a two-element \code{list} containing a vector of parameter names and vector of parameter estimates
#' @export # temporarily for testing
# -------------------------------------------------------------------

format_mirt_params <- function(par.mat, item.names){

  ## initial item-parameter combinations
  par.names <- colnames(par.mat)[grepl("^[ad][1-9]|[ad]$", colnames(par.mat))] # all possible intercept or slope parameters
  item.pars <- paste(rep(item.names, each = length(par.names)), # all item by parameter combinations
                     rep(par.names, times = length(item.names)),
                     sep = ".")

  ## creating vector of parameter estimates and final item-parameter combinations
  par.mat <- par.mat[,colnames(par.mat) %in% par.names] # remove non-slope/intercept params
  par.vec <- c(t(par.mat))       # convert from matrix to vector arranged by item then param (slope, intercepts); contains NAs
  names(par.vec) <- item.pars
  par.vec <- par.vec[!is.na(par.vec)] # remove NA parameters, and therefore, irrelevant item.pars labels
  item.pars <- names(par.vec) # grab correct item.pars
  # names(par.vec) <- NULL # remove to be consistent with other classes?

  return(list(param.names = item.pars,
              parameters = par.vec))
}

#-------------------------------------------------------------------
#' Extract item parameter estimates and their covariance matrix from \code{lavaan}.
#'
#' @param lavaan.object an object of class \code{lavaan}. Expected to be a 1-factor model estimated with \code{std.lv = TRUE}.
#' @return A three-element \code{list}:
#' \itemize{
#' \item vector of parameter names taking the form "item.param"
#' \item vector of item parameter estimates
#' \item covariance matrix of item parameter estimates
#' }
#'
#' @importFrom lavaan lavInspect
#' @export
# -------------------------------------------------------------------

get_lavaan_params <- function(lavaan.object){

  if(!lavaan::lavInspect(lavaan.object, what = "options")$std.lv){
    stop("Models must be run with std.lv = TRUE to obtain complete item parameter vcov matrix for DIF procedure.")
  }

  ## extracting parameters
  pars <- lavaan::lavInspect(lavaan.object, what = "est")

  slopes <- data.frame(item = row.names(pars$lambda),
                       param = "a",
                       est = pars$lambda[,1])

  ints <- data.frame(item = gsub("\\|.*$", "", row.names(pars$tau)),
                     param = gsub("^.*\\|", "", row.names(pars$tau)),
                     est = pars$tau[,1])

  ## binding and sorting parameters
  par.df <- rbind(slopes, ints)
  ov.names <- lavaan::lavNames(lavaan.object, "ov")
  par.df$item <- factor(par.df$item, levels = ov.names) # make sure items in correct order
  par.df <- par.df[with(par.df, order(item, param)), ]  # sorting parameter estimates

  ## extracting and re-organizing vcov matrix
  v <- lavaan::lavInspect(lavaan.object, what = "vcov")
  v2 <- v[grepl("=~|\\|", row.names(v)), grepl("=~|\\|", colnames(v))] # ensure matrix only has intercept/threshold and slope vcovs
  colnames(v2) <- gsub("f1=~", "", colnames(v2))  # remove factor prefix from item names
  row.names(v2) <- gsub("f1=~", "", row.names(v2))
  v3 <- v2[row.names(par.df), row.names(par.df)] # re-arrange by item, then parameter rather than the reverse
  # row.names(v3) <- colnames(v3) <- paste(par.df$item, par.df$param, sep = ".") # do we need to bother renaming?

  return(list(param.names = paste(par.df$item, par.df$param, sep = "."),
       parameters = par.df$est,
       var.cov = v3))
}

#-------------------------------------------------------------------
#' Extract item parameter estimates and their covariance matrix from MPlus.
#'
#' @param mplus.object an object of class \code{mplus.model} from \code{\link[MplusAutomation]{readModels}}. Expected to be a 1-factor model. Mplus .out file must include Tech3.
#' @return A three-element \code{list}:
#' \itemize{
#' \item vector of parameter names taking the form "item.param"
#' \item vector of item parameter estimates
#' \item covariance matrix of item parameter estimates
#' }
#'
#' @export
# -------------------------------------------------------------------

get_mplus_params <- function(mplus.object){

  if(length(mplus2pl.mods$twopl_groupa.out$tech3) == 0){
    stop("mplus.object must include tech3, which contains the covariance matrix of item parameter estimates")
  }


  ## extracting parameters
  pars <- mplus.object$parameters$unstandardized # same as standardized

  # removing unnecessary information
  slopes.all <- pars[grepl("BY$", pars$paramHeader), 2:3]
  ints.all <- pars[pars$paramHeader %in% c("Thresholds","Steps"), 2:3]

  # formatting estimates
  slopes <- data.frame(item = slopes.all$param,
                       param = "a",
                       est = slopes.all$est)

  ints <- data.frame(item = gsub("\\$.*$", "", ints.all$param),
                     param = paste0("d",gsub("^.*\\$", "", ints.all$param)),
                     est = ints.all$est)

  ## binding and sorting parameters
  n.items <- nrow(slopes.all)
  # item.names <- paste0("i", 1:n.items)
  item.names <- slopes.all$param # user name (after Mplus capitalizes it)

  par.df <- rbind(slopes, ints)
  par.df$item <- factor(par.df$item, levels = item.names) # make sure items in correct order
  par.df <- par.df[with(par.df, order(item, param)), ]    # sorting parameter estimates

  ## extracting and re-organizing vcov matrix
  v <- mplus.object$tech3$paramCov      # Mplus order is the same as lavaan (see Tech 1)
  v[upper.tri(v)] <- t(v)[upper.tri(v)] # fill in upper triangle
  v2 <- v[row.names(par.df), row.names(par.df)] # re-arrange by item, then parameter rather than the reverse
  # row.names(v2) <- colnames(v2) <- paste(par.df$item, par.df$param, sep = ".") # do we need to bother renaming?

  return(list(param.names = paste(par.df$item, par.df$param, sep = "."),
              parameters = par.df$est,
              var.cov = v2))
}
