############################################################################
# need to update Mplus
############################################################################
# -------------------------------------------------------------------
#' Extract and format item parameter estimates and their covariance matrix
#'
#' @description
#' Takes a 1-factor model fit or list of 1-factor model fits from \code{\link[mirt]{mirt}} or \code{\link[lavaan]{cfa}}
#' and formats the item parameter estimates and their covariance matrix for use in other \code{robustDIF} functions.
#'
#' @param object model fit from a multigroup analysis or list of model fits for each group for a 1-factor model. See **Details**.
#'
#' @details
#' The function takes a fitted 1-factor multigroup model or list of fitted 1-factor single group models. The factor must be standardized (i.e., variance = 1) and the covariance matrix be asymptotically correct.
#' Currently, the function accepts:
#' \itemize{
#' \item a \code{\link[mirt]{mirt}} object of class \code{SingleGroupClass} or \code{MultipleGroupClass} with \code{SE = TRUE} (to return covariance matrix) and \code{itemtype} of any combination of \code{"2PL", "graded", or "gpcm"}.
#' \item a \code{lavaan} object estimated from \code{\link[lavaan]{cfa}} with \code{std.lv = TRUE}.
#'}
#' It is possible to use fits from other software with \code{robustDIF} functions, but the parameter estimates and their covariance matrices must be formatted identically to what is returned by \code{get_model_parms}. For details, see the documentation for the example dataset \code{\link[robustDIF]{rdif.eg}}.
#'
#' @return A three-element \code{list}:
#' \itemize{
#' \item vector of parameter names taking the form "item.parameter"
#' \item list (one element per group) of vectors of item parameter estimates
#' \item list (one element per group) of covariance matrices of item parameter estimates
#' }
#'
#' @seealso \code{\link[robustDIF]{rdif.eg}}
#' @export
#'
# -------------------------------------------------------------------

get_model_parms <- function(object) {
  if(inherits(object, "list")){
    if(inherits(object[[1]], "SingleGroupClass")){
      temp <- lapply(object, get_mirt_pars)
      out <- NULL

    } else if(inherits(object[[1]], "lavaan")){
      temp <- lapply(object, get_lavaan_pars)
      out <- NULL

    # Need to recode Mplus
    # }  else if(inherits(object[[1]], "mplus.model")){ # list of single group objects
    #   temp <- lapply(object, get_mplus_params)
    #   out <- NULL
    #
    # } else if(inherits(object[[1]], "mplus.inp")){ # multigroup object
    #   out <- get_mplus_params(object)
    #   # groups in alphabetical order
    }

    if(is.null(out)){
      ## groups in order provided in list
      out <- list(par.names = temp[[1]]$par.names,
                  est = lapply(temp, "[[", 2),
                  var.cov = lapply(temp, "[[", 3))
    }

  } else if(inherits(object, "MultipleGroupClass")){
    out <- get_mirt_pars(object)
    # groups in alphabetical order;

  } else if(inherits(object, "lavaan")){
    out <- get_lavaan_pars(object)
    # groups in order observed in data; see lavaan::lavInspect(object, what = "group.label")
  }
 reformat_out(out)
}

#-------------------------------------------------------------------
#' Extract item parameter estimates and their covariance matrix from \code{\link[mirt]{mirt}}.
#'
#' @param mirt.object a \code{\link[mirt]{mirt}} object of class \code{SingleGroupClass} or \code{MultipleGroupClass}. Expected to be a 1-factor model with \code{SE = TRUE} and \code{itemtype} of any combination of \code{"2PL", "graded", or "gpcm"}.
#' @return A three-element \code{list}:
#' \itemize{
#' \item vector of parameter names taking the form "item.parameter"
#' \item list (one element per group) of vectors of item parameter estimates
#' \item list (one element per group) of covariance matrices of item parameter estimates
#' }
#'
#' @importFrom mirt coef
#' @importFrom mirt vcov
#' @export
# -------------------------------------------------------------------

get_mirt_pars <- function(mirt.object){

  ## check for 1-factor model
  if(mirt.object@Model$nfact != 1){
    stop("mirt.object must be a 1-factor model.")
  }
  # need to check itemtype? mirt.object@Model$itemtype

  ## number of items
  n.items <- mirt.object@Data$nitems

  if(inherits(mirt.object, "MultipleGroupClass")) {
    mo <- mirt.object@ParObjects$pars[[1]]
  } else {
    mo <- mirt.object
  }

  ## extract slope and intercept parameter names for each item
  item.pars <- lapply(1:n.items, function(x){
                  parnames <- mo@ParObjects$pars[[x]]@parnames
                  parnames[grepl("^[ad][1-9]|[ad]$", parnames)]})

  original.names <- colnames(mirt.object@Data$data) # original item names
  internal.names <- paste0("item", 1:n.items) # item names for robustDIF

  ## combining item and parameter names into single vector
  internal <- unlist(lapply(1:n.items,
                function(x) paste(internal.names[[x]], item.pars[[x]], sep = ".")))
  original <- unlist(lapply(1:n.items,
                function(x) paste(original.names[[x]], item.pars[[x]], sep = ".")))
  par.names <- list(internal = internal, original = original)

  ## extracting parameter estimates
  est <- mirt.object@Internals$shortpars # call made by mirt::extract.mirt(object, "parvec")

  if(inherits(mirt.object, "MultipleGroupClass")){
    est <- split(est, cut(seq_along(est), mirt.object@Data$ngroups, labels = FALSE))
  }

  ## extracting vcov matrix and removing dimnames
  v <- mirt::vcov(mirt.object)
  if(inherits(mirt.object, "MultipleGroupClass")){

    gs <- paste0("g", 1:mirt.object@Data$ngroups)
    ipg <- paste(rep(par.names[[1]], times = length(gs)), # create all item.par.group names
                 rep(gs, each = length(par.names[[1]])),
                 sep = ".")

    row.names(v) <- colnames(v) <- ipg

    ## subset vcov by group
    v <- lapply(gs, function(x) v[grepl(x, row.names(v)), grepl(x, colnames(v))])
    v <- lapply(v, function(x){row.names(x) <- colnames(x) <- NULL
                               return(x)})
  } else{
    row.names(v) <- colnames(v) <- NULL
  }
  return(list(par.names = par.names,
              est = est,
              var.cov = v))
}


#-------------------------------------------------------------------
#' Extract item parameter estimates and their covariance matrix from \code{lavaan}.
#'
#' @param lavaan.object an object of class \code{lavaan}. Expected to be a 1-factor model estimated with \code{std.lv = TRUE}.
#' @return A three-element \code{list}:
#' \itemize{
#' \item vector of parameter names taking the form "item.parameter"
#' \item list (one element per group) of vectors of item parameter estimates
#' \item list (one element per group) of covariance matrices of item parameter estimates
#' }
#'
#' @importFrom lavaan lavInspect
#' @export
# -------------------------------------------------------------------

get_lavaan_pars <- function(lavaan.object){

  if(!lavaan::lavInspect(lavaan.object, what = "options")$std.lv){
    stop("Models must be run with std.lv = TRUE to obtain complete item parameter vcov matrix for DIF procedure.")
  }

  ngroups <- lavaan::lavInspect(lavaan.object, what = "ngroups")

  ## extracting parameters and item names
  pars <- lavaan::lavInspect(lavaan.object, what = "est")    # groups in order of appearance in dataset
  original.names <- lavaan::lavNames(lavaan.object, "ov")    # original item names
  internal.names <- paste0("item", 1:length(original.names)) # item names to be used for robustDIF functions

  ## formatting parameter vector and defining full parameter names
  if(ngroups == 1){

    par.df <- format_pars(pars, names.vec = original.names, type = "lavaan")
    lav.par.names <- row.names(par.df)     # used for subsetting vcov matrix
    par.names = list(internal = "placeholder",
                     original = paste(par.df$item, par.df$par, sep = "."))
    est <- par.df$est

  } else {

    par.df <- lapply(pars, format_pars, names.vec = original.names, type = "lavaan")
    lav.par.names <- row.names(par.df[[1]])    # used for subsetting vcov matrix
    par.names = list(internal = "placeholder",
                       original = paste(par.df[[1]]$item, par.df[[1]]$par, sep = "."))
    est <- lapply(par.df, "[[", "est")
  }

  ## defining parameter names for internal use
  x <- par.names$original

  for(i in 1:length(original.names)){
    x <- sub(original.names[[i]], internal.names[[i]], x) # replacing original names with generic name
  }
  par.names$internal <- x

  ## extracting matrix
  v <- lavaan::lavInspect(lavaan.object, what = "vcov")
  v2 <- v[grepl("=~|\\|", row.names(v)), grepl("=~|\\|", colnames(v))] # ensure matrix only has intercept/threshold and slope vcovs
  colnames(v2) <- gsub(".*=~", "", colnames(v2))  # remove factor prefix from item names
  row.names(v2) <- gsub(".*=~", "", row.names(v2))

  if(ngroups == 1){
    v3 <- v2[lav.par.names, lav.par.names] # re-arrange by item
    row.names(v3) <- colnames(v3) <- NULL

  } else {

  ## define parameter name for each group using lavaan's names
  ipg <- vector("list", ngroups)
  ipg[[1]] <- lav.par.names # lavaan does not add suffix for g1
    for(i in 2:ngroups){
      ipg[[i]] <- paste0(lav.par.names, ".g", i)
    }
  ## subset vcov by group and reorder by item
  v3 <- lapply(ipg, function(x) v2[x,x])
  v3 <- lapply(v3, function(x){row.names(x) <- colnames(x) <- NULL
                     return(x)})
  }

  return(list(par.names = par.names,
              est = est,
              var.cov = v3))
}

#-------------------------------------------------------------------
# Extract item parameter estimates and their covariance matrix from MPlus.
#
# @param mplus.object an object of class \code{mplus.model} from \code{\link[MplusAutomation]{readModels}}. Expected to be a 1-factor model. Mplus .out file must include TECH1 and TECH3.
# @return A three-element \code{list}:
# \itemize{
# \item vector of parameter names taking the form "item.parameter"
# \item list (one element per group) of vectors of item parameter estimates
# \item list (one element per group) of covariance matrices of item parameter estimates
# }
#
# @export
# -------------------------------------------------------------------
#
# get_mplus_params <- function(mplus.object){
#
#   if(length(mplus.object$tech3) == 0 |
#      length(mplus.object$tech1) == 0){
#     stop("mplus.object must include TECH1 and TECH3.
#          TECH3 contains the covariance matrix of item parameter estimates.
#          TECH1 specifies the order of the estimates in TECH3.")
#   }
#
#   ## extracting parameters and item names
#   pars <- mplus.object$parameters$unstandardized # groups in numeric order based on input
#   original.names <- unique(pars$param[grepl("BY$", pars$paramHeader)]) # user name (after Mplus capitalizes it)
#   internal.names <- paste0("item", 1:length(original.names)) # item names to be used for robustDIF functions
#
#   ## check that factor was standardized
#   if(mean(pars[pars$paramHeader == "Variances",]$est) != 1){ # assumes there are no other variances besides the factor
#     stop("Factor variance must be standardized (e.g., f@1).")
#   }
#
#   ## extracting tech1 with parameter specifications for subsetting and sorting vcov matrix
#   t1 <- mplus.object$tech1$parameterSpecification
#
#   if("Group" %in% names(pars)){
#
#     # split dataframe into list with dataframe for each group
#     par.group <- lapply(unique(pars$Group), function(x) pars[pars$Group == x,])
#     par.df <- lapply(par.group, format_pars, names.vec = original.names, type = "mplus")
#     param.est <- lapply(par.df, "[[", "est") # list with each element a vector of estimates
#
#     # defining parameter names
#     param.names = list(internal = "placeholder",
#                        original = paste(par.df[[1]]$item, par.df[[1]]$par, sep = "."))
#
#     # named vector of parameter numbers in desired order
#     # used for subsetting and ordering vcov matrix
#     param.nums <- lapply(t1, function(x){
#       pn <- c(x$tau, t(x$lambda))
#       names(pn) <- c(colnames(x$tau), rownames(x$lambda))
#       pn <- pn[order(factor(names(pn), levels = par.df[[1]]$label))]
#       return(pn)
#     })
#
#   } else {
#
#     par.df <- format_pars(pars, names.vec = original.names, type = "mplus")
#     param.est <- par.df$est # vector of estimates
#
#     # defining parameter names
#     param.names = list(internal = "placeholder",
#                        original = paste(par.df$item, par.df$par, sep = "."))
#
#     # named vector of parameter numbers in desired order
#     # used for subsetting and ordering vcov matrix
#     param.nums <- c(t1$tau, t(t1$lambda))
#     names(param.nums) <- c(colnames(t1$tau), rownames(t1$lambda))
#     param.nums <- param.nums[order(factor(names(param.nums), levels = par.df$label))]
#
#   }
#
#   ## defining parameter names for internal use
#   x <- param.names$original
#   for(i in 1:length(original.names)){
#     x <- sub(original.names[[i]], internal.names[[i]], x) # replacing original name with generic name
#   }
#   param.names$internal <- x
#
#   ## extracting and re-organizing vcov matrix
#   # Order of parameters differs for multigroup and single group; thresholds then slopes vs slopes then thresholds, respectively
#   # need to use tech1 to determine
#   v <- mplus.object$tech3$paramCov
#   v[upper.tri(v)] <- t(v)[upper.tri(v)] # fill in upper triangle
#
#   if(inherits(param.nums, "list")){
#
#     ## subset vcov by group and reorder by item then parameter
#     v2 <- lapply(param.nums, function(x) v[x,x])
#     v2 <- lapply(v2, function(x){
#       row.names(x) <- colnames(x) <- NULL
#       return(x)
#     })
#
#   } else {
#     v2 <- v[param.nums, param.nums] # reorder by item, then parameter rather than the reverse
#     row.names(v2) <- colnames(v2) <- NULL
#   }
#
#   return(list(par.names = param.names,
#               est = param.est,
#               var.cov = v2))
# }


#-------------------------------------------------------------------
#' Helper function used to format parameters estimates
#'
#' @param pars numeric vector of item parameter estimates
#' @param names.vec character vector item names
#' @param type character; are \code{pars} from \code{"lavaan"} or \code{"mplus"}?
#' @return data.frame of item parameter estimates
#' @seealso [robustDIF::get_lavaan_params()], [robustDIF::get_mplus_params()]
# -------------------------------------------------------------------

format_pars <- function(pars, names.vec, type){

  if(type == "lavaan"){

    slopes <- data.frame(item = row.names(pars$lambda),
                         par = "a",
                         est = pars$lambda[,1])

    ints <- data.frame(item = gsub("\\|.*$", "", row.names(pars$tau)),
                       par = gsub("^.*\\|", "", row.names(pars$tau)),
                       est = pars$tau[,1])

  # } else if(type == "mplus"){
  #
  #   # removing unnecessary information
  #   lambda <- pars[grepl("BY$", pars$paramHeader), 2:3]
  #   tau <- pars[pars$paramHeader %in% c("Thresholds","Steps"), 2:3]
  #
  #   slopes <- data.frame(label = lambda$param, # needed to subset and organize vcov
  #                        item = lambda$param,
  #                        par = "a",
  #                        est = lambda$est)
  #
  #   ints <- data.frame(label = tau$param, # needed to subset and organize vcov
  #                      item = gsub("\\$.*$", "", tau$param),
  #                      par = paste0("d",gsub("^.*\\$", "", tau$param)),
  #                      est = tau$est)
  }

  ## binding and sorting parameters
  par.df <- rbind(slopes, ints)
  par.df$item <- factor(par.df$item, levels = names.vec) # make sure items in correct order
  par.df <- par.df[with(par.df, order(item, par)), ]  # sorting parameter estimates
  return(par.df)
}

reformat_out <- function(out){
  out$par.names$internal <- gsub("a1", "a", out$par.names$internal)
  out$par.names$internal <- gsub("a", "a1", out$par.names$internal)
  out$par.names$internal <- gsub("t", "d", out$par.names$internal)
  out$par.names$internal <- gsub("idem", "item", out$par.names$internal)

  # fix! wont work with > 99 items
  end.item.name.index <- unlist(gregexpr(".", out$par.names$internal, fixed = T))
  items <- unique(substr(out$par.names$internal, 1, end.item.name.index))
  n.items <- length(items)
  new.items <- paste0("item", 1:n.items, ".")
  if (any(items != new.items)) {
    for(i in 1:n.items){
      out$par.names$internal <- gsub(items[i], new.items[i],
                                   out$par.names$internal, fixed = T)
    }
  }
  pars.per.item <- lapply(new.items, function(x)grep(x, out$par.names$internal, fixed = T))
  max.pars.per.item <-max(Reduce(c, lapply(pars.per.item, length)))
  par.mat <- as.data.frame(matrix(NA, nrow = n.items, ncol = max.pars.per.item))
  names(par.mat) <- c("a1", paste0("d", 1:(max.pars.per.item-1)))
  row.names(par.mat) <- paste0("item", 1:n.items)
  n.groups <- length(out$est)

  for(i in 1:n.groups){
    for(j in 1:n.items){
      k <- length(pars.per.item[[j]])
      par.mat[j, 1:k] <- out$est[[i]][pars.per.item[[j]]]
    }
    out$est[[i]] <- par.mat
  }

  names(out$est) <- names(out$var.cov) <- paste0("group.", 1:n.groups)
  out
}

