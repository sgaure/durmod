#'

##' Moore-Penrose generalized inverse
##' 
##' @export
geninv <- function(X, tol=sqrt(.Machine$double.eps)) {
  if (length(dim(X)) > 2L || !is.numeric(X)) 
    stop("'X' must be a numeric matrix")
  if (!is.matrix(X))  X <- as.matrix(X)
  nm <- colnames(X)
  Xsvd <- svd(X)
  Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)
  inv <- if (all(Positive)) 
    Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u))
  else if (!any(Positive)) 
    array(0, dim(X)[2L:1L])
  else {
    im <- Xsvd$v[, Positive, drop = FALSE] %*% ((1/Xsvd$d[Positive]) * 
                                                 t(Xsvd$u[, Positive, drop = FALSE]))
    im[!Positive,!Positive] <- NA
    im
  }
  structure(inv,badvars=nm[!Positive], dimnames=list(nm,nm))
}

nazero <- function(x) ifelse(is.na(x),0,x)

# cbind(val=flatten(dd$par), sd=sqrt(diag(safeinv(dd$fisher))))

# create a matrix(t1mus, t2mus,... probs)


#' Extract the mixed proportional hazard distribution
#'
#' @description
#' 
#' @param pset
#' a parameter set of class \code{"mphcrm.pset"}, typically
#' \code{opt[[1]]$par}, where \code{opt} is returned from \code{\link{mphcrm}}.
#' @param round
#' integer. Number of decimals. In a competing risk model, some of the hazards can be very small, \code{mphdist}
#' can round the value so that e.q. \code{1e-8} becomes zero, mostly useful for nice displays.
#' @return
#' matrix where there is one row for each masspoint. The first consists of the probabilities,
#' the other columns are the hazards for each transition.
#' @export
mphdist <- function(pset,round=20) {
  mus <- round(exp(sapply(pset$parset, function(pp) pp$mu)),round)
  if(is.null(colnames(mus))) {
    mus <- t(mus)
    colnames(mus) <- names(pset$parset)
  }
  rownames(mus) <- sprintf('point %2d',seq_len(nrow(mus)))
  prob <- a2p(pset$parg)
  cbind(prob,mus)
}

#' @rdname mphdist
#' Same as \code{mphdist}, but the logarithms of the hazards.
#' @export
mphdist.log <- function(pset) {
  mus <- sapply(pset$parset, function(pp) pp$mu)
  if(is.null(colnames(mus))) {
    mus <- t(mus)
    colnames(mus) <- names(pset$parset)
  }
  rownames(mus) <- sprintf('point %2d',seq_len(nrow(mus)))
  prob <- a2p(pset$parg)
  cbind(prob,mus)
}

#' Extract moments of the mixed proportional hazard distribution
#' @seealso @mphdist
#' @export
mphmoments <- function(pset,round=20) {
  dist <- mphdist(pset,round)
  mean <- rowSums(apply(dist, 1, function(x) x[1]*x[-1]))
  variance <- rowSums(apply(dist, 1, function(x) x[1]*(x[-1]-mean)^2))
  sd <- sqrt(variance)
  cbind(mean,variance,sd)
}

#' @rdname mphmoments
#' @export
mphmoments.log <- function(pset) {
  dist <- mphdist.log(pset)
  mean <- rowSums(apply(dist, 1, function(x) x[1]*x[-1]))
  variance <- rowSums(apply(dist, 1, function(x) x[1]*(x[-1]-mean)^2))
  sd <- sqrt(variance)
  cbind(mean,variance,sd)
}


#' Extract standard errors of the estimated parameters
#' @param x
#' The Fisher matrix, typically from \code{opt[[1]]$fisher}, where \code{opt} is returned
#' from \code{\link{mphcrm}}.
#' @export
se <- function(x,tol=1e-9) {
  if(is.matrix(x)) return(sqrt(diag(geninv(x))))
  if(!is.null(x$fisher)) return(sqrt(diag(geninv(x$fisher))))
  stop("Can't find a matrix to invert")
}

#' Prettyprint a time interval
#' @description
#' Prints a time in seconds as e.g. \code{"3m4s"}.
#' @param t
#' numeric. time in seconds.
#' @export
timestr <- function(t) {
  dec <- t - as.integer(t)
  t <- as.integer(t-dec)
  s <- t %% 60
  t <- t %/% 60
  m <- t %% 60
  h <- t %/% 60
  if(h > 0) {
    str <- sprintf('%dh%dm',h,m)
  } else if(m > 0) {
    str <- sprintf('%dm%.0fs',m,s)
  } else {
    str <- sprintf('%.1fs',s+dec)
  }
  str
}
