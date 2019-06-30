#'

#' Moore-Penrose generalized inverse
#'
#' @param X matrix
#' @param tol tolerance for determining bad entries
#' @examples
#' # create a positive definite 5x5 matrix
#' x <- crossprod(matrix(rnorm(25),5))
#' # make it singular
#' x[,2] <- x[,3]+x[,5]
#' geninv(x)
#' @export
geninv <- function(X, tol=.Machine$double.eps) {
  if (length(dim(X)) > 2L || !is.numeric(X)) 
    stop("'X' must be a numeric matrix")
  if (!is.matrix(X))  X <- as.matrix(X)
  nm <- colnames(X)
  Xsvd <- svd(X)
  Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)
  inv <- if (all(Positive)) 
    Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u))
  else if (!any(Positive)) 
    array(NA, dim(X)[2L:1L])
  else {
    im <- Xsvd$v[, Positive, drop = FALSE] %*% ((1/Xsvd$d[Positive]) * 
                                                 t(Xsvd$u[, Positive, drop = FALSE]))
    im[!Positive,!Positive] <- NA
    im
  }
  structure(inv,badvars=nm[!Positive], dimnames=list(nm,nm))
}

nazero <- function(x) ifelse(is.na(x),0,x)


#' Extract the mixed proportional hazard distribution
#'
#' @description
#' Various functions for extracting the proportional hazard distribution.
#'
#' \code{mphdist} extracts the hazard distribution.
#' @param pset
#' a parameter set of class \code{"mphcrm.pset"}, typically
#' \code{opt[[1]]$par}, where \code{opt} is returned from \code{\link{mphcrm}}.
#' If given a list of results, extracts the
#' first in the list.
#' @return
#' A matrix.
#' @examples
#' # load a dataset and a precomputed fitted model
#' data(durdata)
#' best <- fit[[1]]
#' mphdist(best)
#' mphmoments(best)
#' mphcov.log(best)
#' @export
mphdist <- function(pset) {
  if(inherits(pset,'mphcrm.list')) pset <- pset[[1]]$par
  if(inherits(pset,'mphcrm.opt')) pset <- pset$par
  mus <- exp(sapply(pset$parset, function(pp) pp$mu))
  if(is.null(colnames(mus))) {
    mus <- t(mus)
    colnames(mus) <- names(pset$parset)
  }
  rownames(mus) <- sprintf('point %2d',seq_len(nrow(mus)))
  prob <- a2p(pset$pargs)
  cbind(prob,mus)
}

#' Extract the mixed proportional log hazard distribution
#' @rdname mphdist
#' @description \code{mphdist.log} extracts the log hazard distribution.
#' @export
mphdist.log <- function(pset) {
  if(inherits(pset,'mphcrm.list')) pset <- pset[[1]]$par
  if(inherits(pset,'mphcrm.opt')) pset <- pset$par
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
#' @rdname mphdist
#' @description
#' \code{mphmoments} returns the first and second moments of the hazard distribution.
#' @export
mphmoments <- function(pset) {
  dist <- mphdist(pset)
  mean <- colSums(dist[,1]*dist[,-1,drop=FALSE])
  if(length(mean) == 1) 
    variance <- sum(dist[,1]*(dist[,-1]-mean)^2) 
  else 
    variance <- rowSums(apply(dist, 1, function(x) x[1]*(x[-1]-mean)^2))
  sd <- sqrt(variance)
  cbind(mean,variance,sd)
}
#' Extract moments of the mixed proportional log hazard distribution
#' @rdname mphdist
#' @description
#' \code{mphmoments.log} returns the first and second moments of the log hazard distribution.
#' @export
mphmoments.log <- function(pset) {
  dist <- mphdist.log(pset)
  mean <- colSums(dist[,1]*dist[,-1,drop=FALSE])
  if(length(mean) == 1) 
    variance <- sum(dist[,1]*(dist[,-1]-mean)^2) 
  else
    variance <- rowSums(as.matrix(apply(dist, 1, function(x) x[1]*(x[-1]-mean)^2)))
  sd <- sqrt(variance)
  cbind(mean,variance,sd)
}

#' Extract covariance matrix of the proportional hazard distribution
#' @rdname mphdist
#' @description
#' \code{mphcov} returns the variance/covariance matrix of the hazard distribution.
#' @export
mphcov <- function(pset) {
  dist <- mphdist(pset)
  mean <- colSums(dist[,1]*dist[,-1,drop=FALSE])
  crossprod(sqrt(dist[,1])*(dist[,-1,drop=FALSE]-rep(mean,each=nrow(dist))))
}



#' Extract covariance matrix of the proportional hazard distribution
#' @rdname mphdist
#' @description
#' \code{mphcov.log} returns the variance/covariance matrix of the log hazard distribution.
#' @export
mphcov.log <- function(pset) {
  dist <- mphdist.log(pset)
  mean <- colSums(dist[,1]*dist[,-1,drop=FALSE])
  crossprod(sqrt(dist[,1])*(dist[,-1,drop=FALSE]-rep(mean,each=nrow(dist))))
}

#' Extract standard errors of the estimated parameters
#' @param x
#' The Fisher matrix, typically from \code{opt[[1]]$fisher}, where \code{opt} is returned
#' from \code{\link{mphcrm}}.
#' @param tol tolerance for \link{geninv}
#' @export
se <- function(x,tol=.Machine$double.eps) {
  if(is.matrix(x)) return(sqrt(diag(geninv(x,tol))))
  if(!is.null(x$fisher)) return(sqrt(diag(geninv(x$fisher,tol))))
  stop("Can't find a matrix to invert")
}

#' Prettyprint a time interval
#' @description
#' Prints a time in seconds as e.g. \code{"3m4s"}.
#' @param t
#' numeric. time in seconds.
#' @examples
#' timestr(1.3)
#' timestr(73)
#' timestr(4684)
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
