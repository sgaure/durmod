#'

#' Moore-Penrose generalized inverse
#'
#' @param X matrix
#' @param tol tolerance for determining bad entries
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


#' Extract the mixed proportional hazard distribution
#'
#' @description
#' Various functions for extracting the proportional hazard distribution.
#' @param pset
#' a parameter set of class \code{"mphcrm.pset"}, typically
#' \code{opt[[1]]$par}, where \code{opt} is returned from \code{\link{mphcrm}}.
#' @return
#' matrix where there is one row for each masspoint. The first consists of the probabilities,
#' the other columns are the hazards for each transition.
#' @export
mphdist <- function(pset) {
  mus <- exp(sapply(pset$parset, function(pp) pp$mu))
  if(is.null(colnames(mus))) {
    mus <- t(mus)
    colnames(mus) <- names(pset$parset)
  }
  rownames(mus) <- sprintf('point %2d',seq_len(nrow(mus)))
  prob <- a2p(pset$parg)
  cbind(prob,mus)
}

#' Extract the mixed proportional log hazard distribution
#' @rdname mphdist
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
#' @rdname mphdist
#' @export
mphmoments <- function(pset) {
  dist <- mphdist(pset)
  mean <- rowSums(apply(dist, 1, function(x) x[1]*x[-1]))
  variance <- rowSums(apply(dist, 1, function(x) x[1]*(x[-1]-mean)^2))
  sd <- sqrt(variance)
  cbind(mean,variance,sd)
}

#' Extract moments of the mixed proportional log hazard distribution
#' @rdname mphdist
#' @export
mphmoments.log <- function(pset) {
  dist <- mphdist.log(pset)
  mean <- rowSums(apply(dist, 1, function(x) x[1]*x[-1]))
  variance <- rowSums(apply(dist, 1, function(x) x[1]*(x[-1]-mean)^2))
  sd <- sqrt(variance)
  cbind(mean,variance,sd)
}

#' Extract covariance matrix of the proportional hazard distribution
#' @rdname mphdist
#' @export
mphcovs <- function(pset) {
  dist <- mphdist(pset)
  mean <- rowSums(apply(dist, 1, function(x) x[1]*x[-1]))
  crossprod(dist[,1]*(dist[,-1]-rep(mean,echo=nrow(dist))))
}

#' Extract covariance matrix of the proportional hazard distribution
#' @rdname mphdist
#' @export
mphcovs.log <- function(pset) {
  dist <- mphdist.log(pset)
  mean <- rowSums(apply(dist, 1, function(x) x[1]*x[-1]))
  crossprod(dist[,1]*(dist[,-1]-rep(mean,echo=nrow(dist))))
}

#' Extract standard errors of the estimated parameters
#' @param x
#' The Fisher matrix, typically from \code{opt[[1]]$fisher}, where \code{opt} is returned
#' from \code{\link{mphcrm}}.
#' @param tol tolerance for \link{geninv}
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
