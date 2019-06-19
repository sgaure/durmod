#'

#' Moore-Penrose generalized inverse
#' 
#' @export
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

#' @export
pdist <- function(pset) {
  mus <- sapply(pset$parset, function(pp) pp$mu)
  if(is.null(colnames(mus))) {
    mus <- t(mus)
    colnames(mus) <- names(pset$parset)
  }
  rownames(mus) <- sprintf('point %2d',seq_len(nrow(mus)))
  prob <- a2p(pset$parg)
  cbind(prob,mus)
}

#' @export
pmoments <- function(pset) {
  dist <- pdist(pset)
  mean <- rowSums(apply(dist, 1, function(x) x[1]*x[-1]))
  variance <- rowSums(apply(dist, 1, function(x) x[1]*(x[-1]-mean)^2))
  sd <- sqrt(variance)
  cbind(mean,variance,sd)
}

#' @export
pmoments.exp <- function(pset) {
  dist <- pdist(pset)
  mean <- rowSums(apply(dist, 1, function(x) x[1]*exp(x[-1])))
  variance <- rowSums(apply(dist, 1, function(x) x[1]*(exp(x[-1])-mean)^2))
  sd <- sqrt(variance)
  cbind(mean,variance,sd)
}

#' @export
se <- function(x,tol=1e-9) {
  if(is.matrix(x)) return(sqrt(diag(geninv(x))))
  if(!is.null(x$fisher)) return(sqrt(diag(geninv(x$fisher))))
  stop("Can't find a matrix to invert")
}

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
