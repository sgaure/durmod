#'
#' @export
safeinv <- function(X,eps=1e-6) {
  inv <- try(solve(X), silent=TRUE)
  if(!inherits(inv,'try-error')) return(inv)
  badvars <- match(names(collin(nazero(X))), colnames(X))
  inv <- matrix(Inf,nrow(X),ncol(X),dimnames=dimnames(X))
  inv[-badvars,-badvars] <- solve(X[-badvars,-badvars])
  structure(inv,badvars=badvars)
}

nazero <- function(x) ifelse(is.na(x),0,x)
# Do a Cholesky to detect multi-collinearities
cholx <- function(mat, eps=1e-6) {
  if(is.null(dim(mat))) dim(mat) <- c(1,1)
  N <- dim(mat)[1]
  if(N == 1) {
      return(structure(sqrt(mat),badvars=if(mat<=0) 1 else NULL))
  }

  # first, try a pivoted one
  tol <- N*eps
  chp <- chol(mat,pivot=TRUE,tol=tol)
  rank <- attr(chp,'rank')
  if(rank == N) return(chol(mat))
#  if(rank == N) return(chp)
  pivot <- attr(chp,'pivot')
  oo <- order(pivot)
  badvars <- pivot[((rank+1):N)]
  ok <- (1:N)[-badvars]
  ch <- chol(mat[ok,ok])
  return(structure(ch,badvars=badvars))
}

# cbind(val=flatten(dd$par), sd=sqrt(diag(safeinv(dd$fisher))))

# create a matrix(t1mus, t2mus,... probs)
pdist <- function(pset) {
  mus <- as.matrix(sapply(pset$parset, function(pp) pp$mu))
  rownames(mus) <- paste0('point ',seq_len(nrow(mus)))
  prob <- a2p(pset$parg)
  cbind(prob,mus)
}

pmoments <- function(pset) {
  dist <- pdist(pset)
  mean <- rowSums(apply(dist, 1, function(x) x[1]*x[-1]))
  variance <- rowSums(apply(dist, 1, function(x) x[1]*(x[-1]-mean)^2))
  sd <- sqrt(variance)
  cbind(mean,variance,sd)
}
pmoments.exp <- function(pset) {
  dist <- pdist(pset)
  mean <- rowSums(apply(dist, 1, function(x) x[1]*exp(x[-1])))
  variance <- rowSums(apply(dist, 1, function(x) x[1]*(exp(x[-1])-mean)^2))
  sd <- sqrt(variance)
  cbind(mean,variance,sd)
}

se <- function(x) {
  if(is.matrix(x)) return(sqrt(diag(safeinv(x))))
  if(!is.null(x$fisher)) return(sqrt(diag(safeinv(x$fisher))))
  stop("Can't find a matrix to invert")
}











