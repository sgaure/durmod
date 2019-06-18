#' @method coef mphcrm.pset
#' @export
coef.mphcrm.pset <- function(object, ...) {
  pars <- unlist(lapply(object$parset, function(pp) {
    c(pp$pars,unlist(pp$facs))
  }))
  mus <- unlist(lapply(object$parset, function(pp) {
    pp$mu
  }))
  probs <- a2p(object$pargs)
  names(probs) <- paste('P',seq_along(probs),sep='')
  c(pars,mus,probs)
}

#' @method print mphcrm.pset
#' @export
print.mphcrm.pset <- function(x, ...) {
  mus <- NULL
  vecs <- NULL
  parset <- x$parset
  for(p in seq_along(parset)) {
    # collect parameters
    tn <- names(parset)[p]
    pp <- parset[[p]]
    vec <- pp$pars
    names(vec) <- paste(names(vec), tn,sep='_')
    facs <- unlist(pp$facs)
    if(!is.null(facs))
      names(facs) <- sapply(seq_along(facs), function(ff) paste(names(facs)[ff],levels(facs[ff]),'_',tn,sep=''))
    vecs <- c(vecs,vec,facs)
    mu <- pp$mu
    names(mu) <- paste(names(mu), tn, sep='_')
    vecs <- c(vecs,mu)
  }
  probs <- a2p(x$pargs)
  names(probs) <- paste('P',seq_along(probs),sep='')
  print(data.frame(vecs))
  print(probs)
}

#' @export
flatten <- function(x, exclude=attr(x,'exclude')) {
  class(x) <- setdiff(class(x),'mphcrm.pset')
  vec <- unlist(x)
  # fix the names so more readable
  # parset.t1.{pars,mu,facs}.x <- t1.x
  newnames <- gsub('^parset\\.(.*)\\.(pars|mu|facs)\\.(.*)','\\1.\\3',names(vec))
  names(vec) <- newnames

  skel <- attr(vec,'skeleton')
  class(skel) <- c('mphcrm.pset',class(skel))
  attr(vec,'skeleton') <- skel
  if(length(exclude) == 0) return(vec)
  # remove the excluded parameters
  exc <- sapply(exclude, function(ee) ee[[1]])
  if(is.character(exc)) exc <- which(names(vec) %in% exc)
  exc <- intersect(exc,seq_along(vec))
  if(length(exc)==0) return(structure(vec,exclude=NULL))
  vals <- vec[exc]
  vec <- vec[-exc]
  attr(vec,'skeleton') <- skel
  exc <- structure(mapply(function(a,b) list(a,b), exc, vals,SIMPLIFY=FALSE), names=names(vals))
  structure(vec,exclude=exc)
}

#' @export
unflatten <- function(flesh, skeleton=attr(flesh, 'skeleton'), exclude=attr(flesh,'exclude')) {
  class(skeleton) <- setdiff(class(skeleton),'mphcrm.pset')
  ne <- length(exclude)
  if(ne == 0) {
    pset <- relist(flesh,skeleton)
    class(pset) <- c('mphcrm.pset',class(pset))
    return(pset)
  }

  # insert the excluded parameters before relisting
  exc <- sapply(exclude, function(ee) ee[[1]])
  vals <- sapply(exclude, function(ee) ee[[2]])
  vec <- numeric(length(flesh)+ne)
  idx <- seq_along(vec)
  ok <- setdiff(idx,exc)
  vec[ok] <- flesh
  vec[exc] <- vals
  names(vec)[ok] <- names(flesh)
  names(vec)[exc] <- names(vals)
  pset <- relist(vec,skeleton)
  class(pset) <- c('mphcrm.pset',class(pset))
  structure(pset,exclude=as.list(exc))
}

exclude <- function(x,exclude,...) {
  # the generic doesn't have a ... argument, so call it directly
  relist(getS3method('unlist','mphcrm.pset')(x,exclude=exclude))
}
