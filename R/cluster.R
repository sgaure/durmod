prepcluster <- function(dset,control) {
  # split the dataset to the cluster
  # evenly on spells, but so that they get approx the same number
  # of observations each
  # make it temporarily 1-based, easier in R
  cluster <- control$cluster
  
  spellidx <- dset$spellidx+1
  K <- length(parallel::clusterCall(cluster, function() 1))
  if(K == 1) {message('cluster of size 1 ignored'); return()}
#  K <- 7
  # we should find a set of spells so that dset is evenly divided
  N <- spellidx[length(spellidx)]
  cumeach <- (1:K)/K*N
  ends <- sapply(cumeach, function(i) tail(which(spellidx < i),1))
  starts <- 1L+c(0L,ends[-K])
  spe <- c(spellidx,N+1)
  csplit <- mapply(function(s,e) spellidx[s]:(spe[e+1]-1), starts,ends,SIMPLIFY=TRUE)

  # csplit is a list. One element for each node.
  # Each element is a vector of observations which
  # that node should compute on
  # The dataset is a list with element data, d, timing, id, spellidx, duration, state and riskset
  #
  # entries timing and riskset should just be kept as is
  # 
  # The vectors d, state, id, and duration should be split as observations
  # a new spellidx (zero based) must be created

  # the data element is a list, one element for each transition
  # i.e. data = list(t1=list(mat=..., faclist=...), t2=list(mat=..., faclist=...))
  # mat is a matrix, there is one column for each observation
  # the columns should be split
  # faclist is a list of factors, each factor possibly has an attribute 'x'
  # the factors should be split, and if present, the x
  # the levels must be kept intact
  dsplit <- lapply(csplit, function(obsvec) {
    d <- dset$d[obsvec]
    duration <- dset$duration[obsvec]
    state <- dset$state[obsvec]
    id <- factor(dset$id[obsvec])
    data <- lapply(dset$data, function(dd) {
      mat <- dd$mat[,obsvec,drop=FALSE]
      faclist <- lapply(dd$faclist, function(ff) {
        f <- ff[obsvec]
        x <- attr(ff,'x',exact=TRUE)
        if(!is.null(x)) attr(f,'x') <- x[obsvec]
        f
      })
      m <- list(mat=mat,faclist=faclist)
      attributes(m) <- attributes(dd)
      m
    })
    list(data=data,d=d, timing=dset$timing, id=id, 
         spellidx=NULL, duration=duration, state=state, riskset=dset$riskset)
  })

  # Then add the id and spellidx in each entry
  for(i in seq_along(dsplit)) {
    sidx <- c(spellidx[starts[i]:ends[i]], spellidx[ends[i]+1])
    dsplit[[i]]$spellidx <- sidx - sidx[1]  # zero based
  }

  # push the different datasets to the global environments on the cluster nodes
  parallel::clusterEvalQ(cluster, library(durmod))
  parallel::clusterApply(cluster,dsplit,function(dset) {assign('dset', dset, 
                                                               environment(durmod::.cloglik)); NULL})
  assign('cluster',cluster,environment(mphloglik))
}

cleancluster <- function(cluster) {
  parallel::clusterEvalQ(cluster, {assign('dset',NULL,environment(durmod::.cloglik)); gc()})
  assign('cluster',NULL,environment(mphloglik))
}

mphloglik <- local({
  cluster <- NULL
  function(...) {
    mc <- match.call()
    if(is.null(cluster)) {
      mc[[1L]] <- cloglik
      return(eval.parent(mc))
    }
    mc[[1L]] <- quote(list)
    # ditch dataset
    mc[[2L]] <- NULL  
    # get the other args
    args <- eval.parent(mc)
    args[['control']] <- args[['control']][c('threads','fishblock')]
    # put back dataset
    mc <- as.call(c(list(quote(durmod::.cloglik)),quote(dset), args))
    ret <- parallel::clusterCall(cluster, eval, mc)
    # collect results
    result <- Reduce(`+`,ret)
    # any gradient or fisher, collect them as well
    att <- attributes(ret[[1]])
    for(at in names(att)) {
      attr(result,at) <- Reduce(`+`,lapply(ret,function(r) attr(r,at)))
    }
    result
  }
})

#' @export
.cloglik <- local({
  dset <- NULL
  function(...) {
    mc <- match.call()
    mc[[1L]] <- cloglik
    mc[[2L]] <- dset
    eval.parent(mc)
  }
})
