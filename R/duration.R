#' Estimate a mixed proportional hazard model
#'
#' @description
#' \code{mphcrm} implements estimation of a mixed proportional hazard competing risk model.
#' The baseline hazard is of the form \eqn{exp(X \beta)} where \eqn{X} is a
#' matrix of covariates, and \eqn{\beta} is a vector of paramaters to estimate.
#' In addition there is an intercept term \eqn{\mu}, i.e. the hazard is \eqn{exp(X \beta + \mu)}.
#' There are several transitions to be made, and a set of (\eqn{X}, \eqn{\beta}, and \eqn{\mu}) for
#' each possible transition.
#'
#' Each individual may have several observations, with either a transition at the end of the
#' observation, or not a transition. It is a competing risk, there can be more than one possible
#' transition for an observation, but only one is taken at the end of the period.
#' 
#' For each individual \eqn{i} there is a log likelihood as a function of \eqn{\mu}, called \eqn{M_i(\mu)}.
#'
#' The mixture is that the \eqn{\mu}'s are stochastic. I.e. 
#' we have probabilities \eqn{p_j}, and a vector of \eqn{\mu_j} of masspoints (one for each transition),
#' for each such \eqn{j}.
#'
#' So the full likelihood for an individual is \eqn{L_i = \sum_j p_j M_i(\mu_j)}.
#'
#' The \code{mphcrm()} function maximizes the likelihood \eqn{\sum_i L_i} over
#' \eqn{p_j}, \eqn{mu_j}, and \eqn{\beta}.
#'
#' In addition to the parameters \eqn{\beta}, a variable which records the duration of each
#' observation must be specified.
#'
#' In some datasets it is known that not all risks are present at all times. Like, losing your
#' job when you do not have one. In this case it should be specified which risks are present.
#'
#' The estimation starts out with one masspoint, maximizes the likelihood, tries to add another
#' point, and continues in this fashion.
#' 
#' @param formula
#' A formula specifying the covariates.
#' In a formula like \code{d ~ x1 + x2},
#' the \code{d} is the transition which is taken, coded as an integer where \code{0} means
#' no transition, and otherwise \code{d} is the number of the transition which is taken.
#' \code{d} can also be a factor, in which case the level which is no transition must be named
#' \code{"0"} or \code{"none"}.
#'
#' The \code{x1+x2} part is like in \code{\link{lm}}, i.e. ordinary covariates or factors.
#'
#' If the covariates differ among the transitions, a multi-part formula can be specified.
#' E.g. \code{d ~ x1+x2 | x3+x4 | x5+x6}, in which case the first part (\code{x1+x2}) is common
#' for all the transitions, the second part, \code{x3+x4}, is specific to the first transition,
#' the third part is specific for the second transition, and so on. I.e. the covariates
#' for transition 1 are \code{x1+x2+x3+x4}, whereas the covariates for transition 2 are \code{x1+x2+x5+x6}.
#' 
#' @param data
#' A data frame which contains the covariates.
#' @param id
#' The name of the covariate which identifies an individual.
#' @param durvar
#' The name of the covariate which contains the length of the observation period. Can be a
#' numeric, in which case all observations are assumed to have the same length. If missing, it
#' is assumed to be 1.
#' @param state
#' The name of a variable which is an index into \code{risksets}, specifying for each observation
#' which risks are present.
#' @param risksets
#' A list of integer vectors. Each vector is a list of transitions, which risks are present for
#' the observation.
#' @param subset
#' For specifying a subset of the dataset, similar to \code{\link{lm}}.
#' @param na.action
#' For handling of missing cases, similar to \code{\link{lm}}.
#' @param control
#' List of control parameters for the estimation. See \code{\link{mphcrm.control}}.
#' @return
#' [TBS]
#' @details
#' [TBS]
#' @export
mphcrm <- function(formula,data,id,durvar,state,risksets=NULL,
                   timing=c('exact','interval','logit'),
                   subset, na.action, control=mphcrm.control()) {
  timing <- match.arg(timing)
  F <- Formula::as.Formula(formula)
  if(missing(id)) stop('id must be specified')
  environment(F) <- environment(formula)
  # d ~ x + y | z | w ...

  # create model frame
  mf <- match.call(expand.dots=TRUE)
  m <- match(c('formula','data','subset','na.action'), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(model.frame)
  mf[[2L]] <- F

  # add id, durvar and state to formula to get them into model frame to make subset work
  form <- F
  id.v <- bquote(I(.(substitute(id))))
  form <- update(form, as.formula(bquote(. ~ . + .(id.v))))
  durvar.v <- bquote(I(.(substitute(durvar))))
  form <- update(form, as.formula(bquote(. ~ . + .(durvar.v))))
  if(!missing(state)) {
    state.v <- bquote(I(.(substitute(state))))
    form <- update(form, as.formula(bquote(. ~ . + .(state.v))))
  }
  mf[[2L]] <- form
  mf <- eval.parent(mf)
  
  N <- nrow(mf)

  spec <- parseformula(F,mf)

  attr(spec,'timing') <- timing
  id <- factor(eval(as.name(deparse(id.v)),mf,environment(F)))

  if(length(id) != N) stop('id must have length ',N,' (',as.character(id),')')
  if(length(unique(id)) != length(rle(as.integer(id))$values)) {
    stop('dataset must be sorted on id')
  }
  attr(spec,'id') <- id
  # zero-based index of beginning of spells. padded with one after the last observation
  attr(spec,'spellidx') <- c(0,which(diff(as.integer(id))!=0),length(id))
  attr(spec,'maxspellen') <- max(diff(attr(spec,'spellidx')))
  
  if(missing(durvar)) {
    duration <- 1
  } else {
    duration <- as.numeric(eval(as.name(deparse(durvar.v)),mf,environment(F)))
  }

  if(length(duration) != N && length(duration) != 1) {
    stop('Duration must either be a constant or a vector of length ',N,' (',as.character(durvar),')')
  }
  if(length(duration) == 1) duration <- rep(duration,N)
  attr(spec,'duration') <- duration

  hasriskset <- !is.null(risksets)
  if(hasriskset && missing(state)) warning("riskset is specified, but no state")
  if(missing(state)) {
    state <- 0L
    hasriskset <- FALSE
  } else {
    state <- eval(as.name(deparse(state.v)),mf,environment(F))
  }
  
  if(hasriskset) {
    srange <- range(state)
    if(srange[1] != 1) stop('smallest state must be 1 (index into riskset)')
    if(srange[2] != length(risksets)) {
      stop(sprintf('max state is %d, but there are %d risksets',state[2],length(risksets)))
    }
    # check if transitions are taken which are not in the riskset
    if(0 %in% unlist(risksets)) stop("0 can't be in risk set")
    d <- attr(spec,'d')
    badtrans <- mapply(function(dd,r) dd != 0 && !(dd %in% r), d, risksets[state])
    if(any(badtrans)) {
      n <- which(badtrans)[1]
      stop(sprintf("In observation %d(id=%d), a transition to %d is taken, but the riskset of the state(%d) does not allow it",
                   n, id[n], d[n], state[n]))
    }
    attr(spec,'state') <- state
    attr(spec,'riskset') <- risksets
  } else {
    attr(spec,'state') <- 0
    attr(spec,'riskset') < list()
  }
  
    

  pset <- makeparset(spec,1)
  pointiter(spec,pset,control)
}

#' Control parameters for mphcrm
#'
#' @description
#' Modify the default estimation parameters for \code{\link{mphcrm}}
#' @param ...
#' parameters to change \itemize{
#'   \item \code{callback}
#'   A user-specified \code{function(fromwhere, opt, spec)} which is called after each optimization step.
#'   It can be used to report what is happening, check whatever it wants, and optionally stop
#'   the estimation by calling stop(). In this case, \code{mphcrm()} will return with the
#'   most recently estimated set of parameters. See the help on \code{\link{callback_default}} for information
#'   on the argument.
#' }
#' @return
#' list of control parameters suitable for the \code{control} argument of \code{\link{mphcrm}}
#' @details
#' [TBS]
#' @export
mphcrm.control <- function(...) {
  ctrl <- list(iters=12,threads=getOption('durmod.threads'),gradient=TRUE, fisher=TRUE, hessian=FALSE, 
               gdiff=TRUE, minprob=1e-6, eqtol=1e-3, callback=callback_default)
  args <- list(...)
  nm <- names(ctrl)[names(ctrl) %in% names(args)]
  ctrl[nm] <- args[nm]
  ctrl
}

#' Default callback function for mphcrm
#' @export
callback_default <- function(fromwhere,opt,spec) {
  
  if(is.logical(opt$convergence) && opt$convergence != 0) {
    message(date(), ' ', fromwhere, ' optimization failed to converge: ',opt$message,'(',opt$convergence,')')
  }
#  if(!identical(fromwhere, 'full')) return()
  if(!identical(fromwhere, 'full')) {message(date(), ' optimized ',fromwhere);return(TRUE)}
  if(!is.null(opt$fisher)) {
    rc <- rcond(opt$fisher)
  } else {
    rc <- NA
  }
  if(!is.null(opt$gradient)) {
    grd <- sqrt(sum(opt$gradient^2))
  } else {
    grd <- NA
  }
  probs <- a2p(opt$par$pargs)
  message(date(), sprintf(' loglik %.4f grad %.3g minprob %.4g rcond %.2g (%d points)',
                          -opt$value, grd, min(probs), rc, length(opt$par$pargs)+1))
  if(rc < 1e-12) {message(sprintf('  ***bad condition: %.2e',rc)); return(FALSE)}
  return(TRUE)
}


pointiter <- function(spec,pset,control) {
  # optimize null model first
  arg0 <- flatten(pset)
  arg0[] <- 0
  pset <- unflatten(arg0)
  LL0 <- function(arg,spec,pset) {
    for(i in seq_along(pset$parset)) {
      pset$parset[[i]]$mu[] <- arg[i]
    }
    -cloglik(spec,pset,nthreads=control$threads)
  }
  opt0 <- optim(runif(length(pset$parset),-4,-1), LL0,method='BFGS',spec=spec,pset=pset)
#  message('zero model: '); print(opt0$par)
  for(i in seq_along(pset$parset)) {
    pset$parset[[i]]$mu[] <- opt0$par[i]
  }
  opt0$par <- pset
  control$callback('nullmodel',opt0,spec)
  intr <- FALSE
  iopt <- opt <- opt0 
  tryCatch(
    for(i in seq_len(control$iters)) {
      opt <- ml(spec,pset,control)
      control$callback('full',opt, spec)
      pset <- addpoint(spec,opt$par,opt$value,control)
    },
    error=function(...) {
      intr <<- TRUE; iopt <<- opt; 
      warning('error occured, returning most recent estimate: ',...)
    },
    interrupt=function(...) {
      intr <<- TRUE
      iopt <<- opt
      warning("estimation interrupted, returning most recent estimate")
    },
    finally=if(intr) opt <- iopt)
  opt
}

optprobs <- function(spec,pset,control) {
  pfun <- function(a) {
    pset$pargs[] <- a
    -cloglik(spec,pset,nthreads=control$threads)
  }
#  message('optimize probs')
  aopt <- optim(pset$pargs,pfun,method='BFGS',control=list(trace=0,REPORT=1,factr=10))
#  message('probs:',sprintf(' %.7f',a2p(aopt$par)), ' value: ',aopt$value)
  pset$pargs[] <- aopt$par
  control$callback('prob',aopt,spec)
  pset
}

optdist <- function(spec,pset,control) {
  dfun <- function(a) {
    n <- length(pset$pargs)+1
    pset$pargs[] <- a[seq_len(n-1)]
    pos <- n
    for(i in seq_along(pset$parset)) {
      pset$parset[[i]]$mu[] <- a[pos:(pos+n-1)]
      pos <- pos+n
    }
    -cloglik(spec,pset,nthreads=control$threads)
  }
  args <- c(pset$pargs,sapply(pset$parset, function(pp) pp$mu))
#  message('optimize dist')
  dopt <- optim(args,dfun,method='BFGS',control=list(trace=0,REPORT=100))
  n <- length(pset$pargs)
#  message('new probs',sprintf(' %.7f',a2p(dopt$par[1:n])), ' value: ',dopt$value)
  pset$pargs[] <- dopt$par[seq_len(n)]
  pos <- n+1
  for(i in seq_along(pset$parset)) {
    pset$parset[[i]]$mu[] <- dopt$par[pos:(pos+n)]
    pos <- pos+n+1
  }
  dopt$par <- pset
  control$callback('dist',dopt,spec)
  pset
}

newpoint <- function(spec,pset,value,control) {
  gdiff <- control$gdiff
  pr <- a2p(pset$pargs)
  newprob <- if(gdiff) 0 else min(0.1*pr,1e-4)
  newpr <- c( (1 - newprob)*pr, newprob)
#  newpr <- newpr/sum(newpr)
  #  newpr <- c(pr,0)
  newset <- makeparset(spec,length(pr)+1,pset)
  newset$pargs[] = p2a(newpr)
  np <- length(newpr)
  ntrans <- length(pset$parset)
  fun <- function(mu) {
    for(i in seq_along(mu)) {
      newset$parset[[i]]$mu[np] <- mu[i]
    }
    -cloglik(spec,newset,gdiff=gdiff,nthreads=control$threads)
  }
  args <- runif(length(newset$parset),-10,1)
  muopt <- nloptr::nloptr(args, fun, lb=rep(-10,ntrans), ub=rep(1,ntrans),
                          opts=list(algorithm='NLOPT_GN_ISRES',stopval=if(gdiff) 0 else value-1,
                                    maxtime=120,maxeval=10000,population=10*length(args)))
  muopt$convergence <- muopt$status
  control$callback('newpoint',muopt,spec)
  if(!(muopt$status %in% c(0,2))) 
    message(muopt$status, ' ', muopt$message,' (',muopt$objective,') iter ',muopt$iterations)
  
  for(i in seq_along(newset$parset)) {
    newset$parset[[i]]$mu[np] <- muopt$solution[i]
  }
  if(gdiff) {
    newpr <- c((1-1e-5)*pr,1e-5)
    #  newpr <- newpr/sum(newpr)
    newset$pargs[] = p2a(newpr)
  }
  optprobs(spec,newset,control)
}

badpoints <- function(pset,control) {
  # are any masspoints equal?
  np <- length(pset$pargs)+1
  p <- a2p(pset$pargs)
  okpt <- p > control$minprob
  for(i in seq_len(np)) {
    mui <- sapply(pset$parset, function(pp) pp$mu[i])
    for(j in seq_len(i-1)) {
      if(!okpt[j]) next
      muj <- sapply(pset$parset, function(pp) pp$mu[j])
      if(max(abs(exp(muj) - exp(mui))) < control$eqtol) {
        message(sprintf('points %d and %d are equal',j,i))
        mumat <- rbind(muj,mui)
        rownames(mumat) <- c(j,i)
        colnames(mumat) <- names(pset$parset)
        print(cbind(prob=c(p[j],p[i]),mumat))
        okpt[i] <- FALSE
      }
    }
  }
  if(!all(okpt)) message('remove points', sprintf(' %.2e',p[!okpt]))
  p <- p[okpt]

  pset$pargs <- p2a(p)
  for(i in seq_along(pset$parset)) {
    pset$parset[[i]]$mu <- pset$parset[[i]]$mu[okpt]
  }
  
  structure(pset,badremoved=!all(okpt))
}

addpoint <- function(spec,pset,value,control) {
  # find a new point
  # optimize probabilities
  newset <- pset
  repeat {
    newset <- newpoint(spec,newset,value,control)
    newset <- optdist(spec,newset,control)
    newset <- badpoints(newset,control)
    if(!isTRUE(attr(newset,'badremoved'))) break
  } 
  newset
}

ml <- function(spec,pset,control) {
  args <- flatten(pset)
  skel <- attr(args,'skel')
  LL <- function(args,spec,skel) {
    pset <- unflatten(args,skel)
    -cloglik(spec,pset,nthreads=control$threads)
  }
  gLL <- function(args,spec,skel) {
    pset <- unflatten(args,skel)
    -attr(cloglik(spec,pset,dogradient=TRUE,nthreads=control$threads),'gradient')
  }

  fLL <- function(args, spec, skel) {
    pset <- unflatten(args,skel)
    -attr(cloglik(spec,pset,dofisher=TRUE,nthreads=control$threads),'fisher')
  }
  hLL <- function(args, spec, skel) {
    numDeriv::jacobian(gLL,args,spec=spec,skel=skel)
  }

  jdLL <- function(args, spec, skel) {
    numDeriv::jacobian(dLL,args,spec=spec,skel=skel)
  }
  dhLL <- function(args, spec, skel) {
    numDeriv::hessian(LL,args,spec=spec,skel=skel)
  }

  dLL <- function(args,spec,skel) {
    numDeriv::grad(LL,args,spec=spec,skel=skel)
  }

#  args[16:18] <- rnorm(3)
#  args[5:7] <- rnorm(3)

#  message('Hessian matrix:'); print(hess)
#  message('Jacobian matrix:'); print(jac)
#  print(system.time(hess <- dhLL(args,spec,skel)))
#  message('hess eigs:');print(eigen(hess,only.values=TRUE)$values)
#  print(system.time(jac <- hLL(args,spec,skel)))
#  message('jac eigs:');print(eigen(jac,only.values=TRUE)$values)
#  print(system.time(jacnum <- jdLL(args,spec,skel)))
#  message('jacnum eigs:');print(eigen(jacnum,only.values=TRUE)$values)
#  print(system.time(fish <- fLL(args,spec,skel)))
#  message('fish eigs:');print(eigen(fish,only.values=TRUE)$values)
#  message('Fisher matrix:'); print(system.time(print(fLL(args,spec,skel))))
  message('gLL: '); print(system.time(print(gLL(args,spec,skel))))
  message('dLL: '); print(system.time(print(dLL(args,spec,skel))))
#  stop('debug')
#  dLL <- function(args,spec,skel) numDeriv::grad(LL,args,method='simple',spec=spec,skel=skel)
  opt <- optim(args,LL,gLL,method='BFGS',
               control=list(trace=0,REPORT=10,maxit=500,lmm=60,abstol=1e-4,reltol=1e-14),
               spec=spec,skel=skel)
  sol <- unflatten(opt$par)

  # reorder the masspoints, highest probability first
  probs <- a2p(sol$pargs)
  oo <- order(probs,decreasing=TRUE)
  sol$pargs <- p2a(probs[oo])
  sol$parset <- lapply(sol$parset, function(pp) {
    nm <- names(pp$mu)
    pp$mu <- pp$mu[oo]
    names(pp$mu) <- nm
    pp
  })
  opt$par <- sol

  # create a covariance matrix, we use the inverse of the (negative) fisher matrix, it's fast to compute
  nm <- names(flatten(opt$par))

  if(isTRUE(control$gradient)) {
    opt$gradient <- gLL(flatten(opt$par),spec,skel)
  }
  if(isTRUE(control$fisher)) {
    opt$fisher <- fLL(flatten(opt$par),spec,skel)
    dimnames(opt$fisher) <- list(nm,nm)
  }
  if(isTRUE(control$hessian)) {
    opt$hessian <- hLL(flatten(opt$par),spec,skel)
    dimnames(opt$hessian) <- list(nm,nm)
  }
#  opt$vcov <- try(solve(opt$fisher), silent=TRUE)

  # standard errors (we should adjust for the dof, #spells or #observations? I think spells.)
#  opt$se <- try(unflatten(sqrt(diag(opt$vcov)), skel=skel), silent=TRUE)
#  opt$hvcov <- solve(dhLL(opt$par,spec,skel))
#  opt$hse <- unflatten(sqrt(diag(opt$hvcov)), skel=skel)
  opt
}

#' @export
a2p <- function(a) {b <- c(0,a); p <- exp(b)/sum(exp(b)); ifelse(is.na(p),1,p)}
##' @export
#a2logp <- function(a) {b <- c(0,a); logp <- b - logsumofexp(b,0); ifelse(is.na(logp),0,logp)}

#' @export
p2a <- function(p) log(p/p[1])[-1]

makeparset <- function(spec,npoints,oldset) {
  parset <- lapply(spec, function(tt) {
    list(pars=structure(rep(0,nrow(tt$mat)),names=rownames(tt$mat)),
         facs=lapply(tt$faclist, function(ff) structure(rep(0,nlevels(ff)), names=levels(ff))),
         mu=structure(rep(0,npoints), names=paste('mu',1:npoints,sep=''))
         )
  })
  newset <- as.relistable(list(parset=parset, pargs=rnorm(npoints-1)))
  if(!missing(oldset)) {
    oldset <- flatten(oldset)
    newset <- flatten(newset)
    common <- intersect(names(oldset),names(newset))
    newset[common] <- oldset[common]
    newset <- unflatten(newset)
  }
  class(newset) <- c('mphcrm.pset',class(newset))
  newset
}


parseformula <- function(formula,mf) {
  orig.d <- model.response(mf)
  df <- as.factor(orig.d)
  # put the zero/none/null/0 level first for conversion to integer
  # i.e. ensure that the no-transition is zero
  L <- levels(df)
  nlpos <- L %in% c('0','none')
  ztrans <- sum(nlpos) > 0
  if(sum(nlpos) > 1) stop("There can't be both a 0 and \"none\" in the outcome")
  newlev <- if(ztrans) c(L[nlpos],L[!nlpos]) else L
  df <- factor(df,levels=newlev)
  if(!is.factor(orig.d)) levels(df) <- paste('t',levels(df),sep='')
  d <- as.integer(df)-ztrans

  nrhs <- length(formula)[2]
  if(nrhs > nlevels(df)) stop("There are more parts in the formula than there are transitions")


  # analyze the formula to figure out which
  # covariates explain which transitions
  # split off factors

  transitions <- nlevels(df)-ztrans
  mt <- terms(mf)
  cls <- attr(mt,'dataClasses')

  # for each transition, make a model matrix for the numeric covariates,
  # and a list of factors
  spec <- lapply(2:(transitions+1), function(t) {
    rhs <- c(1,t)
    if(t > nrhs) rhs <- 1
    ff <- Formula::as.Formula(formula(formula,rhs=rhs,lhs=0))
    mt <- terms(ff,keep.order=TRUE)
    fact <- attr(mt,'factors')
    # now, filter out those terms which are neither factors, nor interactions with factors

    keep <- unlist(lapply(colnames(fact), function(lab) {
      contains <- rownames(fact)[fact[,lab] > 0]
      if(!any(cls[contains]=='factor')) return(lab)
      return(NULL)
    }))

    attr(mt,'factors') <- fact[,keep,drop=FALSE]
    attr(mt,'intercept') <- 0
    mat <- model.matrix(mt,mf)

# then the factor related stuff

    faclist <- lapply(colnames(fact), function(term) {
      codes <- fact[,term]

      contains <- rownames(fact)[codes > 0]
      # are there any factors in this term?
      isfac <- cls[contains] == 'factor'
      if(!any(isfac)) return(NULL)

      # interact all the factors in the term
      # but some are with contrasts, some are not. The code is 2 if no contrast, 1 otherwise
      # make a list of the factors

      codes <- codes[contains[isfac]]
      flist <- eval(as.call(lapply(c('list',contains[isfac]),as.name)), mf,environment(formula))
#      names(flist) <- contains[isfac]
      # remove a reference level if contrasts
      fl <- mapply(function(f,useall) {
        if(useall) return(f)
        factor(f,levels=c(NA,levels(f)[-1]))
      }, flist, codes==2, SIMPLIFY=FALSE)

      if(sum(isfac) > 1) {
        iaf <- Reduce(':',fl)
#        iaf <- eval(as.call(lapply(c(':',contains[isfac]),as.name)), mf,environment(formula))
      } else {
        iaf <- fl[[1]] #eval(as.name(contains[isfac]), mf)
      }
      if(any(!isfac)) {
        attr(iaf,'x') <- as.matrix(Reduce('*',eval(as.call(lapply(c('list',contains[!isfac]),as.name)),
                                                   mf,environment(formula))))
      } else {
        attr(iaf,'x') <- numeric(0)
      }
      iaf
    })
    names(faclist) <- colnames(fact)
    faclist <- Filter(Negate(is.null), faclist)
    list(mat=t(mat),faclist=faclist)
  })
  names(spec) <- levels(df)[-1]
  attr(spec,'d') <- d
  spec
}
