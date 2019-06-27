#' Estimate a mixed proportional hazard model
#'
#' @description
#' \code{mphcrm} implements estimation of a mixed proportional hazard competing risk model.
#' The baseline hazard is of the form \eqn{exp(X \beta)} where \eqn{X} is a
#' matrix of covariates, and \eqn{\beta} is a vector of parameters to estimate.
#' In addition there is an intercept term \eqn{\mu}, i.e. the hazard is \eqn{exp(X \beta + \mu)}.
#' There are several transitions to be made, and a set of \eqn{X}, \eqn{\beta}, and \eqn{\mu} for
#' each possible transition.
#'
#' Each individual may have several observations, with either a transition at the end of the
#' observation, or not a transition. It is a competing risk, there can be more than one possible
#' transition for an observation, but only one is taken at the end of the period.
#' 
#' For each individual \eqn{i} there is a log likelihood as a function of \eqn{\mu}, called \eqn{M_i(\mu)}.
#'
#' The mixture part is that the \eqn{\mu}'s are stochastic. I.e. 
#' we have probabilities \eqn{p_j}, and a vector of \eqn{\mu_j} of masspoints (one for each transition),
#' for each such \eqn{j}.
#'
#' So the full likelihood for an individual is \eqn{L_i = \sum_j p_j M_i(\mu_j)}.
#'
#' The \code{mphcrm()} function maximizes the likelihood \eqn{\sum_i L_i} over
#' \eqn{p_j}, \eqn{mu_j}, and \eqn{\beta}.
#'
#' In addition to the covariates specified by a formula, a variable which records the duration of each
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
#' In a formula like \code{d ~ x1 + x2 + ID(id) + D(dur) + C(job,alpha) + S(state)},
#' the \code{d} is the transition which is taken, coded as an integer where \code{0} means
#' no transition, and otherwise \code{d} is the number of the transition which is taken.
#' \code{d} can also be a factor, in which case the level which is no transition must be named
#' \code{"0"} or \code{"none"}. If \code{d} is an integer, the levels for transitions will be named
#' \code{"t1"}, \code{"t2"}, .... And \code{"none"}.
#'
#' The \code{x1+x2} part is like in \code{\link{lm}}, i.e. ordinary covariates or factors.
#'
#' The \code{D()} specifies the covariate which holds the duration of each observation. The
#' transition in \code{d} is assumed to be taken at the end of this period.
#'
#' The \code{ID()} part specifies the covariate which holds the individual identification.
#'
#' The \code{S()} specifies the covariate which holds an index into the \code{risksets} list.
#'
#' These three special symbols are replaced with \code{I()}, so it is possible to have calculations
#' inside them.
#' 
#' If the covariates differ among the transitions, one can specify covariates conditional on
#' the transition taken. If e.g. the covariates \code{alpha} and \code{x3} should only explain transition to
#' job, specify \code{C(job, alpha+x3)}. This comes in addition to the ordinary covariates.
#' Then name \code{job} refers to a level in \code{d}, the transition taken. 
#' @param data
#' A data frame which contains the covariates.
#' @param timing
#' character. The timing in the duration model. Can be one of
#' \itemize{
#' \item "exact" The timing is exact, the transition occured at the end of the observation interval.
#' \item "interval" The transition occured some time during the observation interval. This model
#'   can be notoriously hard to estimate due to unfavourable numerics. It could be worthwhile to
#'   lower the control parameter \code{ll.improve} to 0.0001 or something, to ensure it does not
#'   give up to early.
#' \item "none" There is no timing, the transition occured, or not. A logit model is used.
#' }
#' @param risksets
#' A list of character vectors. Each vector is a list of transitions, i.e. which risks are present for
#' the observation. The elements of the vectors must be levels of the covariate which is the
#' left hand side of the \code{formula}.
#' If the state variable in the formula is a factor, the \code{risksets} argument should be a named list, with
#' names matching the levels of \code{state}.
#' @param subset
#' For specifying a subset of the dataset, similar to \code{\link{lm}}.
#' @param na.action
#' For handling of missing cases, similar to \code{\link{lm}}.
#' @param control
#' List of control parameters for the estimation. See \code{\link{mphcrm.control}}.
#' @return
#' A list, one entry for each iteration. Ordered in reverse order. Ordinarily you will be
#' interested in the first entry.
#' @details
#' The estimation starts by estimating the null-model, i.e. all parameters set to 0, only one
#' intercept for each transition is estimated.
#'
#' Then it estimates the full model, still with one intercept in each transition.
#'
#' After the initial model has been estimated, it tries to add a masspoint to the mixing
#' distribution, then estimates the model with this new distribution.
#'
#' It continues to add masspoints in this way until eithr it can not improve the likelihood, or
#' the number of iterations as specified in control$iters are reached.
#'
#' The result of every iteration is returned in a list.
#'
#' If you interrupt \code{mphcrm} it will catch the interrupt and return with the
#' estimates it has found so far.
#' @note
#' The algorithm is not fully deterministic. New points are searched for randomly, there
#' is no canonical order in which they can be found. It can happen that a point is found
#' early which makes the rest of the estimation hard, so it terminates early. In particular
#' when using interval timing. One should
#' make a couple of runs to ensure they yield reasonably equal results.
#' @examples
#' data(durdata)
#' head(durdata)
#' risksets <- list(c('job','program'), c('job'))
#' Fit <- mphcrm(d ~ x1+x2 + C(job,alpha) + ID(id) + D(duration) + S(alpha+1), data=durdata, 
#'      risksets=risksets, control=mphcrm.control(threads=1,iters=2))
#' best <- Fit[[1]]
#' summary(best)
#' @seealso A description of the dataset is available in \code{\link{datagen}} and \code{\link{durdata}},
#' and in the vignette \code{vignette("mphcrmwhat")}
#' @export
mphcrm <- function(formula,data,risksets=NULL,
                   timing=c('exact','interval','none'),
                   subset, na.action, control=mphcrm.control()) {

  timing <- match.arg(timing)
  F <- formula

  # create model frame
  mf <- match.call(expand.dots=TRUE)
  m <- match(c('formula','data','subset','na.action'), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(model.frame)

  dataset <- mymodelmatrix(F,mf)

  dataset$timing <- timing
  id <- dataset$id

  if(length(unique(id)) != length(rle(as.integer(id))$values)) {
    stop('dataset must be sorted on ID')
  }

  # zero-based index of beginning of spells. padded with one after the last observation
  dataset$spellidx <- c(0L,which(diff(as.integer(id))!=0),length(id))
  dataset$nspells <- length(dataset$spellidx)-1L

  state <- dataset$state
  hasriskset <- !is.null(risksets)
  if(hasriskset && is.null(state)) warning("riskset is specified, but no state S()")
  if(is.null(state)) hasriskset <- FALSE

  if(hasriskset) {
    if(is.factor(state)) {
      m <- match(levels(state), names(risksets))
      if(anyNA(m)) stop(sprintf('level %s of state is not in the riskset names\n',levels(state)[is.na(m)]))
      # recode state to integer
      state <- m[state]
      dataset$state <- state
    } else {
      srange <- range(state)
      if(srange[1] != 1) stop('smallest state must be 1 (index into riskset)')
      if(srange[2] > length(risksets)) {
        stop(sprintf('max state is %d, but there are only %d risksets\n',srange[2],length(risksets)))
      }
    }
    # recode risksets from transition level names to integers
    tlevels <- dataset$tlevels
    risksets <- lapply(risksets, function(set) {
      ind <- match(set,tlevels)
      if(anyNA(ind)) stop(sprintf('Non-existent transition %s in risk set\n',set[is.na(ind)]))
      ind
    })

    # check if transitions are taken which are not in the riskset
    d <- dataset[['d']]

    badtrans <- mapply(function(dd,r) dd != 0 && !(dd %in% r), d, risksets[state])
    if(any(badtrans)) {
      n <- which(badtrans)[1]
      stop(sprintf("In observation %d(id=%d), a transition to %s is taken, but the riskset of the state(%d) does not allow it",
                   n, id[n], tlevels[d[n]], state[n]))
    }

    dataset$riskset <- risksets
  } else {
    dataset$state <- 0L
    dataset$riskset <- list()
  }

  pset <- makeparset(dataset,1)

  # reset timer
  if(is.null(control$callback)) control$callback <- function(...) {}
  if(!is.null(control$cluster)) {prepcluster(dataset,control); on.exit(cleancluster(control$cluster))}

  z <- pointiter(dataset,pset,control)

  structure(z,call=match.call())
}

#' Control parameters for mphcrm
#'
#' @description
#' Modify the default estimation parameters for
#'   \code{\link{mphcrm}}
#' @param ...
#' parameters to change \itemize{
#' \item threads integer. The number of threads to use. Defaults to \code{getOption('durmod.threads')}
#' \item iters integer. How many iterations should we maximally run. Defaults to 12.
#' \item ll.improve numeric. How much must the be log-likelihood improve from the last iteration before
#'   termination. Defaults to 0.001
#' \item newpoint.maxtime numeric. For how many seconds should a global search for a new point
#'   improving the likelihood be conducted before we continue with the best we have found. Defaults to
#'   120.
#' \item callback. A
#'   user-specified \code{function(fromwhere, opt, dataset, control,
#'   ...)} which is called after each optimization step.  It can be
#'   used to report what is happening, check whatever it wants, and
#'   optionally stop the estimation by calling stop(). In this case,
#'   \code{mphcrm()} will return with the most recently estimated set
#'   of parameters. See the help on \code{\link{mphcrm.callback}} for
#'   information on the argument.
#' \item trap.interrupt logical. Should interrupts be trapped so that \code{mphcrm} returns gracefully?
#' In this case the program will continue. Defaults to \code{interactive()}.
#' \item method character. The method used for optimization. Currently 'BFGS' is supported, but it
#' might work with 'L-BFGS-B' also, if memory is constrained.
#' }
#' @note
#' There are other parameters too, I may document them later. Instead of cluttering
#' the source code with constants and stuff required by various optimization routines, they
#' have been put in this control list. You should perhaps not change them.
#' @return
#' list of control parameters suitable for the \code{control}
#'   argument of \code{\link{mphcrm}}
#' @export
mphcrm.control <- function(...) {
  ctrl <- list(iters=15,threads=getOption('durmod.threads'),gradient=TRUE, fisher=TRUE, hessian=FALSE, 
               method='BFGS', gdiff=FALSE, minprob=1e-20, eqtol=1e-4, newprob=1e-4, jobname='mphcrm', 
               ll.improve=1e-3, e.improve=1e-3,
               trap.interrupt=interactive(),
               tspec='%T', newpoint.maxtime=120,
               callback=mphcrm.callback,
               cluster=NULL)
  args <- list(...)
  bad <- !(names(args) %in% names(ctrl))
  if(any(bad)) {message("unknown control parameter(s): ", paste(names(args)[bad], collapse=' '))}
  ctrl[names(args)] <- args
  ctrl
}

#' Default callback function for mphcrm
#'
#' The default callback function prints a line whenever estimation with a masspoint is
#' completed.
#' @param fromwhere
#' a string which identifies which step in the algorithm it is called from. \code{fromwhere=='full'} means
#' that it is a full estimation of all the parameters. There are also other codes, when adding a point,
#' when removing duplicate points. When some optimization is completed it is called with the
#' return status from \code{\link{optim}} (and in some occasions from \code{\link[nloptr]{nloptr}}.
#'
#' @param opt
#' Typically the result of a call to \code{\link{optim}}.
#' @param dataset
#' The dataset in a structured form.
#' @param control
#' The \code{control} argument given to \code{\link{mphcrm}}
#' @param ...
#' other arguments
#' @details
#' If you write your own callback function it will replace the default function, but you can
#' of course call the default callback from your own callback function, and in addition print your
#' own diagnostics, or save the intermediate \code{opt} in a file, or whatever.
#' 
#' @note
#' Beware that
#' \code{control} contains a reference to the callback function, which may contain a reference
#' to the top-level environment, which may contain the full dataset. So if you save \code{control}
#' you may end up saving the entire dataset.
#' @export
mphcrm.callback <- local({
  lastfull <- Sys.time()
  function(fromwhere, opt, dataset, control, ...) {
    now <- Sys.time()
    jobname <- control$jobname
    if(fromwhere == 'removepoints') {
      p <- a2p(opt$pargs)
      bad <- list(...)[['remove']]
      cat(jobname,  format(Sys.time(),control$tspec), 'remove probs', sprintf(' %.2e',p[bad]),'\n')
    } else {
      if(is.numeric(opt$convergence) && opt$convergence != 0) {
      if(fromwhere != 'newpoint' || !(opt$convergence %in% c(2))) 
        cat(jobname, format(now,control$tspec), fromwhere, 'convstatus:"',opt$convergence, opt$message,'"\n')
      }
    }

    if(fromwhere != 'full') return()
    tdiff <- timestr(difftime(now,lastfull,units='s'))
    assign('lastfull',now,environment(sys.function())) # 
    rc <- if(!is.null(opt$fisher)) rcond(opt$fisher) else NA
    grd <- if(!is.null(opt$gradient)) sqrt(sum(opt$gradient^2)) else NA
    p <- a2p(opt$par$pargs)
    cat(jobname,format(now,control$tspec), 
        sprintf('i:%d pts:%d L:%.4f g:%.3g mp:%.5g rcond:%.2g e:%.4f t:%s\n',
                control$mainiter, length(opt$par$pargs)+1, opt$value, grd, min(p), rc, -sum(p*log(p)),
                tdiff))
    #  if(rc < sqrt(.Machine$double.eps)) {message(sprintf('%s  ***bad condition: %.2e',jobname,rc)); return(FALSE)}
  }
})


pointiter <- function(dataset,pset,control) {
  # optimize null model first
  assign('lastfull',Sys.time(),envir=environment(mphcrm.callback))

  arg0 <- flatten(pset)
  arg0[] <- 0
  pset <- unflatten(arg0)
  LL0 <- function(arg,dataset,pset) {
    for(i in seq_along(pset$parset)) {
      pset$parset[[i]]$mu[] <- arg[i]
    }
    -mphloglik(dataset,pset,control=control)
  }
  opt0 <- optim(runif(length(pset$parset),-4,-1), LL0,method='BFGS',dataset=dataset,pset=pset)
#  message('zero model: '); print(opt0$par)
  for(i in seq_along(pset$parset)) {
    pset$parset[[i]]$mu[] <- opt0$par[i]
  }
  class(pset) <- 'mphcrm.pset'
  opt0$value <- -opt0$value
  opt0$par <- pset
  opt0$mainiter <- 0
  opt0$nobs <- dataset$nobs
  opt0$nspells <- dataset$nspells
  class(opt0) <- 'mphcrm.opt'

  control$callback('nullmodel',opt0,dataset,control)
  intr <- FALSE
  iopt <- opt <- list(nullmodel=opt0)
  tryCatch(
    {
      i <- 0; done <- FALSE
      while(!done) {
        i <- i+1
        control$mainiter <- i
        opt <- c(list(ml(dataset,pset,control)), opt)
        opt[[1]]$mainiter <- i
        names(opt)[1] <- sprintf('iter%d',i)
        control$callback('full',opt[[1]], dataset,control)
        ll.improve <- opt[[1]]$value - opt[[2]]$value
        e.improve <- abs(opt[[1]]$entropy - opt[[2]]$entropy)
        done <- (ll.improve < control$ll.improve && e.improve < control$e.improve)|| i >= control$iters
        if(!done) pset <- addpoint(dataset,opt[[1]]$par,opt[[1]]$value,control)
      }
    },
    error=function(e) {
      assign('intr',TRUE,environment(sys.function()))
      assign('iopt',structure(opt,status=conditionMessage(e)),environment(sys.function()))
      if(!control$trap.interrupt) stop(e)
      warning(e, 'Error occured, returning most recent estimate')
    },
    interrupt=function(e) {
      assign('intr',TRUE,environment(sys.function()))
      assign('iopt',structure(opt,status='interrupted'),environment(sys.function()))
      if(!control$trap.interrupt) stop('interrupt')
 #     try(control$callback('interrupt', iopt, dataset, control))
      warning(e, "returning most recent estimate ")
    },
    finally=if(intr) opt <- iopt)

  if(!intr) {
    badset <- badpoints(opt[[1]]$par,control)
    if(isTRUE(attr(badset, 'badremoved'))) {
      cat(sprintf('redo bad prob, ll=%.6f\n',opt[[1]]$value))
      opt[[1]]$par <- optfull(dataset,badset,control)
      cat(sprintf('new ll=%.6f\n',opt[[1]]$value))
    }
  }

# rescale the parameters
  scales <- unlist(lapply(dataset$data, function(dd) attr(dd,'recode')[[2]]))
  opt <- lapply(opt, function(op) {
    nm <- names(op$gradient)
    scalevec <- rep(1,length(nm))
    names(scalevec) <- nm
    scalevec[names(scales)] <- scales
    scalemat <- scalevec * rep(scalevec,each=length(scalevec))
    if(!is.null(op$gradient)) op$gradient <- op$gradient / scalevec
    if(!is.null(op$numgrad)) op$numgrad <- op$numgrad / scalevec
    if(!is.null(op$fisher)) op$fisher <- op$fisher / scalemat
    if(!is.null(op$hessian)) op$hessian <- op$hessian / scalemat
    for(tr in names(dataset$data)) {
      offset <- attr(dataset$data[[tr]],'recode')$offset
      scale <- attr(dataset$data[[tr]],'recode')$scale
      op$par$parset[[tr]]$pars[] <- op$par$parset[[tr]]$pars*scale
      muadj <- sum(op$par$parset[[tr]]$pars*offset)
      op$par$parset[[tr]]$mu <- op$par$parset[[tr]]$mu - muadj
    }
    op
  })

  structure(opt,class='mphcrm.list')
}

optprobs <- function(dataset,pset,control) {
  pfun <- function(a) {
    pset$pargs[] <- a
    -mphloglik(dataset,pset,control=control)
  }
#  message('optimize probs')
  aopt <- optim(pset$pargs,pfun,method='BFGS',control=list(trace=0,REPORT=1,factr=10))
#  message('probs:',sprintf(' %.7f',a2p(aopt$par)), ' value: ',aopt$value)
  pset$pargs[] <- aopt$par
  control$callback('prob',aopt,dataset,control)
  pset
}

optdist <- function(dataset,pset,control) {
  dfun <- function(a) {
    n <- length(pset$pargs)+1
    pset$pargs[] <- a[seq_len(n-1)]
    pos <- n
    for(i in seq_along(pset$parset)) {
      pset$parset[[i]]$mu[] <- a[pos:(pos+n-1)]
      pos <- pos+n
    }
    -mphloglik(dataset,pset,control=control)
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
  control$callback('dist',dopt,dataset,control)
  pset
}

newpoint <- function(dataset,pset,value,control) {
  gdiff <- control$gdiff
  pr <- a2p(pset$pargs)
  newprob <- if(gdiff) 0 else control$newprob
  newpr <- c( (1 - newprob)*pr, newprob)
#  newpr <- newpr/sum(newpr)
  #  newpr <- c(pr,0)
  newset <- makeparset(dataset,length(pr)+1,pset)
  newset$pargs[] = p2a(newpr)
  np <- length(newpr)
  ntrans <- length(pset$parset)
  fun <- function(mu,gdiff) {
    for(i in seq_along(mu)) {
      newset$parset[[i]]$mu[np] <- mu[i]
    }
    -mphloglik(dataset,newset,gdiff=gdiff,control=control)
  }
  args <- runif(length(newset$parset),-10,1)
  muopt <- nloptr::nloptr(args, fun, lb=rep(-10,ntrans), ub=rep(1,ntrans),gdiff=gdiff,
                          opts=list(algorithm='NLOPT_GN_ISRES',stopval=if(gdiff) 0 else -value-newprob,
                                    maxtime=control$newpoint.maxtime,maxeval=10000,population=10*length(args)))
  if(!gdiff && !(muopt$status %in% c(0,2))) {
    # that one failed, try gdiff instead
    newset$pargs[] <- p2a(c(pr,0))
    gdiff <- TRUE
    muopt <- nloptr::nloptr(args, fun, lb=rep(-10,ntrans), ub=rep(1,ntrans),gdiff=TRUE,
                            opts=list(algorithm='NLOPT_GN_ISRES',stopval=0,
                                      maxtime=control$newpoint.maxtime,maxeval=10000,population=10*length(args)))
  }
  muopt$objective <- -muopt$objective
  muopt$convergence <- muopt$status
  control$callback('newpoint',muopt,dataset,control)
  
  for(i in seq_along(newset$parset)) {
    newset$parset[[i]]$mu[np] <- muopt$solution[i]
  }
  if(gdiff) {
    #  newpr <- newpr/sum(newpr)
    newset$pargs[] = p2a(c((1-1e-5)*pr,1e-5))
  }
  optprobs(dataset,newset,control)
}

badpoints <- function(pset,control) {
  # are any masspoints equal?
  np <- length(pset$pargs)+1
  p <- a2p(pset$pargs)
  okpt <- p > control$minprob
  jobname <- control$jobname
  for(i in seq_len(np)) {
    mui <- sapply(pset$parset, function(pp) pp$mu[i])
    for(j in seq_len(i-1)) {
      if(!okpt[j]) next
      muj <- sapply(pset$parset, function(pp) pp$mu[j])
      if(max(abs(exp(muj) - exp(mui))) < control$eqtol) {
        mumat <- rbind(muj,mui)
        rownames(mumat) <- c(j,i)
        colnames(mumat) <- names(pset$parset)
        control$callback('equalpoints',pset,NULL,control,cbind(prob=c(p[j],p[i]),exp(mumat)))
#        cat(sprintf('%s %s points %d and %d are equal\n',jobname, format(Sys.time(),control$tspec),j,i),'\n')
#        print()
        okpt[i] <- FALSE
        p[j] <- p[j]+p[i]
      }
    }
  }

  if(!all(okpt)) {
    control$callback('removepoints',pset,NULL,control,remove=!okpt)
  }

  p <- p[okpt]
  p <- p/sum(p)
  pset$pargs <- p2a(p)
  for(i in seq_along(pset$parset)) {
    pset$parset[[i]]$mu <- pset$parset[[i]]$mu[okpt]
  }
  structure(pset,badremoved=!all(okpt))
}

# some functions for use in optim
#negative log likelihood
LL <- function(args,skel, dataset,ctrl) {
  pset <- unflatten(args,skel)
  -mphloglik(dataset,pset,control=ctrl)
}

# gradient of negative log likelihood
gLL <- function(args,skel, dataset,ctrl) {
  pset <- unflatten(args,skel)
  -attr(mphloglik(dataset,pset,dogradient=TRUE,control=ctrl),'gradient')
}

# fisher matrix
fLL <- function(args, skel, dataset, ctrl) {
  pset <- unflatten(args,skel)
  -attr(mphloglik(dataset,pset,dofisher=TRUE,control=ctrl),'fisher')
}

# numerical gradient, slow
dLL <- function(args, skel, dataset, ctrl) {
  numDeriv::grad(LL, args, dataset=dataset, skel=skel, ctrl=ctrl)
}

# hessian, numerically from gradient, slow
hLL <- function(args, skel, dataset, ctrl) {
  numDeriv::jacobian(gLL,args,dataset=dataset,skel=skel,ctrl=ctrl)
}


optfull <- function(dataset, pset, control) {
  args <- flatten(pset)
  opt <- optim(args,LL,gLL,method=control$method,
        control=list(trace=0,REPORT=10,maxit=20*length(args),lmm=60,
                     abstol=1e-4,reltol=1e-14),
        skel=attr(args,'skeleton'), ctrl=control,dataset=dataset)
  opt$par <- unflatten(opt$par)
  opt
}

addpoint <- function(dataset,pset,value,control) {
  # find a new point
  # optimize probabilities
  newset <- pset
  control$minprob <- 0
  for(i in 1:5) {
    newset <- newpoint(dataset,newset,value,control)
    newset <- optdist(dataset,newset,control)
    newset <- badpoints(newset,control)
    if(!isTRUE(attr(newset,'badremoved'))) break
    # optimize dist after removing point(s)
#    newset <- unflatten(optfull(dataset,newset,control)$par)
  } 
  newset
}

ml <- function(dataset,pset,control) {
  opt <- optfull(dataset,pset,control)
#  control$minprob <- 0
  opt$par <- badpoints(opt$par, control)
  sol <- opt$par

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

  skel <- attr(flatten(sol),'skeleton')



  opt$value <- -opt$value
  opt$par <- sol
  p <- a2p(sol$pargs)
  opt$entropy <- -sum(p*log(p))

  # create a covariance matrix, we use the inverse of the (negative) fisher matrix, it's fast to compute
  nm <- names(flatten(opt$par))

  if(isTRUE(control$gradient)) {
    opt$gradient <- gLL(flatten(opt$par),skel,dataset,control)
    names(opt$gradient) <- nm
  }
  if(isTRUE(control$fisher)) {
    opt$fisher <- fLL(flatten(opt$par),skel,dataset,control)
    dimnames(opt$fisher) <- list(row=nm,col=nm)
  }
  if(isTRUE(control$numgrad)) {
    opt$numgrad <- dLL(flatten(opt$par), skel, dataset, control)
    names(opt$numgrad) <- nm
  }
  if(isTRUE(control$hessian)) {
    opt$hessian <- hLL(flatten(opt$par),skel,dataset,control)
    dimnames(opt$hessian) <- list(nm,nm)
  }

  opt$nobs <- dataset$nobs
  opt$nspells <- dataset$nspells
  structure(opt,class='mphcrm.opt')
}

#' Convert probability parameters to probabilities
#'
#' @description
#' \code{\link{mphcrm}} parametrizes the probabilities that it optimizes.
#' For \eqn{n+1} probabilities there are \eqn{n} parameters \eqn{a_j}, such that
#' probability \eqn{P_i = \frac{a_i}{\sum_j \exp(a_j)}}, where we assume that \eqn{a_0 = 0}.
#'
#' @param a
#' a vector of parameters
#' 
#' @examples
#' # Draw 5 parameters
#' a <- rnorm(5)
#' a
#' # make 6 probabilities
#' p <- a2p(a)
#' p
#' # convert back
#' p2a(p)
#' 
#' @export
a2p <- function(a) {b <- c(0,a); p <- exp(b)/sum(exp(b)); ifelse(is.na(p),1,p)}
##' @export
#a2logp <- function(a) {b <- c(0,a); logp <- b - logsumofexp(b,0); ifelse(is.na(logp),0,logp)}

#' @rdname a2p
#' @param p
#' a vector of probabilities with sum(p) = 1
#' @export
p2a <- function(p) log(p/p[1])[-1]

makeparset <- function(dataset,npoints,oldset) {
  parset <- lapply(dataset$data, function(tt) {
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


mymodelmatrix <- function(formula,mf) {

  #### Handle the specials ####
  mt <- terms(formula, specials=c('ID','D', 'S', 'C'))
  spec <- attr(mt,'specials')
  # replace with I() in formula
  F <- formula
  splist <- lapply(unlist(spec), function(sp) attr(mt,'variables')[sp+1L])
  remove <- function(f, r) update(f, as.formula(substitute(. ~ . - R, list(R=r[[1]]))))
  pureF <- Reduce(remove, splist, F)

  Slist <- splist[names(splist) %in% c('ID','D','S')]
  Ilist <- lapply(Slist, function(ss) substitute(I(R)(), list(R=ss[[1]][[2]])))
  addI <- function(f,r) update(f, as.formula(substitute(. ~ . + I(R), list(R=r[[1]][[2]]))))
  IF <- Reduce(addI,Slist,pureF)

  # then the C-parts, it needs special treatment
  Clist <- splist[!(names(splist) %in% c('ID','D','S'))]
  names(Clist) <- sapply(Clist, function(cc) eval(cc[[1L]], list(C=function(a,b) as.character(substitute(a)))))
  cvars <- lapply(Clist, function(cc) eval(cc[[1]], list(C=function(a,b) substitute(b))))
  addc <- function(f,r) update(f, as.formula(substitute(. ~ . + R, list(R=r))))
  IF <- Reduce(addc,cvars,IF)
  # Now, IF is a suitable formula with all specials removed.
  mf[[2L]] <- IF
  mf <- eval.parent(mf)

  # Pick up the special covariates ID, D, and S
  id <- as.integer(mf[[as.character(Ilist$ID)]])

  if(length(Ilist$D))
    duration <- as.numeric(mf[[as.character(Ilist$D)]])
  else
    duration <- rep(1,nrow(mf))

  state <- NULL
  if(length(Ilist$S)) {
    state <- mf[[as.character(Ilist$S)]]
    if(!is.factor(state)) state <- as.integer(state)
  }
  # and the repsonse, convert to factor with appropriate levels
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
  # make a zero-based d for use in C
  d <- as.integer(df)-ztrans

  # analyze the formula to figure out which
  # covariates explain which transitions
  # split off factors

  # first check that the conditional covariates specifies existing transitions.
  tlevels <- levels(df)
  if(ztrans) tlevels <- tlevels[-1]
  if(length(Clist)) {
    trnames <- names(Clist)
    enm <- match(trnames,tlevels)
    if(anyNA(enm)) stop(sprintf('transition to %s specified in conditional covariate, but no such transition exists.\n',
                                trnames[is.na(enm)]))
  }  

  transitions <- nlevels(df)-ztrans
  cls <- attr(terms(mf),'dataClasses')


  # for each transition, make a model matrix for the numeric covariates,
  # and a list of factors
  data <- lapply(seq_len(transitions), function(t) {
    
    # find the formula for
    # this transition. Add in all the conditional covariates for this transition
    ff <- pureF
    tnam <- levels(df)[t+ztrans]
    if(tnam %in% names(cvars)) {
      # awkward, cvars may have duplicate names if C(foo, a+b) + C(foo,d+c)
      for(i in (which(names(cvars) %in% tnam)))
        ff <- update(ff, as.formula(bquote(. ~ . + .(cvars[[i]]))))
    } 

    mt <- terms(ff,keep.order=TRUE)
    fact <- attr(mt,'factors')
    # now, filter out those terms which are neither factors, nor interactions with factors
    keep <- colSums(fact) > 0
    fact <- fact[,keep,drop=FALSE]
    attr(mt,'term.labels') <- attr(mt,'term.labels')[keep]

    keep <- unlist(lapply(colnames(fact), function(lab) {
      contains <- rownames(fact)[fact[,lab] > 0]
      if(!any(cls[contains]=='factor')) return(lab)
      return(NULL)
    }))

    attr(mt,'factors') <- fact[,keep,drop=FALSE]
    attr(mt,'intercept') <- 0
    # create model matrix from the constructed terms
    mat <- t(model.matrix(mt,mf))

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
      # remove a reference level if contrasts
      fl <- mapply(function(f,useall) {
        if(useall) return(f)
        factor(f,levels=c(NA,levels(f)[-1]))
      }, flist, codes==2, SIMPLIFY=FALSE)

      iaf <- Reduce(`:`, fl)
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

    # force covariates to be between 0 and 1.
    # i.e. subtract minimum, divide by span
    # the estimated betas must also be divided by span
    # beta <- estbeta / scale
    # the gradient and Fisher matrix must also be adjusted
    # the intercept estimates, the mus, must be adjusted
    # mu <- estmu - sum(rang[1,]*beta)
    # all this happens in the function ml()

    offset <- apply(mat,1,mean)
#    scale <- 1/apply(mat,1,sd)
    scale <- 1/apply(mat,1,function(x) max(abs(x-mean(x))))
    scalemat <- (mat-offset)*scale
    structure(list(mat=scalemat,faclist=faclist), recode=list(offset=offset,scale=scale))
  })
  names(data) <- tlevels
  dataset <- list(data=data, d=d, nobs=nrow(mf), tlevels=tlevels,duration=duration,id=id,state=state)
  dataset
}
