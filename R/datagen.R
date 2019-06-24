#' Generate example data
#'
#' @description
#' Generate a data table with example data
#' @details
#' The dataset simulates a labour market program. People entering the dataset are without a job,
#' there are two transitions, either to \code{"exit"}, in which case the person leaves the dataset,
#' or \code{"treatment"}. Those who have received treatment stays in the dataset, but may not
#' receive a new treatment. This process continues until either the person exits, or the
#' \code{"censor"} time is reached. There are two covariates, \code{"x1"} and \code{"x2"}, each with
#' their own effect on the hazards for the transitions. In addition there is some hidden
#' individual heterogeneity. It is this hidden heterogeneity which is modeled by the mixed proportional
#' hazard distribution.
#' @param N integer.
#' The number of individuals in the dataset
#' @param censor numeric. The total observation period.
#' @note
#' The example illustrates how data(durdata) was generated.
#' @examples
#' data.table::setDTthreads(1)
#' dataset <- datagen(5000,80)
#' print(dataset)
#' risksets <- list(untreated=c(1,2), treated=c(1))
#' # just two iterations to save time
#' Fit <- mphcrm(d ~ x1+x2|alpha, data=dataset, id=id, durvar=duration,
#'           state=alpha+1,risksets=risksets,
#'           control=mphcrm.control(threads=1,iters=2))
#' best <- Fit[[1]]
#' print(best)
#' summary(best)
#' @export
datagen <- function(N,censor=80) {
  x1 <- rnorm(N)
  x2 <- rnorm(N)
  id <- seq_len(N)
# two transitions, to exit(1) and to program(2)
# exit is absorbing, program is not, but only exit is allowed afterwards
# if on program one can only exit
  means <- c(-2.5,-3)
  cov2 <- matrix(c(1,0.5,0.5,1),2)
  persons <- data.table(id,x1,x2)
# draw correlated unobserved characteristics
  `:=` <- .N <- ve <- vp <- NULL #avoid check NOTEs
  persons[, c('ve','vp') := {
    vv <- mvtnorm::rmvnorm(.N, mean=means, sigma=cov2)
    list(vv[,1],vv[,2])
  }]
  # create spells
  spells <- persons[,{
    genspell(x1,x2,ve,vp,censor)
  }, by=id]
  spells$d <- factor(spells$d,levels=0:2)
  levels(spells$d) <- c('none','exit','treatment')
  spells
}
