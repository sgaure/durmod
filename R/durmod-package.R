#' A package for estimating a mixed proportional hazard competing risk model.
#'
#' The main method of the package is \code{\link{mphcrm}}. It has an interface
#' somewhat similar to \code{\link{lm}}.
#'
#' The package may use more than one cpu, the default is taken from \code{getOption("durmod.threads")}
#' which is initialized from the environment variable \env{DURMOD_THREADS}, \env{OMP_THREAD_LIMIT},
#' \env{OMP_NUM_THREADS} or \env{NUMBER_OF_PROCESSORS}, or parallel::detectCores() upon loading the package.
#'
#' @name durmod-package
#' @aliases durmod
#' @docType package
#' @useDynLib durmod, .registration=TRUE
#' @importFrom Rcpp evalCpp
#' @import 'stats'
#' @import 'utils'
#' @importFrom Formula as.Formula
#' @importFrom nloptr nloptr
#' @importFrom numDeriv grad hessian jacobian
#' @importFrom parallel detectCores
#' @importFrom data.table data.table
#' @importFrom mvtnorm rmvnorm
NULL

