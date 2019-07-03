#' A package for estimating a mixed proportional hazard competing risk model with the NPMLE.
#'
#' The main method of the package is \code{\link{mphcrm}}. It has an interface
#' somewhat similar to \code{\link{lm}}.  There is an example of use in \code{\link{datagen}}, with
#' a generated dataset similiar to the ones in \cite{Gaure et al. (2007)}. For those who have
#' used the program used in that paper, a mixture of R, Fortran, C, and python,
#' this is an entirely new self-contained package, written from scratch with 12 years of experience.
#' Currently not all functionality from that behemoth has been implemented, but most of it.
#'
#' The package may use more than one cpu, the default is taken from \code{getOption("durmod.threads")}
#' which is initialized from the environment variable \env{DURMOD_THREADS}, \env{OMP_THREAD_LIMIT},
#' \env{OMP_NUM_THREADS} or \env{NUMBER_OF_PROCESSORS}, or parallel::detectCores() upon loading the package.
#' 
#' There is a vignette (\code{vignette("whatmph")}) with more details about \pkg{durmod}.
#' @references
#' Gaure, S., K. Røed and T. Zhang (2007) \cite{Time and causality: A Monte-Carlo Assessment
#'   of the timing-of-events approach}, Journal of Econometrics 141(2), 1159-1195.
#' \url{https://doi.org/10.1016/j.jeconom.2007.01.015}
#' @name durmod-package
#' @aliases durmod
#' @docType package
#' @useDynLib durmod, .registration=TRUE
#' @importFrom Rcpp evalCpp
#' @import 'stats'
#' @import 'utils'
#' @importFrom nloptr nloptr
#' @importFrom numDeriv grad hessian jacobian
#' @importFrom parallel detectCores
#' @importFrom data.table data.table
#' @importFrom mvtnorm rmvnorm
NULL

