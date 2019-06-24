#' Duration data
#'
#' The dataset simulates a labour market program. People entering the dataset are without a job,
#' there are two transitions, either to \code{"exit"}, in which case the person leaves the dataset,
#' or \code{"treatment"}. Those who have received treatment stays in the dataset, but may not
#' receive a new treatment. This process continues until either the person exits, or the
#' \code{"censor"} time is reached. There are two covariates, \code{"x1"} and \code{"x2"}, each with
#' their own effect on the hazards for the transitions. In addition there is some hidden
#' individual heterogeneity. It is this hidden heterogeneity which is modeled by the mixed proportional
#' hazard distribution.
#'
#' The dataset has already been fitted in the \code{fit} object.
#'
#' @docType data
#' @aliases fit
#' @usage data(durdata)
#'
#' @format A data.frame
#'
#' @keywords datasets
#'
#' @examples
#' data(durdata)
#' print(durdata)
#' print(fit)
#' summary(fit[[1]])
"durdata"
