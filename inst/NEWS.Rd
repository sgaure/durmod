\name{NEWS}
\title{durmod news}
\section{Changes in version 1.1}{
  \itemize{
    \item Remove covariates and factor levels which are constant for a
    transition. They will act as intercepts, we don't want that. If you
    really want such a thing, create it with some neglible variation,
    e.g. as \code{rnorm(N,1,10*.Machine$double.eps)}.
    \item Added generics for \code{vcov()}, \code{coef()}.
  }
}
