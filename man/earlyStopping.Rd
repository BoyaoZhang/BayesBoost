% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/earlyStopping.R
\name{earlyStopping}
\alias{earlyStopping}
\title{cAIC-based stopping criteria}
\usage{
earlyStopping(cAIC, pat = 3, min_iters = 20)
}
\arguments{
\item{cAIC}{a serices of cAICs.}

\item{pat}{patient paramter, default is 3.}

\item{min_iters}{the least number of boosting iterations the algorithm should perform,
default is 20. The actual minimum is combined with pat parameter, i.e., "min_iters + pat".}
}
\description{
Determine the stopping iteration of cAICs based on the existence of
  lower values in the following "patience" iterations.
}
