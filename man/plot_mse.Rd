% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_summary.R
\name{plot_mse}
\alias{plot_mse}
\title{Plot mean squared error (MSE)}
\usage{
plot_mse(dt, HR = 1, driftHR = 1, pred = "none")
}
\arguments{
\item{dt}{a \code{data.frame} containing summary statistics for the posterior samples from each simulation generated with \code{get_summary}.}

\item{HR}{pre-specified HR between treatment and control arm in the internal trial for which the MSE should be plotted. Must be within \code{unique(dt$HR)}.}

\item{driftHR}{pre-specified HR between external control arm and internal control arm for which the MSE should be plotted. Must be within \code{unique(dt$driftHR)}.}

\item{pred}{predictors to use when fitting exponential distribution in MCMC for which the MSE should be plotted. Must be within \code{unique(dt$pred)}.}
}
\value{
a bar plot of class \code{ggplot} containing the MSE for each prior distribution.
}
\description{
Plot mean squared error (MSE) for each prior distribution according to selected simulation parameters
}
\keyword{method}
