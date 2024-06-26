% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_int.R
\name{calc_eqdist}
\alias{calc_eqdist}
\title{Equilibrium distribution from movement matrix}
\usage{
calc_eqdist(x, nr = dim(x)[2], start = rep(1/nr, nr), nit = 20)
}
\arguments{
\item{x}{Movement matrix, a square matrix with rows corresponding to origin (sum to 1), and columns corresponding to destination}

\item{nr}{Number of regions}

\item{start}{The initial distribution. Vector of length \code{nr}}

\item{nit}{Integer, the number of times the movement matrix will be applied}
}
\value{
Numeric vector of length \code{nr}
}
\description{
Applies the movement matrix several times in order to obtain the equilibrium spatial distribution of a movement matrix.
Not used in the model but useful for reporting.
}
