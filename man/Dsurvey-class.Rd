% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/define-s4-class.R
\docType{class}
\name{Dsurvey-class}
\alias{Dsurvey-class}
\title{Dsurvey S4 object}
\description{
S4 class that organizes the various data inputs for the MARS model. \code{MARSdata} simply inherits the slots from 6 component classes:
\code{Dmodel}, \code{Dstock}, \code{Dfishery}, \code{Dsurvey} \code{DCKMR}, and \code{Dtag}, where the \code{D}- prefix denotes an object for model data.
}
\details{
For convenience, most arrays and matrices have the associated dimensions in the variable name. For example, \code{Cobs_ymfr} represents
observed catch with the dimension following the underscore, following this template:

\tabular{ll}{
\code{y} \tab Year \cr
\code{m} \tab Season \cr
\code{a} \tab Age \cr
\code{r} \tab Region \cr
\code{f} \tab Fishery \cr
\code{i} \tab Index \cr
\code{s} \tab Stock
}
}
\section{Slots inherited from Dsurvey}{

\describe{
\item{\code{ni}}{Integer, number of indices of abundance. Zero is possible.}
\item{\code{Iobs_ymi}}{Observed indices}
\item{\code{Isd_ymi}}{Lognormal standard deviation of the observed indices}
\item{\code{unit_i}}{Character vector, units of the index. Set to \code{"B"} to use stock weight at age (default) or \code{"N"} for abundance (numbers).}
\item{\code{IAAobs_ymai}}{Survey age composition}
\item{\code{IALobs_ymli}}{Survey length composition}
\item{\code{icomp_like}}{Character, likelihood for the composition data. See \code{\link[=like_comp]{like_comp()}} for options}
\item{\code{IAAN_ymi}}{Sample size of the index age composition by season if using the multinomial or Dirichlet-multinomial likelihoods}
\item{\code{IALN_ymi}}{Sample size of the index length composition by season if using the multinomial or Dirichlet-multinomial likelihoods}
\item{\code{IAAtheta_i}}{Index age composition dispersion parameter if using the Dirichlet-multinomial likelihood}
\item{\code{IALtheta_i}}{Index length composition dispersion parameter if using the Dirichlet-multinomial likelihood}
\item{\code{samp_irs}}{Boolean array that specifies the regions and stocks sampled by the index. \code{samp[i, r, s]} indicates whether index \code{i} operates in region \code{r} and catches stock \code{s}.}
\item{\code{sel_i}}{Character matrix of functional forms for selectivity. See \code{"type"} argument in \code{\link[=conv_selpar]{conv_selpar()}} for options.}
\item{\code{delta_i}}{The elapsed fraction of time in the seasonal time step when the index samples the population.}
}
}

\seealso{
\link{MARSdata-class} \code{\link[=check_data]{check_data()}} \link{Dmodel-class} \link{Dstock-class} \link{Dfishery-class} \link{Dsurvey-class} \link{DCKMR-class} \link{Dtag-class}
}
