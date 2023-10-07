#' @section Slots inherited from Dmodel:
#' \describe{
#' \item{`ny`}{Integer, number of years}
#' \item{`nm`}{Integer, number of seasons}
#' \item{`na`}{Integer, number of ages}
#' \item{`nl`}{Integer, number of length bins. Set to zero if lengths are not modeled.}
#' \item{`nr`}{Integer, number of spatial regions}
#' \item{`ns`}{Integer, number of stocks}
#' \item{`lbin`}{Vector of lower boundary of length bins. Length `nl + 1`}
#' \item{`lmid`}{Vector of midpoint of length bins. Length `nl`}
#' \item{`Fmax`}{Numeric, maximum allowable instantaneous fishing mortality rate (units of per season). Defaults to 3.}
#' \item{`nitF`}{Integer, number of iterations to solve Baranov catch equation from observed catch. Defaults to 5.}
#' \item{`dist_type`}{Character. Whether to estimate seasonal stock distribution with a movement matrix `"mov"` or re-distributing total abundance by a distribution vector `"dist"`}
#' \item{`y_phi`}{Integer, the year from which to obtain values of natural mortality and fecundity for the unfished stock-recruit replacement line (`phi`). Relevant if natural mortality or fecundity are time-varying. Defaults to 1.}
#' }
