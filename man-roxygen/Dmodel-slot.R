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
#' \item{`y_phi`}{Integer, the year from which to obtain values of natural mortality and fecundity for the unfished stock-recruit replacement line (`phi`). Relevant if natural mortality or fecundity are time-varying. Defaults to 1.}
#' \item{`scale_s`}{Vector, length `ns`. Multiplicative scaling factor that informs relative stock size to aid parameter estimation. Larger values implies larger stocks. Default set to 1. See [make_parameters()].}
#' \item{`nyinit`}{Integer, number of years of spool-up to calculate equilibrium unfished and starting conditions for the population model to account for seasonal and spatial dynamics. The numerical spool-up is not needed when both `nm = 1` and `nr = 1`, i.e., `nyinit = 1`. Otherwise, set to `1.5 * na` by default.}
#' \item{`condition`}{Character, either to specify the model estimates fishing mortality as a parameter (`"F"`, default) or equal to the catch (`"catch"`).}
#' \item{`nitF`}{Integer, number of iterations to solve Baranov catch equation from observed catch if `condition = "catch"`. Defaults to 5.}
#' \item{`y_Fmult_f`}{Integer vector by fleet, the year in which to directly estimate F. Choose a year/season/region combination when the catch is average relative to the time series. Only used if `condition = "F"`.}
#' \item{`m_Fmult_f`}{Integer vector by fleet, the season in which to directly estimate F. Choose a year/season/region combination when the catch is average relative to the time series. Only used if `condition = "F"`.}
#' \item{`r_Fmult_f`}{Integer vector by fleet, the region in which to directly estimate F. Choose a year/season/region combination when the catch is average relative to the time series. Only used if `condition = "F"`.}
#' }
