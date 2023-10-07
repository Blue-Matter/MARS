#' @section Slots inherited from Dfishery:
#' \describe{
#' \item{`nf`}{Integer, number of fleets}
#' \item{`Cobs_ymfr`}{Total fishery catch}
#' \item{`fwt_yamfs`}{Fishery weight at age. Set to 1 when fleet catch is in units of abundance. Set to stock weight at age by default.}
#' \item{`CAAobs_ymafr`}{Fishery catch at age composition}
#' \item{`CALobs_ymlfr`}{Fishery catch at length composition}
#' \item{`fcomp_like`}{Character, likelihood for the fishery composition data. See `type` argument in \link{like_comp} for options}
#' \item{`CAAN_ymfr`}{Sample size of the catch at age vector by season if using the multinomial or Dirichlet-multinomial likelihoods}
#' \item{`CALN_ymfr`}{Sample size of the catch at length vector by season if using the multinomial or Dirichlet-multinomial likelihoods}
#' \item{`CAAtheta_f`}{Catch at age dispersion parameter if using the Dirichlet-multinomial likelihood. Default set to 1.}
#' \item{`CALtheta_f`}{Catch at length dispersion parameter if using the Dirichlet-multinomial likelihood. Default set to 1.}
#' \item{`sel_block_yf`}{Index of dummy fleets to model time blocks of selectivity}
#' \item{`sel_f`}{Character matrix of functional form for selectivity}
#' \item{`SC_ymafrs`}{Stock composition data.}
#' \item{`SC_aa`}{Boolean matrix that groups together age classes for the stock composition data}
#' \item{`SC_ff`}{Boolean matrix that groups together fleets for the stock composition data}
#' \item{`SC_like`}{Character, likelihood for the stock composition data. See `type` argument of \link{like_comp} for options}
#' \item{`SCN_ymafr`}{Sample size of the stock composition vector if using the multinomial or Dirichlet-multinomial likelihoods}
#' \item{`SCtheta_f`}{Stock composition dispersion parameter if using the Dirichlet-multinomial likelihood. Default set to 1.}
#' \item{`SCstdev_f`}{Stock composition standard deviation if using the lognormal likelihood. Default set to 0.1.}
#' }
