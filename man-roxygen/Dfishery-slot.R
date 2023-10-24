#' @section Slots inherited from Dfishery:
#' \describe{
#' \item{`nf`}{Integer, number of fleets}
#' \item{`Cobs_ymfr`}{Total fishery catch}
#' \item{`fwt_yamfs`}{Fishery weight at age. Set to 1 when fleet catch is in units of abundance. Set to stock weight at age by default.}
#' \item{`CAAobs_ymafr`}{Fishery catch at age composition}
#' \item{`CALobs_ymlfr`}{Fishery catch at length composition}
#' \item{`fcomp_like`}{Character, likelihood for the fishery composition data. See `type` argument of [like_comp()] for options}
#' \item{`CAAN_ymfr`}{Sample size of the catch at age vector by season if using the multinomial or Dirichlet-multinomial likelihoods}
#' \item{`CALN_ymfr`}{Sample size of the catch at length vector by season if using the multinomial or Dirichlet-multinomial likelihoods}
#' \item{`CAAtheta_f`}{Catch at age dispersion parameter if using the Dirichlet-multinomial likelihood. Default set to 1.}
#' \item{`CALtheta_f`}{Catch at length dispersion parameter if using the Dirichlet-multinomial likelihood. Default set to 1.}
#' \item{`sel_block_yf`}{Index of dummy fleets to model time blocks of selectivity}
#' \item{`sel_f`}{Character vector of the functional form for selectivity. Choose between: `"logistic_length", "dome_length", "logistic_age", "dome_age", "SB", "B"`}
#' \item{`Cinit_mfr`}{Equilibrium seasonal catch prior to the first year. One way to initialize the abundance at the start of the first year
#' in the model. Default of zero.}
#' \item{`SC_ymafrs`}{Stock composition data.}
#' \item{`SC_aa`}{Boolean matrix that aggregates age classes for the stock composition data. See example.}
#' \item{`SC_ff`}{Boolean matrix that aggregates fleets for the stock composition data. See example.}
#' \item{`SC_like`}{Character, likelihood for the stock composition data. See `type` argument of [like_comp()] for options}
#' \item{`SCN_ymafr`}{Sample size of the stock composition vector if using the multinomial or Dirichlet-multinomial likelihoods}
#' \item{`SCtheta_f`}{Stock composition dispersion parameter if using the Dirichlet-multinomial likelihood. Default set to 1.}
#' \item{`SCstdev_f`}{Stock composition standard deviation if using the lognormal likelihood. Default set to 0.1.}
#' }
#' @examples
#' # Aggregate stock composition for ages 1-4 and 5-10 across all fleets
#' na <- 10
#' na_SC <- 2
#' SC_aa <- matrix(0, na_SC, na) # Assumes dim(SC_ymafrs)[3] = na_SC
#' SC_aa[1, 1:4] <- SC_aa[2, 5:10] <- 1
#'
#' nf <- 3
#' nf_SC <- 1
#' SC_ff <- matrix(1, nf_SC, nf) # Assumes dim(SC_ymafrs)[4] <- nf_SC
#'
