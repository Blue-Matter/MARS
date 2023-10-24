#' @section Slots inherited from Dsurvey:
#' \describe{
#' \item{`ni`}{Integer, number of indices of abundance. Zero is possible.}
#' \item{`Iobs_ymi`}{Observed indices}
#' \item{`Isd_ymi`}{Lognormal standard deviation of the observed indices}
#' \item{`unit_i`}{Character vector, units of the index. Set to `"B"` to use stock weight at age (default) or `"N"` for abundance (numbers).}
#' \item{`IAAobs_ymai`}{Survey age composition}
#' \item{`IALobs_ymli`}{Survey length composition}
#' \item{`icomp_like`}{Character, likelihood for the composition data. See [like_comp()] for options}
#' \item{`IAAN_ymi`}{Sample size of the index age composition by season if using the multinomial or Dirichlet-multinomial likelihoods}
#' \item{`IALN_ymi`}{Sample size of the index length composition by season if using the multinomial or Dirichlet-multinomial likelihoods}
#' \item{`IAAtheta_i`}{Index age composition dispersion parameter if using the Dirichlet-multinomial likelihood}
#' \item{`IALtheta_i`}{Index length composition dispersion parameter if using the Dirichlet-multinomial likelihood}
#' \item{`samp_irs`}{Boolean array that specifies the regions and stocks sampled by the index. `samp[i, r, s]` indicates whether index `i` operates in region `r` and catches stock `s`.}
#' \item{`sel_i`}{Character matrix of functional forms for selectivity. See `"type"` argument in `[conv_selpar()]` for options.}
#' \item{`delta_i`}{The elapsed fraction of time in the seasonal time step when the index samples the population.}
#' }

