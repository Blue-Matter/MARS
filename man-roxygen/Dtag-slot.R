#' @section Slots inherited from Dtag:
#' \describe{
#' \item{`tag_ymarrs`}{Array. Number of tags that move between regions. Informs movement matrices of stocks between time steps.}
#' \item{`tag_ymars`}{Array. Number of tags distributed among regions. Informs stock distribution (within time step).}
#' \item{`tag_yy`}{Boolean matrix that aggregates years for the tag data. Only used for the tag movement array `tag_ymarrs`.}
#' \item{`tag_aa`}{Boolean matrix that aggregates ages for the tag data.}
#' \item{`tag_like`}{Character. Likelihood for the tagging data, either the vector of proportions by region of origin for `tag_ymarrs`, or
#' by region of stock distribution for `tag_ymars`. See `type` argument of [like_comp()] for options}
#' \item{`tagN_ymars`}{Array. Sample size of the tag movement vectors if using the multinomial or Dirichlet-multinomial likelihoods.}
#' \item{`tagN_ymas`}{Array. Sample size of the tag distribution vectors if using the multinomial or Dirichlet-multinomial likelihoods.}
#' \item{`tagtheta_s`}{Array. Tag dispersion parameter (by stock) if using the Dirichlet-multinomial likelihoods. Default set to 1.}
#' \item{`tagstdev_s`}{Array. Tag standard deviation (by stock) if using the lognormal likelihood. Default set to 0.1.}
#' }
