#' @section Slots inherited from DCKMR:
#' \describe{
#' \item{`POP_s`}{A list by stock of data frames for parent-offspring pairs. Each row in the data frame corresponds to a "sampling unit" defined by the columns:
#' \tabular{ll}{
#' `a` \tab Capture year of parent \cr
#' `t` \tab Age at capture of parent \cr
#' `y` \tab Birth year of offspring \cr
#' `n` \tab Number of pairwise comparisons \cr
#' `m` \tab Number of POPs \cr
#' }
#' }
#' \item{`HSP_s`}{A list by stock of data frames for half-sibling pairs. Each row in the data frame corresponds to a "sampling unit" defined by the columns:
#' \tabular{ll}{
#' `yi` \tab Birth year of older sibling \cr
#' `yj` \tab Birth year of younger sibling \cr
#' `n` \tab Number of pairwise comparisons \cr
#' `m` \tab Number of HSPs
#' }
#' }
#' \item{`CKMR_like`}{Character, likelihood for the POP and HSP sampling units. See `type` argument in [like_CKMR()] for options.}
#' }
