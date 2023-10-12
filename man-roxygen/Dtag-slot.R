#' @section Slots inherited from Dtag:
#' \describe{
#' \item{`tag_ymrr`}{Array. Number of tags that move between regions. Used to inform movement matrices between time steps.}
#' \item{`tag_ymr`}{Array. Number of tags among regions. Used to inform stock distribution within time step.}
#' \item{`tag_like`}{Character. Likelihood for the tagging data, either the vector of proportions by region of origin for `tag_ymrr`, or
#' by region of stock presence for `tag_ymr`. See `type` argument of [like_comp()] for options}
#' }
