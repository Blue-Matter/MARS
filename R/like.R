
#' Likelihood for length or age composition vectors
#'
#' Returns the log-likelihood for composition data, with various statistical distributions supported.
#'
#' @param obs A vector of observed values
#' @param pred A vector of predicted values. Same length as `obs`
#' @param N Numeric, the sample size corresponding to `obs`
#' @param type Character for the desired distribution
#' @param ... Other arguments depending on `type`. See below.
#' @return Numeric representing the log-likelihood.
#'
#' @details
#' For `type = "dirmult1"` or `"dirmult2"`, provide `p` which represents the linear or saturating parameter, respectively.
#'
#' @references
#' Thorson et al. 2017. Model-based estimates of effective sample size in stock assessment models using the
#' Dirichlet-multinomial distribution. Fish. Res. 192:84-93. \doi{10.1016/j.fishres.2016.06.005}
#' @examples
#' M <- 0.1
#' age <- seq(1:10)
#' pred <- exp(-M * age)
#' obs <- pred * rlnorm(10, sd = 0.05)
#' like_comp(obs, pred, N = 10, type = "multinomial", N = 10)
#' like_comp(obs, pred, N = 100, type = "multinomial", N = 100)
#' like_comp(obs, pred, N = 10, type = "dirmult1", p = 1)
#' like_comp(obs, pred, N = 10, type = "dirmult1", p = 20)
#' @export
like_comp <- function(obs, pred, N = sum(obs), type = c("multinomial", "dirmult1", "dirmult2"), ...) {
  stopifnot(length(obs) == length(pred))
  type <- match.arg(type)
  dots <- list(...)

  if (all(is.na(obs)) || !sum(obs)) {

    v <- 0

  } else if (type == "multinomial") {

    v <- dmultinom(x = N * obs/sum(obs), prob = pred, log = TRUE)

  } else if (type == "dirmult1") {

    alpha <- dots$p * N * pred/sum(pred)
    v <- ddirmnom(obs, size = N, alpha = alpha, log = TRUE)

  } else if (type == "dirmult2") {

    alpha <- dots$p * pred/sum(pred)
    v <- ddirmnom(obs, size = N, alpha = alpha, log = TRUE)

  }
  return(v)
}

ddirmnom <- function(x, size, alpha, log = FALSE) {
  x <- size * x/sum(x)
  val <- lgamma(sum(alpha)) + lgamma(sum(x) + 1) - lgamma(sum(alpha) + sum(x))
  val2 <- lgamma(x + alpha) - lgamma(alpha) - lgamma(x + 1)

  log_res <- val + sum(val2)

  if (log) log_res else exp(log_res)
}

