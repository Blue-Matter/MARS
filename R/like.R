
#' Likelihood for composition vectors
#'
#' Returns the log-likelihood for composition data, e.g., length, age, or stock composition,
#' with various statistical distributions supported.
#'
#' @param obs A vector of observed values
#' @param pred A vector of predicted values. Same length as `obs`
#' @param type Character for the desired distribution
#' @param N Numeric, the sample size corresponding to `obs` for multinomial or Dirichlet multinomial distributions.
#' @param theta Numeric, the linear (`type = "dirmult1"`) or saturating (`type = "dirmult2"`) Dirichlet-multinomial parameter, respectively. See Thorson et al. (2017)
#' @param stdev Numeric or vectorized for `obs`, the likelihood standard deviation for lognormal or logit-normal distributions.
#' @return Numeric representing the log-likelihood.
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
#' like_comp(obs, pred, N = 10, type = "dirmult1", theta = 1)
#' like_comp(obs, pred, N = 10, type = "dirmult1", theta = 20)
#' @export
like_comp <- function(obs, pred, type = c("multinomial", "dirmult1", "dirmult2", "lognormal", "logitnormal"),
                      N = sum(obs), theta, stdev) {

  stopifnot(length(obs) == length(pred))
  type <- match.arg(type)

  if (all(is.na(obs)) || !sum(obs)) {

    v <- 0

  } else if (type == "multinomial") {

    pobs <- obs/sum(obs)
    v <- dmultinom(x = N * pobs, prob = pred, log = TRUE)

  } else if (type == "dirmult1") {

    alpha <- theta * N * pred/sum(pred)
    v <- ddirmnom(obs, size = N, alpha = alpha, log = TRUE)

  } else if (type == "dirmult2") {

    alpha <- theta * pred/sum(pred)
    v <- ddirmnom(obs, size = N, alpha = alpha, log = TRUE)

  } else if (type == "lognormal") {

    pobs <- obs/sum(obs)
    ppred <- pred/sum(pred)

    v <- dnorm(log(ppred/pobs), 0, stdev, log = TRUE) %>% sum()

  } else if (type == "logitnormal") {

    lobs <- qlogis(obs/sum(obs))
    lpred <- qlogis(pred/sum(pred))

    v <- dnorm(lobs, lpred, stdev, log = TRUE) %>% sum()
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

#' Likelihood for CKMR
#'
#' Returns the log-likelihood for a set of pairwise comparisons. For a parent-offspring pair, a comparison
#' is defined by the capture year of parent, capture age of parent, and birth year of offspring.
#' For a half-sibling pair, a comparison is defined by the cohort year of each sibling.
#' Binomial and Poisson distributions are supported (Conn et al. 2020).
#'
#' @param n The number of pairwise comparisons
#' @param m The number of kinship matches
#' @param p The probability of kinship match
#' @return Numeric representing the log-likelihood.
#' @seealso [calc_POP()] [calc_HSP()]
#' @references
#' Conn, P.B. et al. 2020. Robustness of close-kin mark–recapture estimators to dispersal
#' limitation and spatially varying sampling probabilities. Ecol. Evol. 10: 5558-5569. \doi{10.1002/ece3.6296}
#' @export
like_CKMR <- function(n, m, p, type = c("binomial", "poisson")) {
  type <- match(type)
  if (type == "binomial") {
    v <- dbinom(m, n, p, log = TRUE)
  } else if (type == "poisson") {
    v <- dpois(m, n * p, log = TRUE)
  }
  return(v)
}
