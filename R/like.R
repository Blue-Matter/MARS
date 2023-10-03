
#' Likelihood for composition vectors
#'
#' Returns the log-likelihood for composition data, e.g., length, age, or stock composition,
#' with various statistical distributions supported.
#'
#' @param obs A vector of observed values
#' @param pred A vector of predicted values. Same length as `obs`
#' @param type Character for the desired distribution
#' @param N Numeric, the sample size corresponding to `obs` for multinomial or Dirichlet multinomial distributions.
#' @param p Numeric, the linear (`type = "dirmult1"`) or saturating (`type = "dirmult2"`) Dirichlet-multinomial parameter, respectively. See Thorson et al. (2017)
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
#' like_comp(obs, pred, N = 10, type = "dirmult1", p = 1)
#' like_comp(obs, pred, N = 10, type = "dirmult1", p = 20)
#' @export
like_comp <- function(obs, pred,
                      type = c("multinomial", "dirmult1", "dirmult2", "lognormal", "logitnormal"),
                      N = sum(obs), p, stdev) {

  #if (is.matrix(obs)) obs <- colSums(obs)
  #if (is.matrix(pred)) pred <- colSums(pred)

  stopifnot(length(obs) == length(pred))
  type <- match.arg(type)
  dots <- list(...)

  if (all(is.na(obs)) || !sum(obs)) {

    v <- 0

  } else if (type == "multinomial") {

    pobs <- obs/sum(obs)
    v <- dmultinom(x = N * pobs, prob = pred, log = TRUE)

  } else if (type == "dirmult1") {

    alpha <- dots$p * N * pred/sum(pred)
    v <- ddirmnom(obs, size = N, alpha = alpha, log = TRUE)

  } else if (type == "dirmult2") {

    alpha <- dots$p * pred/sum(pred)
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
