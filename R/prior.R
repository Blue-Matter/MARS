

prior_sel <- function() {

}


prior_R0x <- function() {

}

alphaconv <- function(m, sd) m * (((m * (1 - m))/(sd^2)) - 1)
betaconv <- function(m, sd) (1 - m) * (((m * (1 - m))/(sd^2)) - 1)
sdconv <- function(m, sd) (log(1 + ((sd^2)/(m^2))))^0.5

#' @name prior
#'
#' @title Priors for MARS model
#' @description
#' Priors in MARS are set by providing character strings which are then parsed into an expression and evaluated in the model environment
#' (see example). This provides flexibility to set a prior for any desired model parameter or variable. See list of parameters in the
#' documentation for `[check_parameters()]` for options (note that priors for `log_rdev_ys` and `log_initrdev_as` are not needed as they're hard-coded into the model).
#' Several functions below generate the character string for the prior for important dynamics parameters, such as natural mortality and steepness.
#'
#' @param MARSdata Data object. Class \linkS4class{MARSdata}
#' @return Character.
#' @examples
#' # Add M and steepness prior to model
#'
#' dat <- new("MARSdata")
#' dat@Dmodel@ns <- 1
#' dat@Dstock@SRR_s <- "BH"
#'
#' pr_M <- prior_M(dat, s = 1, log(0.25), 0.3)
#' pr_h <- prior_h(dat, s = 1, 0.8, 0.15)
#' dat@Dmodel@prior <- c(pr_M, pr_h)
NULL

#' @name prior
#' @details - `prior_h` returns the log prior for stock-recruit steepness. Beta distribution for the Beverton-Holt function and normal distribution for Ricker function.
#' @param s Integer for stock
#' @param m Mean in un-transformed space
#' @param stdev Standard deviation in un-transformed space
#' @export
prior_h <- function(MARSdata, s = 1, m, stdev) {
  SRR <- match.arg(MARSdata@Dstock@SRR_s[s], choices = c("BH", "Ricker"))

  if (SRR == "BH") {

    m_beta <- 1.25 * m - 0.25
    stdev_beta <- 1.25 * stdev

    a <- alphaconv(m_beta, stdev_beta) %>% round(3)
    b <- betaconv(1.25 * m - 0.25, 1.25 * stdev) %>% round(3)

    if (a <= 0) stop("Beta distribution alpha parameter < 0. Try reducing the prior SD.", call. = FALSE)
    if (b <= 0) stop("Beta distribution beta parameter < 0. Try reducing the prior SD.", call. = FALSE)

    txt <- paste0("dbeta(plogis(p$t_h_s[", s, "]), ", a, ", ", b, ", log = TRUE) + log(plogis(p$t_h_s[", s, "]) - plogis(p$t_h_s[", s, "])^2)")

  } else {
    txt <- paste0("dnorm(exp(p$t_h_s[", s, "]) + 0.2, ", round(m, 3), ", ", round(stdev, 3), ", log = TRUE) + p$t_h_s[", s, "]")
  }
  return(txt)
}

#' @name prior
#' @details - `prior_M` returns the log prior for natural mortality. Lognormal distribution.
#' @aliases prior_M
#' @param meanlog Mean of the lognormal distribution on the log scale
#' @param sdlog Standard of the lognormal distribution on the log scale
#' @export
prior_M <- function(MARSdata, s = 1, meanlog, sdlog) {
  txt <- paste0("dnorm(p$log_M_s[", s, "], ", round(meanlog, 3), ", ", round(sdlog, 3), ", log = TRUE)")
  return(txt)
}

#' @name prior
#' @details - `prior_q` returns the log prior for index catchability. Lognormal distribution.
#' @aliases prior_q
#' @param i Integer for the corresponding index
#' @export
prior_q <- function(MARSdata, i = 1, meanlog, sdlog) {
  txt <- paste0("dnorm(log(q_i[", i, "]), ", round(meanlog, 3), ", ", round(sdlog, 3), ", log = TRUE) - log(q_i[", i, "])")
  return(txt)
}
