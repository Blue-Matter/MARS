
like_comp <- function(obs, pred, type = c("multinominal"), ...) {
  nbin <- length(obs)
  stopifnot(length(obs) == length(pred))
  dots <- list(...)

  if (type == "multinomial") {
    dmultinom(x = obs, size = dots$N, prob = pred, log = TRUE)
  }

}
