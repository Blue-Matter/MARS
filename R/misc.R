
#' Quadratic penalty function
#'
#' Taped penalty function if `x < eps`
#'
#' @param x Numeric, the parameter
#' @param eps Numeric, the threshold below which a penalty will be applied
#'
#' @return
#' The penalty value is
#'
#' \deqn{
#' \textrm{penalty} =
#' \begin{cases}
#' 0.1 (x - \varepsilon)^2 & x \le \varepsilon\\
#' 0 & x > \varepsilon
#' \end{cases}
#' }
#'
#' @export
posfun <- function(x, eps) CondExpGe(x, eps, 0, 0.01 * (x - eps) * (x - eps))

#' Softmax function
#'
#' Takes a vector of real numbers and returns the corresponding vector of probabilities
#'
#' @param eta Vector
#' @param log Logical, whether to return the value of the logarithm
#'
#' @return A vector equal to length of `eta`: \eqn{\exp(\eta)/\sum\exp(\eta)}
#' @details Uses `MARS:::logspace.add` for numerical stability
#' @export
softmax <- function(eta, log = FALSE) {
  den <- Reduce(logspace.add, eta)
  v <- eta - den

  if (log) {
    v
  } else {
    exp(v)
  }
}

logspace.add <- function(lx, ly) pmax(lx, ly) + log1p(exp(-abs(lx - ly)))

#' Calculate covariance matrix
#'
#' Uses Cholesky factorization to generate a covariance matrix (or any symmetric positive definite matrix).
#'
#' @param sigma Numeric vector of marginal standard deviations (all greater than zeros)
#' @param lower_diag Numeric vector to populate the lower triangle of the correlation matrix. All real numbers.
#' Length `sum(1:(length(sigma) - 1))`
#' @examples
#' set.seed(23)
#' n <- 5
#' sigma <- runif(n, 0, 2)
#' lower_diag <- runif(sum(1:(n-1)), -10, 10)
#' Sigma <- conv_Sigma(sigma, lower_diag)
#' Sigma/t(Sigma) # Is symmetric matrix? All ones
#' cov2cor(Sigma)
#' @export
conv_Sigma <- function(sigma, lower_diag) {
  n <- length(sigma)
  stopifnot(length(lower_diag) == sum(1:(n-1)))

  # Parameterizes correlation matrix of X in terms of Cholesky factors
  # https://github.com/kaskr/RTMB/blob/master/tmb_examples/sdv_multi.R
  L <- diag(n)
  L[lower.tri(L)] <- lower_diag
  row_norms <- apply(L, 1, function(row) sqrt(sum(row * row)))
  L <- L / row_norms
  R <- L %*% t(L)  # Correlation matrix of X (guaranteed positive definite)

  V <- diag(sigma)
  Sigma <- V %*% R %*% V
  return(Sigma)
}


#' @importFrom stats nlminb rnorm
optimize_RTMB_model <- function(obj, use_hessian = FALSE, restart = 0, do_sd = TRUE,
                                control = list(), lower = -Inf, upper = Inf) {
  old_warn <- options()$warn
  options(warn = 2)
  on.exit(options(warn = old_warn))

  restart <- as.integer(restart)

  if (is.null(obj$env$random) && use_hessian) h <- obj$he else h <- NULL

  opt <- tryCatch(
    nlminb(obj$par, obj$fn, obj$gr, h, control = control, lower = lower, upper = upper),
    error = function(e) as.character(e)
  )

  if (do_sd) {
    SD <- get_sdreport(obj)

    if (!SD$pdHess && restart > 0) {
      if (!is.character(opt)) obj$par <- opt$par * exp(rnorm(length(opt$par), 0, 1e-3))
      Recall(obj, use_hessian, restart - 1, do_sd, control, lower, upper)
    } else {
      return(list(opt = opt, SD = SD))
    }
  } else {
    return(list(opt = opt, SD = NULL))
  }
}

check_det <- function(h, abs_val = 0.1, is_null = TRUE) {
  if (is.null(h)) return(is_null)
  det_h <- det(h) %>% abs()
  !is.na(det_h) && det_h < abs_val
}

#' @importFrom stats optimHess cov2cor
get_sdreport <- function(obj, getReportCovariance = FALSE, bias.correct = FALSE) {
  old_warn <- options()$warn
  options(warn = -1)
  on.exit(options(warn = old_warn))

  par.fixed <- obj$env$last.par.best

  if (is.null(obj$env$random)) {
    h <- obj$he(par.fixed)
    if (any(is.na(h)) || any(is.infinite(h)) || det(h) <= 0) {
      h <- NULL
    } else {
      res <- sdreport(obj, par.fixed = par.fixed, hessian.fixed = h,
                      bias.correct = bias.correct,
                      getReportCovariance = getReportCovariance)
      if (!res$pdHess) h <- NULL
    }
  } else {
    par.fixed <- par.fixed[-obj$env$random]
    #par.fixed <- par.fixed[obj$env$lfixed()]
    h <- NULL
  }

  if (check_det(h)) {  # If hessian doesn't exist or marginal positive-definite cases, with -0.1 < det(h) <= 0
    h <- optimHess(par.fixed, obj$fn, obj$gr)
    res <- sdreport(obj, par.fixed = par.fixed, hessian.fixed = h,
                    bias.correct = bias.correct, getReportCovariance = getReportCovariance)
  }

  if (check_det(h) && !res$pdHess && requireNamespace("numDeriv", quietly = TRUE)) {
    h <- numDeriv::jacobian(obj$gr, par.fixed)
    res <- sdreport(obj, par.fixed = par.fixed, hessian.fixed = h,
                    bias.correct = bias.correct, getReportCovariance = getReportCovariance)
  }

  if (all(is.na(res$cov.fixed)) && res$pdHess) {
    ch <- try(chol(h), silent = TRUE)
    if (!is.character(ch)) res$cov.fixed <- chol2inv(ch)
  }

  func <- attr(obj$env$data, "func")
  obj2 <- MakeADFun(func, obj$env$parameters, type = "ADFun", ADreport = TRUE, silent = obj$env$silent)
  gr <- obj2$gr(obj$env$last.par.best)

  res$env$corr.fixed <- cov2cor(res$cov.fixed) %>%
    round(3) %>%
    structure(dimnames = list(names(res$par.fixed), names(res$par.fixed)))

  return(res)
}

#' @importFrom TMB summary.sdreport
sdreport_int <- function(object, select = c("all", "fixed", "random", "report"), p.value = FALSE, ...) {
  if (is.character(object)) return(object)
  select <- match.arg(select, several.ok = TRUE)
  if ("all" %in% select) select <- c("fixed", "random", "report")
  if ("report" %in% select) {
    AD <- TMB::summary.sdreport(object, "report", p.value = p.value) %>% cbind("Gradient" = NA_real_)
  } else AD <- NULL

  if ("fixed" %in% select) {
    fix <- TMB::summary.sdreport(object, "fixed", p.value = p.value) %>% cbind("Gradient" = as.vector(object$gradient.fixed))
  } else fix <- NULL

  if (!is.null(object$par.random) && "random" %in% select) {
    random <- TMB::summary.sdreport(object, "random", p.value = p.value) %>% cbind("Gradient" = rep(NA_real_, length(object$par.random)))
  } else {
    random <- NULL
  }

  out <- rbind(AD, fix, random)
  out <- cbind(out, "CV" = ifelse(abs(out[, "Estimate"]) > 0, out[, "Std. Error"]/abs(out[, "Estimate"]), NA_real_))
  rownames(out) <- make_unique_names(rownames(out))
  return(out)
}


make_unique_names <- function(par_names) {
  unique_names <- unique(par_names)

  par_new <- lapply(unique_names, function(x) {
    ind <- par_names == x
    if (sum(ind) > 1) {
      paste0(x, "_", 1:sum(ind))
    } else x
  })
  do.call(c, par_new)
}



message <- function(...) {
  if (requireNamespace("usethis", quietly = TRUE)) {
    dots <- list(...)
    do.call(c, dots) %>% paste0(collapse = "") %>% usethis::ui_done()
  } else {
    base::message(...)
  }
}

warning <- function(...) {
  if (requireNamespace("usethis", quietly = TRUE)) {
    dots <- list(...)
    do.call(c, dots) %>% paste0(collapse = "") %>% usethis::ui_warn()
  } else {
    base::warning(...)
  }
}


stop <- function(..., call. = TRUE, domain = NULL) {
  if (requireNamespace("usethis", quietly = TRUE)) {
    dots <- list(...)
    do.call(c, dots) %>% paste0(collapse = "") %>% usethis::ui_stop()
  } else {
    base::stop(..., call. = call., domain = domain)
  }
}

