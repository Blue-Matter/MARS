
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

logspace.add <- function(lx, ly) CondExpGt(lx, ly, lx, ly) + log1p(exp(-abs(lx - ly)))

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

conv_steepness <- function(x, SRR = c("BH", "Ricker")) {
  SRR <- match.arg(SRR)
  switch(
    SRR,
    "BH" = 0.8 * plogis(x),
    "Ricker" = exp(x)
  ) + 0.2
}

conv_mat <- function(x, na) {
  a50 <- na * plogis(x[1])
  a95 <- a50 + exp(x[2])

  a <- seq(1, na)
  m <- 1/(1 + exp(-log(19) * (a - a50)/(a95 - a50)))
  return(m)
}


#' Optimize RTMB model
#'
#' A convenient function that fits a RTMB model and calculates standard errors.
#'
#' @param obj The list returned by [RTMB::MakeADFun()]
#' @param hessian Logical, whether to pass the Hessian function `obj$he` to [stats::nlminb()]. Only used if
#' there are no random effects in the model.
#' @param restart Integer, the maximum number of additional attempts to fit the model. See details.
#' @param do_sd Logical, whether to calculate standard errors through [get_sdreport()]
#' @param control List of options passed to [stats::nlminb()]
#' @param lower Lower bounds of parameters passed to [stats::nlminb()]
#' @param upper Upper bounds of parameters passed to [stats::nlminb()]
#' @param silent Logical, whether to report progress to console
#' @return A named list: "opt" is the output of [stats::nlminb()] and "SD" is the output of [get_sdreport()]
#' @details
#' Argument `restart` allows for recursive model fitting to obtain convergence, through the following procedure:
#' 1. Optimize model with [stats::nlminb()].
#' 2. Determine convergence, defined by [TMB::sdreport()] by whether the Cholesky decomposition of the covariance matrix is possible.
#' 3. If convergence is not achieved, jitter parameter estimates with multiplicative factor `rlnorm(mean = 0, sd = 1e-3)` and return to step 1.
#' @importFrom stats nlminb rnorm
#' @seealso [get_sdreport()]
#' @keywords internal
#' @export
optimize_RTMB <- function(obj, hessian = FALSE, restart = 0, do_sd = TRUE,
                          control = list(iter.max = 2e+05, eval.max = 4e+05),
                          lower = -Inf, upper = Inf, silent = FALSE) {
  old_warn <- options()$warn
  options(warn = -1)
  on.exit(options(warn = old_warn))

  restart <- as.integer(restart)

  if (is.null(obj$env$random) && hessian) h <- obj$he else h <- NULL

  if (!silent) message_info("Fitting model..")
  opt <- tryCatch(
    nlminb(obj$par, obj$fn, obj$gr, h, control = control, lower = lower, upper = upper),
    error = function(e) as.character(e)
  )

  if (do_sd) {
    SD <- get_sdreport(obj, silent = silent)

    if (!SD$pdHess && restart > 0) {
      if (!is.character(opt)) obj$par <- opt$par * exp(rnorm(length(opt$par), 0, 1e-3))
      Recall(obj, hessian, restart - 1, do_sd, control, lower, upper, silent)
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

#' Calculate standard errors
#'
#' A wrapper function to return standard errors. Various numerical techniques are employed to obtain
#' a positive-definite covariance matrix.
#' @inheritParams optimize_RTMB
#' @param getReportCovariance Logical, passed to [TMB::sdreport()]
#' @param silent Logical, whether to report progress to console. See details.
#' @param ... Other arguments to [TMB::sdreport()] besides `par.fixed, hessian.fixed, getReportCovariance`
#' @details
#' In marginal cases where the determinant of the Hessian matrix is less than 0.1, several steps are utilized to
#' obtain a positive-definite covariance matrix, in the following order:
#' 1. Invert hessian returned by `h <- obj@he(obj$env.last.par.best)` (skipped in models with random effects)
#' 2. Invert hessian returned by `h <- stats::optimHess(obj$env.last.par.best, obj$fn, obj$gr)`
#' 3. Invert hessian returned by `h <- numDeriv::jacobian(obj$gr, obj$env.last.par.best)`
#' 4. Calculate covariance matrix from `chol2inv(chol(h))`
#' @return
#' Object returned by [TMB::sdreport()]. A correlation matrix is generated and stored in: `get_sdreport(obj)$env$corr.fixed`
#' @importFrom stats optimHess
#' @export
get_sdreport <- function(obj, getReportCovariance = FALSE, silent = FALSE, ...) {
  old_warn <- options()$warn
  options(warn = -1)
  on.exit(options(warn = old_warn))

  old_comparison <- TapeConfig()["comparison"]
  on.exit(TapeConfig(comparison = old_comparison), add = TRUE)
  TapeConfig(comparison = "tape")

  par.fixed <- obj$env$last.par.best

  if (is.null(obj$env$random)) {
    h <- obj$he(par.fixed)
    if (any(is.na(h)) || any(is.infinite(h)) || check_det(h)) {
      h <- NULL
    } else {
      if (!silent) message_info("Calculating standard errors with hessian from obj$he()..")
      res <- sdreport(obj, par.fixed = par.fixed, hessian.fixed = h,
                      getReportCovariance = getReportCovariance, ...)
    }
  } else {
    par.fixed <- par.fixed[-obj$env$random]
    #par.fixed <- par.fixed[obj$env$lfixed()]
    h <- NULL
  }

  if (is.null(h) || check_det(h)) {  # If hessian doesn't exist or marginal positive-definite cases, with -0.1 < det(h) <= 0
    if (!silent) message_info("Calculating standard errors with hessian from stats::optimHess()..")
    h <- optimHess(par.fixed, obj$fn, obj$gr)
    res <- sdreport(obj, par.fixed = par.fixed, hessian.fixed = h,
                    getReportCovariance = getReportCovariance, ...)
  }

  if (check_det(h) && !res$pdHess && requireNamespace("numDeriv", quietly = TRUE)) {
    if (!silent) message_info("Calculating standard errors with hessian from numDeriv::jacobian()..")
    h <- numDeriv::jacobian(obj$gr, par.fixed)
    res <- sdreport(obj, par.fixed = par.fixed, hessian.fixed = h,
                    getReportCovariance = getReportCovariance, ...)
  }

  if (all(is.na(res$cov.fixed)) && res$pdHess) {
    if (!silent) message_info("Calculating standard errors from chol2inv(chol(h))..")
    ch <- try(chol(h), silent = TRUE) # Not needed, this is the test for convergence in sdreport
    if (!is.character(ch)) res$cov.fixed <- chol2inv(ch)
  }

  fixed.names <- make_unique_names(res, select = "fixed")

  res$env$corr.fixed <- cov2cor(res$cov.fixed) %>% round(3) %>%
    structure(dimnames = list(fixed.names, fixed.names))

  res$env$hessian <- round(h, 3) %>%
    structure(dimnames = list(fixed.names, fixed.names))

  if (!res$pdHess) {
    if (!silent) {
      message_oops("Check convergence. Covariance matrix is not positive-definite.")
      message_info("Maximum gradient is ", round(max(abs(res$gradient.fixed)), 5))
      if (exists("h", inherits = FALSE) && !is.null(h)) {
        message_info("Determinant of Hessian is ", round(det(h), 5))
      }
    }
  }

  return(res)
}

#' @importFrom TMB summary.sdreport
sdreport_int <- function(object, select = c("all", "fixed", "random", "report"), p.value = FALSE, ...) {
  if (is.character(object)) return(object)
  select <- match.arg(select, several.ok = TRUE)
  if ("all" %in% select) select <- c("fixed", "random", "report")
  if ("report" %in% select) {
    AD <- TMB::summary.sdreport(object, "report", p.value = p.value) %>% cbind("Gradient" = NA_real_)
    ADnames <- make_unique_names(object, select = "report")
  } else AD <- ADnames <- NULL

  if ("fixed" %in% select) {
    fix <- TMB::summary.sdreport(object, "fixed", p.value = p.value) %>% cbind("Gradient" = as.vector(object$gradient.fixed))
    fixnames <- make_unique_names(object, select = "fixed")
  } else fix <- fixnames <- NULL

  if (!is.null(object$par.random) && "random" %in% select) {
    random <- TMB::summary.sdreport(object, "random", p.value = p.value) %>% cbind("Gradient" = rep(NA_real_, length(object$par.random)))
    randomnames <- make_unique_names(object, select = "random")
  } else {
    random <- randomnames <- NULL
  }

  out <- rbind(AD, fix, random)
  out <- cbind(out, "CV" = ifelse(abs(out[, "Estimate"]) > 0, out[, "Std. Error"]/abs(out[, "Estimate"]), NA_real_))
  rownames(out) <- c(ADnames, fixnames, randomnames)
  return(out)
}


#' Retrieve data object used to fit model
#'
#' A convenient function to retrieve the data object used to fit the model. The object is embedded in an environment
#' within the RTMB object.
#'
#' @param MARSassess [MARSassess-class] object returned by `fit_MARS()`
#' @return [MARSdata-class] object
#' @export
get_MARSdata <- function(MARSassess) {
  func <- attr(MARSassess@obj$env$data, "func")
  MARSdata <- get("MARSdata", envir = environment(func), inherits = FALSE)
  return(MARSdata)
}

make_unique_names <- function(x, select = c("fixed", "random", "report")) {
  select <- match.arg(select)

  if (select == "fixed") {
    par_names <- unique(names(x$par.fixed))
    par_list <- as.list(x, what = "Estimate", report = FALSE)
  } else if (select == "report") {
    par_names <- unique(names(x$value))
    par_list <- as.list(x, what = "Estimate", report = TRUE)
  } else {
    par_names <- unique(names(x$random))
    par_list <- as.list(x, what = "Estimate", report = FALSE)
  }


  par_dims <- lapply(par_names, function(y) {
    dim_y <- dim(par_list[[y]])
    if (is.null(dim_y)) dim_y <- length(par_list[[y]])
    ind_y <- lapply(dim_y, function(i) seq(1, i))

    est_grid <- do.call(expand.grid, ind_y)

    if (select != "report") {
      map_y <- attr(x$env$parameters[[y]], "map")

      if (!is.null(map_y)) {
        est_y <- map_y >= 0 & !duplicated(map_y)
        est_grid <- est_grid[est_y, ]
      }
    }

    dim_char <- sapply(1:nrow(est_grid), function(i) {
      paste0(
        "[",
        paste0(est_grid[i, ], collapse = ", "),
        "]"
      )
    })
    paste0(y, dim_char)
  })

  do.call(c, par_dims)
}


make_yearseason <- function(year, nm = 4) {
  if (nm <= 1) return(year)
  year_long <- lapply(year, function(y) y + (1:nm - 1)/nm)
  do.call(c, year_long)
}

# x must be a three dimensional array due to rbind()
collapse_yearseason <- function(x, MARGIN = c(1, 2)) {

  if (length(dim(x)) > 2) {
    ny <- dim(x)[MARGIN[1]]
    nm <- dim(x)[MARGIN[2]]

    xout <- lapply(1:ny, function(y) {
      out <- array(NA_real_, c(nm, dim(x)[-MARGIN]))
      comma <- rep(", ", length(dim(x)) - 1) %>% paste0(collapse = "")
      eval(parse(text = paste0("out[] <- x[y", comma, "]")))
      return(out)
    })

    do.call(rbind, xout)
  } else {
    as.numeric(t(x))
  }
}

message <- function(...) {
  if (requireNamespace("usethis", quietly = TRUE)) {
    dots <- list(...)
    do.call(c, dots) %>% paste0(collapse = "") %>% usethis::ui_done()
  } else {
    base::message(...)
  }
}


message_info <- function(...) {
  if (requireNamespace("usethis", quietly = TRUE)) {
    dots <- list(...)
    do.call(c, dots) %>% paste0(collapse = "") %>% usethis::ui_info()
  } else {
    base::message(...)
  }
}

message_oops <- function(...) {
  if (requireNamespace("usethis", quietly = TRUE)) {
    dots <- list(...)
    do.call(c, dots) %>% paste0(collapse = "") %>% usethis::ui_oops()
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

