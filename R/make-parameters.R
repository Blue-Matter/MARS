
#' Make list of parameters for RTMB
#'
#' @description
#' Sets up the list of parameters, map of parameters (see `map` argument in [TMB::MakeADFun()]), and identifies some random effects parameters
#' based on the input data and some user choices on model configuration.
#'
#' These functions provide a template for the parameter and map setup that can be adjusted for alternative configurations. [check_parameters()]
#' checks whether custom made parameter lists are of the correct dimension.
#' @seealso [MARSdata-class]
#'
#' @section Parameters:
#'
#' Generally parameter names will have up to three components, separated by underscores.
#' For example, `log_M_s` represents the natural logarithm of natural mortality, and is a vector by stock `s`.
#'
#' The first component describes the transformation from the estimated parameter space to the normal parameter space,
#' frequently, `log` or `logit`. Prefix `t` indicates some other custom transformation that is described below.
#'
#' Second is the parameter name, e.g., `M` for natural mortality, `rdev` for recruitment deviates, etc.
#'
#' Third is the dimension of the parameter variable and the indexing for the vectors, matrices, and arrays, e.g., `y` for year, `s` for stock.
#' See [MARSdata-class]. Here, an additional index `p` represents some other number of parameters that is described below.
#'
#' \describe{
#' \item{`t_R0_s`}{Vector by `s`. Unfished recruitment, i.e., intersection of unfished replacement line and average stock recruit function,
#' is represented as: `R0_s <- exp(t_R0_s) * MARSdata@Dmodel@scale_s`. By default, `t_R0_s = 3`}
#' \item{`t_h_s`}{Vector by `s`. Steepness of the stock-recruit function. Logit space for Beverton-Holt and log space for Ricker functions.
#' Default steepness value of 0.8}
#' \item{`mat_ps`}{Matrix `[2, s]`. Maturity parameters (can be estimated or specified in data object). Logistic functional form. The
#' parameter in the first row is the age of 50 percent maturity in logit space: `a50_s <- plogis(mat_ps[1, ] * na)`.
#' In the second row is the age of 95 percent maturity as a logarithmic offset: `a95_s <- exp(a50_s + mat_ps[2, ])`.
#' Default `a50_s <- 0.5 * na` and `a95_s <- a50_s + 1`}
#' \item{`log_M_s`}{Vector by `s`. Natural logarithm of natural mortality (can be estimated or specified in data object).
#' Default parameter value for all stocks: `M <- -log(0.05)/MARSdata@Dmodel@na`}
#' \item{`log_rdev_ys`}{Matrix `[y, s]`. Log recruitment deviations. By default, all start values are at zero.}
#' \item{`log_sdr_s`}{Vector by `s`. log-Standard deviation of the log recruitment deviations. Default SD = 0.4}
#' \item{`log_q_fs`}{Matrix `[f, s]`. The natural logarithm of `q_fs`, the relative fishing efficiency of `f` for stock `s`.
#' Equal values imply equal catchability of all stocks. See equations in [calc_F()]. Default sets all values to zero.}
#' \item{`sel_pf`}{Matrix `[3, f]`. Fishery selectivity parameters in logit or log space. See equations [conv_selpar()], where `sel_pf` is the `x` matrix.}
#' \item{`sel_pi`}{Matrix `[3, i]`. Index selectivity parameters in logit or log space. See equations [conv_selpar()], where `sel_pi` is the `x` matrix.}
#' \item{`mov_x_marrs`}{Array `[m, a, r, r, s]`. Base movement matrix. Set to -1000 to effectively exclude movements from region pairs.
#' See equations in [conv_mov()]}
#' \item{`mov_g_ymars`}{Array `[y, m, a, r, s]`. Attractivity term in gravity model for movement. If `x` and `v` are zero,
#' this matrix specifies the distribution of total stock abundance into the various regions. See equations in [conv_mov()]}
#' \item{`mov_v_ymas`}{Array `[y, m, a, s]`. Viscosity term in gravity model for movement. See equations in [conv_mov()]}
#' \item{`log_sdg_rs`}{Array `[r, s]`. Marginal log standard deviation in the stock distribution (`mov_g_ymars`) among regions for stock `s`.
#' Only used when `est_mov = "dist_random"`. Default SD of 0.1.}
#' \item{`t_corg_ps`}{Array `[sum(1:(nr - 1)), s]`. Lower triangle of the correlation matrix for `mov_g_ymars`, to be obtained with the
#' Cholesky factorization. Only used when `est_mov = dist_random`. Default values of zero.}
#' \item{`log_initF_mfr`}{Array `[m, f, r]`. Initial F corresponding to the equilibrium catch.}
#' \item{`log_initrdev_as`}{Array `[a, s]`. Recruitment deviations for the initial abundance-at-age vector.}
#' }
#'
#' @param MARSdata S4 data object
#' @param start An optional list of parameters. Named list of parameters with the associated dimensions and transformations below.
#' Overrides default values created by [make_parameters()].
#' @param silent Logical, whether to report messages to the console
#' @param ... Various arguments for [make_map()] (could be important!)
#' @return
#' [make_parameters()] returns a list of parameters (`"p"`) concatenated with the output of [make_map()].
#' @export
make_parameters <- function(MARSdata, start = list(), silent = FALSE, ...) {

  getAllS4(MARSdata@Dmodel)

  p <- start

  # Stock parameters ----
  if (is.null(p$t_R0_s)) p$t_R0_s <- rep(3, ns)
  if (is.null(p$t_h_s)) {
    p$t_h_s <- local({
      h <- 0.8
      ifelse(MARSdata@Dstock@SRR_s == "BH", qlogis((h - 0.2)/0.8), log(h - 0.2))
    })
  }
  if (is.null(p$mat_ps)) {
    p$mat_ps <- sapply(1:ns, function(s) {
      a50 <- 0.5 * na
      a95 <- a50 + 1
      logit_a50 <- qlogis(a50/na)
      log_diff <- log(a95 - a50)
      c(logit_a50, log_diff)
    })
  }

  if (is.null(p$log_M_s)) p$log_M_s <- rep(-log(0.05)/na, ns)
  if (is.null(p$log_rdev_ys)) p$log_rdev_ys <- matrix(0, ny, ns)
  if (is.null(p$log_sdr_s)) p$log_sdr_s <- rep(0.4, ns)

  if (is.null(p$mov_x_marrs)) p$mov_x_marrs <- array(0, c(nm, na, nr, nr, ns))
  if (is.null(p$mov_g_ymars)) p$mov_g_ymars <- array(0, c(ny, nm, na, nr, ns))
  if (is.null(p$mov_v_ymas)) p$mov_v_ymas <- array(0, c(ny, nm, na, ns))
  if (is.null(p$log_sdg_rs)) p$log_sdg_rs <- array(log(0.1), c(nr, ns))
  if (is.null(p$t_corg_ps)) p$t_corg_ps <- array(0, c(sum(1:(nr - 1)), ns))

  # Fleet parameters ----
  if (is.null(p$log_q_fs)) p$log_q_fs <- matrix(0, nf, ns)
  if (is.null(p$sel_pf)) {
    p$sel_pf <- local({
      sel <- matrix(NA, 3, max(MARSdata@Dfishery@sel_block_yf))
      sel[1, ] <- 0 # Apical sel between 0 and max age/length
      sel[2, ] <- log(0.01) # Knife edge ascending selectivity
      sel[3, ] <- log(0.1) # More sloping descending limb of selectivity
      sel
    })
  }

  # Index parameters ----
  if (is.null(p$sel_pi)) {
    p$sel_pi <- local({
      sel <- matrix(NA, 3, MARSdata@Dsurvey@ni)
      sel[1, ] <- 0 # Apical sel between 0 and max age/length
      sel[2, ] <- log(0.01) # Knife edge ascending selectivity
      sel[3, ] <- log(0.1) # More sloping descending limb of selectivity
      sel
    })
  }

  # Initial conditions ----
  if (is.null(p$log_initF_mfr)) {
    p$log_initF_mfr <- ifelse(MARSdata@Dfishery@Cinit_mfr < 1e-8, -1000, log(0.1))
  }
  if (is.null(p$log_initrdev_as)) p$log_initrdev_as <- matrix(0, na, ns)

  do_map <- make_map(p, MARSdata, ...)
  out <- c(list(p = p), do_map)
  out$p <- check_parameters(out$p, out$map, MARSdata, silent)

  return(out)
}

#' @rdname make_parameters
#' @aliases make_map
#' @param p List of parameters, e.g., returned by [make_parameters()]
#' @param est_M Logical, estimate natural mortality?
#' @param est_h Logical, estimate steepness?
#' @param est_mat Logical, estimate maturity?
#' @param est_sdr Logical, estimate standard deviation of recruitment deviates?
#' @param est_mov Character describing structure of stock movement parameters. See details below.
#' @param est_qfs Logical, estimate relative catchability of stocks by each fleet? Fix `log_q_fs` for the first stock if `TRUE`
#' @section Movement setup for `make_map()`:
#' If a single region model or `est_mov = "none"`: no movement parameters are estimated.
#'
#' If `est_mov = "dist_random"`: fix all values for `mov_x_marrs` and `mov_v_ymas`. Fix `mov_g_ymars` for the first region for each year,
#' season, age, and stock. `mov_g_ymars` are random effects.
#'
#' If `est_mov = "gravity_fixed"`: fix all values for `mov_x_marrs`. Fix `mov_g_ymars` for the first region for each year,
#' season, age, and stock. Estimate all `mov_v_ymas`. Both `mov_g_ymars` and `mov_v_ymas` are fixed effects.
#'
#' @importFrom dplyr filter
#' @importFrom rlang .data .env
#' @return
#' [make_map()] returns a named list containing parameter mappings (`"map"`) and a character vector of random effects (`"random"`).
#' @export
make_map <- function(p, MARSdata,
                     est_M = FALSE, est_h = FALSE, est_mat = FALSE, est_sdr = FALSE,
                     est_mov = c("none", "dist_random", "gravity_fixed"),
                     est_qfs = FALSE,
                     silent = FALSE) {

  est_mov <- match.arg(est_mov)
  getAllS4(MARSdata)
  getAllS4(MARSdata@Dmodel)

  random <- NULL
  map <- list()

  # Stock parameters ----
  if (est_M) {
    if (!silent) message_info("Estimating natural mortality")
  } else {
    map$log_M_s <- factor(rep(NA, ns))
    if (!length(Dstock@Md_yas)) stop("Natural mortality is not estimated. Need M values in the Dstock data object.")
  }
  if (est_h) {
    if (!silent) message_info("Estimating steepness")
  } else {
    map$t_h_s <- factor(rep(NA, ns))
  }
  if (est_mat) {
    if (!silent) message_info("Estimating maturity ogive")
  } else {
    map$mat_ps <- factor(array(NA, dim(p$mat_ps)))
    if (!length(Dstock@matd_yas)) stop("Maturity ogive is not estimated. Need maturity at age values in the Dstock data object.")
  }

  if (est_sdr) {
    if (!silent) message_info("Estimating sigma_R (SD of recruitment deviates)")
  } else {
    map$log_sdr_s <- factor(rep(NA, ns))
  }

  if (nr == 1 || est_mov == "none") {
    map$mov_x_marrs <- factor(array(NA, dim(p$mov_x_marrs)))
    map$mov_g_ymars <- factor(array(NA, dim(p$mov_g_ymars)))
    map$mov_v_ymas <- factor(array(NA, dim(p$mov_v_ymas)))

    map$log_sdg_rs <- factor(array(NA, dim(p$log_sdg_rs)))
    map$t_corg_ps <- factor(array(NA, dim(p$t_corg_ps)))

    #if (!silent) message_info("No movement parameters are estimated")

  } else if (est_mov == "dist_random") {

    map$mov_x_marrs <- factor(array(NA, dim(p$mov_x_marrs)))

    gval <- expand.grid(r = 2:nr, y = 1:ny, m = 1:nm, a = 1:na, s = 1:ns)
    gval$g <- seq(1, nrow(gval))
    map$mov_g_ymars <- sapply2(1:ns, function(s) {
      sapply2(1:na, function(a) {
        sapply2(1:nm, function(m) {
          sapply(1:ny, function(y) {
            gval_ymas <- filter(gval,
              .data$y == .env$y, .data$m == .env$m, .data$a == .env$a, .data$s == .env$s
            )
            g <- rep(NA, nr)
            g[gval_ymas$r] <- gval_ymas$g
          })
        })
      })
    }) %>%
      aperm(c(2, 3, 4, 1, 5)) %>%
      factor()

    map$mov_v_ymas <- factor(array(NA, dim(p$mov_v_ymas)))
    random <- c(random, "mov_g_ymars")

    if (!silent) message_info("Stock distribution is estimated as a random effect")

  } else {
    map$mov_x_marrs <- factor(array(NA, dim(p$mov_x_marrs)))

    gval <- expand.grid(r = 2:nr, y = 1:ny, m = 1:nm, a = 1:na, s = 1:ns)
    gval$g <- seq(1, nrow(gval))
    map$mov_g_ymars <- sapply2(1:ns, function(s) {
      sapply2(1:na, function(a) {
        sapply2(1:nm, function(m) {
          sapply(1:ny, function(y) {
            gval_ymas <- filter(
              gval, .data$y == .env$y, .data$m == .env$m, .data$a == .env$a, .data$s == .env$s
            )
            g <- rep(NA, nr)
            g[gval_ymas$r] <- gval_ymas$g
          })
        })
      })
    }) %>%
      aperm(c(2, 3, 4, 1, 5)) %>%
      factor()

    map$log_sdg_rs <- factor(array(NA, dim(p$log_sdg_rs)))
    map$t_corg_ps <- factor(array(NA, dim(p$t_corg_ps)))

    if (!silent) message_info("Stock movement is estimated as a fixed effect")
  }

  # Fleet parameters ----
  if (!est_qfs) {
    map$log_q_fs <- factor(matrix(NA, nf, ns))
  } else { # Fix q for first stock
    map$log_q_fs <- local({
      np <- ns * (nf - 1)
      m <- matrix(NA, nf, ns)
      m[, -1] <- 1:np
      factor(m)
    })
  }

  ## Fix dome parameter if selectivity is logistic or all parameters if mirrored to maturity
  if (any(!grepl("dome", Dfishery@sel_f))) {
    map$sel_pf <- local({
      n <- max(Dfishery@sel_block_yf)
      m <- sapply(1:n, function(f) {
        sel_f <- Dfishery@sel_f[f]
        vec <- rep(TRUE, 3)
        if (sel_f %in% c("logistic_age", "logistic_length")) vec[3] <- NA
        if (sel_f %in% c("B", "SB")) vec[] <- NA
        return(vec)
      })
      m[!is.na(m)] <- 1:sum(m, na.rm = TRUE)
      factor(m)
    })
  }

  # Survey parameters ----
  ## Fix dome parameter or all parameters or all parameters if mirrored to fleet or maturity
  if (any(!grepl("dome", Dsurvey@sel_i))) {
    old_warn <- options()$warn
    options(warn = -1)
    on.exit(options(warn = old_warn))

    map$sel_pi <- local({
      m <- sapply(1:Dsurvey@ni, function(i) {
        sel_i <- Dsurvey@sel_i[i]
        vec <- rep(TRUE, 3)
        if (sel_i %in% c("logistic_age", "logistic_length")) vec[3] <- NA
        if (sel_i %in% c("B", "SB") || !is.na(as.integer(sel_i))) vec[] <- NA
        return(vec)
      })
      m[!is.na(m)] <- 1:sum(m, na.rm = TRUE)
      factor(m)
    })
  }

  # Initial conditions ----
  if (any(Dfishery@Cinit_mfr < 1e-8)) {
    map$log_initF_mfr <- local({
      m <- ifelse(Dfishery@Cinit_mfr < 1e-8, NA, TRUE)
      m[!is.na(m)] <- 1:sum(m, na.rm = TRUE)
      factor(m)
    })
  }
  map$log_initrdev_as <- factor(matrix(NA, na, ns))

  return(list(map = map, random = random))
}

#' @rdname make_parameters
#' @param map List of mapped parameters. Used by [check_parameters()] only to count parameters.
#' @return
#' [check_parameters()] invisibly returns the parameter list if no problems are encountered.
#' @export
check_parameters <- function(p = list(), map, MARSdata, silent = FALSE) {

  if (!silent && !missing(map)) {
    vars <- names(p)
    npar <- lapply(vars, function(x) {
      if (is.null(map[[x]])) {
        nparam <- length(p[[x]])
      } else {
        m_x <- map[[x]]
        nparam <- length(unique(m_x[!is.na(m_x)]))
      }
      if (nparam) {
        data.frame(Variable = x, Parameters = nparam)
      } else {
        NULL
      }
    })
    message_info("Estimated parameters:")
    print(do.call(rbind, npar))
  }

  return(invisible(p))
}
