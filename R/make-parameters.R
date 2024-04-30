
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
#' In the second row is the age of 95 percent maturity as a logarithmic offset: `a95_s <- a50_s + exp(mat_ps[2, ])`.
#' Default `a50_s <- 0.5 * na` and `a95_s <- a50_s + 1`}
#' \item{`log_M_s`}{Vector by `s`. Natural logarithm of natural mortality (can be estimated or specified in data object).
#' Default parameter value for all stocks: `M <- -log(0.05)/MARSdata@Dmodel@na`}
#' \item{`log_rdev_ys`}{Matrix `[y, s]`. Log recruitment deviations. By default, all start values are at zero.}
#' \item{`log_sdr_s`}{Vector by `s`. log-Standard deviation of the log recruitment deviations. Default SD = 0.4}
#' \item{`log_q_fs`}{Matrix `[f, s]`. The natural logarithm of `q_fs`, the relative fishing efficiency of `f` for stock `s`.
#' Equal values imply equal catchability of all stocks. See equations in [calc_F()]. Default sets all values to zero.}
#' \item{`log_Fdev_ymfr`}{Array `[y, m, f, r]`. Fishing mortality parameters. For each fleet, the log of F is estimated directly for the
#' reference year, season, region. For other strata, F is an offset from this value:
#' \deqn{
#' F_{y,m,f,r} = \begin{cases}
#' \exp(x^{\textrm{Fmult}}_f) \quad & y = y_{\textrm{ref}}, m = m_{\textrm{ref}}, r = r_{\textrm{ref}}\\
#' \exp(x^{\textrm{Fmult}}_f + x^{\textrm{Fdev}}_{y,m,r}) \quad & \textrm{otherwise}
#' \end{cases}
#' }
#' }
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
#' @section Start list:
#' Users can provide `R0_s` and `h_s` in the start list. [make_parameters()] will make the appropriate transformation for the starting values
#' of `t_R0_s` and `t_h_s`, respectively.
#'
#' @param MARSdata S4 data object
#' @param start An optional list of parameters. Named list of parameters with the associated dimensions and transformations below.
#' Overrides default values created by [make_parameters()].
#' @param silent Logical, whether [make_map()] reports messages to the console
#' @param ... Various arguments for [make_map()] (could be important!)
#' @return
#' [make_parameters()] returns a list of parameters (`"p"`) concatenated with the output of [make_map()].
#' @importFrom stats approx
#' @export
make_parameters <- function(MARSdata, start = list(), map = list(), silent = FALSE, ...) {

  getAllS4(MARSdata@Dmodel)
  nf <- MARSdata@Dfishery@nf

  p <- start

  # Stock parameters ----
  if (!is.null(start$R0_s)) {
    p$t_R0_s <- log(start$R0_s/MARSdata@Dmodel@scale_s)
    p$R0_s <- NULL
  } else if (is.null(start$t_R0_s)) {
    p$t_R0_s <- rep(3, ns)
  }

  if (is.null(start$h_s)) {
    start$h_s <- rep(0.8, ns)
  } else {
    p$h_s <- NULL
  }
  p$t_h_s <- ifelse(MARSdata@Dstock@SRR_s == "BH", qlogis((start$h_s - 0.2)/0.8), log(start$h_s - 0.2))

  if (is.null(p$mat_ps)) {
    p$mat_ps <- sapply(1:ns, function(s) {
      a50 <- 0.5*na
      a95 <- a50 + 1
      logit_a50 <- qlogis(a50/na)
      log_diff <- log(a95 - a50)
      c(logit_a50, log_diff)
    })
  }

  if (is.null(p$log_M_s)) p$log_M_s <- rep(log(-log(0.05)/na), ns)
  if (is.null(p$log_rdev_ys)) p$log_rdev_ys <- matrix(0, ny, ns)
  if (is.null(p$log_sdr_s)) p$log_sdr_s <- rep(log(0.4), ns)

  if (is.null(p$log_recdist_rs)) p$log_recdist_rs <- matrix(0, nr, ns)

  if (is.null(p$mov_x_marrs)) {
    p$mov_x_marrs <- array(0, c(nm, na, nr, nr, ns))
    if (any(!MARSdata@Dstock@presence_rs)) {
      for(s in 1:ns) {
        presence_r <- MARSdata@Dstock@presence_rs[, s]
        if (any(!presence_r)) p$mov_x_marrs[, , presence_r, presence_r, s] <- -1000
      }
    }
  }
  if (is.null(p$mov_g_ymars)) p$mov_g_ymars <- array(0, c(ny, nm, na, nr, ns))
  if (is.null(p$mov_v_ymas)) p$mov_v_ymas <- array(0, c(ny, nm, na, ns))
  if (is.null(p$log_sdg_rs)) p$log_sdg_rs <- array(log(0.1), c(nr, ns))
  if (is.null(p$t_corg_ps)) p$t_corg_ps <- array(0, c(sum(1:(nr - 1)), ns))

  # Fleet parameters ----
  if (is.null(p$log_q_fs)) {
    p$log_q_fs <- matrix(0, nf, ns)
  }
  if (is.null(p$log_Fdev_ymfr)) {
    if (condition == "F") {
      p$log_Fdev_ymfr <- sapply2(1:nr, function(r) {
        sapply2(1:nf, function(f) {
          sapply(1:nm, function(m) {
            sapply(1:ny, function(y) {
              Fmult_y <- y == y_Fmult_f[f]
              Fmult_m <- m == m_Fmult_f[f]
              Fmult_r <- r == r_Fmult_f[f]
              if (Fmult_y && Fmult_m && Fmult_r) {
                log(-log(0.05)/na/nm)
              } else {
                0
              }
            })
          })
        })
      })
      p$log_Fdev_ymfr[MARSdata@Dfishery@Cobs_ymfr < 1e-8] <- -1000
    } else {
      p$log_Fdev_ymfr <- array(0, c(ny, nm, nf, nr))
    }
  }

  if (is.null(p$sel_pf)) {
    p$sel_pf <- sapply(unique(MARSdata@Dfishery@sel_block_yf), function(b) {
      sel_b <- MARSdata@Dfishery@sel_f[b]
      val <- numeric(3)
      if (grepl("length", sel_b) && length(MARSdata@Dfishery@CALobs_ymlfr)) {
        f_yb <- MARSdata@Dfishery@sel_block_yf == b

        CAL <- sapply2(1:nf, function(f) {
          sapply(1:MARSdata@Dmodel@ny, function(y) {
            if (f_yb[y, f]) {
              apply(MARSdata@Dfishery@CALobs_ymlfr[y, , , f, , drop = FALSE], 3, sum)
            } else {
              rep(0,  MARSdata@Dmodel@nl)
            }
          })
        }) %>% apply(1, sum)

        if (sum(CAL)) {
          LFS <- min(MARSdata@Dmodel@lmid[which.max(CAL)], 0.75 * max(MARSdata@Dmodel@lmid))
          L5 <- approx(cumsum(CAL)/sum(CAL), MARSdata@Dmodel@lmid, 0.05)$y

          if (L5 < LFS) {
            sigma_asc <- min((LFS - L5)/sqrt(-log(0.05, 2)), 0.25 * diff(range(MARSdata@Dmodel@lmid)))
            val[2:3] <- log(sigma_asc)
            val[1] <- qlogis(LFS/max(0.95 * MARSdata@Dmodel@lmid))
          }
        }
      }
      if (all(!val)) {
        # Apical sel between 0 and max age/length
        # Knife edge ascending selectivity, ascend SD = 0.1
        # More sloping descending limb of selectivity, descend SD = 2
        val[] <- c(0, log(0.1), log(2))
      }
      return(val)
    })
  }

  # Index parameters ----
  ni <- MARSdata@Dsurvey@ni
  if (ni > 0 && is.null(p$sel_pi)) {
    p$sel_pi <- sapply(1:ni, function(i) {
      sel_i <- MARSdata@Dsurvey@sel_i[i]
      val <- numeric(3)
      if (grepl("length", sel_i) && length(MARSdata@Dsurvey@IALobs_ymli)) {

        IAL <- apply(MARSdata@Dsurvey@IALobs_ymli[, , , i, drop = FALSE], 3, sum)

        if (sum(IAL)) {
          LFS <- min(MARSdata@Dmodel@lmid[which.max(IAL)], 0.75 * max(MARSdata@Dmodel@lmid))
          L5 <- approx(cumsum(IAL)/sum(IAL), MARSdata@Dmodel@lmid, 0.05)$y

          if (L5 < LFS) {
            sigma_asc <- min((LFS - L5)/sqrt(-log(0.05, 2)), 0.25 * diff(range(MARSdata@Dmodel@lmid)))
            val[2:3] <- log(sigma_asc)
            val[1] <- qlogis(LFS/max(0.95 * MARSdata@Dmodel@lmid))
          }
        }
      }
      if (all(!val)) {
        # Apical sel between 0 and max age/length
        # Knife edge ascending selectivity, ascend SD = 0.1
        # More sloping descending limb of selectivity, descend SD = 2
        val[] <- c(0, log(0.1), log(2))
      }
      return(val)
    })
  }

  # Initial conditions ----
  if (is.null(p$log_initF_mfr)) {
    p$log_initF_mfr <- ifelse(MARSdata@Dfishery@Cinit_mfr < 1e-8, -1000, log(0.1))
  }
  if (is.null(p$log_initrdev_as)) {
    p$log_initrdev_as <- matrix(0, na, ns)
  }

  do_map <- make_map(p, MARSdata, map = map, silent = silent, ...)
  out <- c(list(p = p), do_map)
  out$p <- check_parameters(out$p, out$map, MARSdata, silent)

  return(out)
}

#' @rdname make_parameters
#' @aliases make_map
#' @param p List of parameters, e.g., returned by [make_parameters()]
#' @param map List of mapped parameters. Overrides following `est_*` arguments
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
#' By default `p$mov_x_marrs` is zero. Set to -1000 for areas for which there is no abundance of a particular stock.
#'
#' @importFrom dplyr filter
#' @importFrom rlang .data .env
#' @return
#' [make_map()] returns a named list containing parameter mappings (`"map"`) and a character vector of random effects (`"random"`).
#' @export
make_map <- function(p, MARSdata, map = list(),
                     est_M = FALSE, est_h = FALSE, est_mat = FALSE, est_sdr = FALSE,
                     est_mov = c("none", "dist_random", "gravity_fixed"),
                     est_qfs = FALSE,
                     silent = FALSE) {

  est_mov <- match.arg(est_mov)
  getAllS4(MARSdata)
  getAllS4(MARSdata@Dmodel)

  nf <- MARSdata@Dfishery@nf

  random <- NULL
  #map <- list()

  # Stock parameters ----
  if (!silent) {
    R0_s <- signif(exp(p$t_R0_s) * scale_s, 4)
    if (is.null(map$t_R0_s)) {
      message_info("Estimating t_R0_s, starting R0 = ", paste(R0_s, collapse = ", "))
    } else {
      message_info("Starting R0 = ", paste(R0_s, collapse = ", "))
      if (any(is.na(map$t_R0_s))) {
        message_info("Fixed for stock ", paste(which(is.na(map$t_R0_s)), collapse = ", "))
      }
    }

    if (is.null(map$t_h_s)) {
      h_s <- sapply(1:ns, function(s) conv_steepness(p$t_h_s[s], Dstock@SRR_s[s])) %>% signif(4)
      message_info("Starting steepness = ", paste(h_s, collapse = ", "))
      if (any(is.na(map$t_h_s))) {
        message_info("Fixed for stock ", paste(which(is.na(map$t_h_s)), collapse = ", "))
      }
    }
  }

  if (is.null(map$mat_ps)) {
    if (est_mat) {
      if (!silent) message_info("Estimating maturity ogive parameters")
    } else {
      map$mat_ps <- factor(array(NA, dim(p$mat_ps)))
      if (!silent) message_info("Fixed maturity to values in data slot 'matd_yas'")
      if (!length(Dstock@matd_yas)) stop("Maturity ogive is not estimated. Need maturity at age values in the data slot 'matd_yas'.")
    }
  } else if (!silent && any(!is.na(map$mat_ps))) {

    message_info("Estimating maturity ogive parameters. Start values:")
    sapply(1:ns, function(s) {
      map_s <- matrix(map$mat_ps, 2, ns)[, s]

      if (any(!is.na(map_s))) {
        a50 <- signif(na * plogis(p$mat_ps[1, s]), 3)
        a95 <- signif(a50 + exp(p$mat_ps[2, s]), 3)
        message_info("Stock ", s, ": a50 = ", a50, ", a95 = ", a95)
      }
    })
  }

  if (is.null(map$log_M_s)) {
    if (est_M) {
      if (!silent) message_info("Estimating natural mortality for all stocks")
    } else {
      map$log_M_s <- factor(rep(NA, ns))
      if (!silent) message_info("Fixed natural mortality to values in data slot 'Md_yas'")
      if (!length(Dstock@Md_yas)) stop("Natural mortality is not estimated. Need M values in the data slot 'Md_yas'.")
    }
  } else if (!silent && any(!is.na(map$log_M_s))) {
    message_info("Estimating natural mortality for stock ", paste(which(!is.na(map$log_M_s)), collapse = ", "))
  }

  if (!silent) {
    if (is.null(map$log_rdev_ys)) {
      message_info("Estimating recruitment deviates for all years and stocks")
    } else if (all(is.na(map$log_rdev_ys))) {
      message_info("No recruitment deviates are estimated")
    } else {
      message_info("Estimating recruitment deviates for")
      sapply(1:ns, function(s) {
        map_s <- matrix(map$log_rdev_ys, ny, ns)[, s]
        message_info("Stock ", s, ": ", sum(!is.na(map_s)), " out of ", ny, " years")
      })
    }
  }

  if (is.null(map$log_sdr_s)) {
    if (est_sdr) {
      if (!silent) message_info("Estimating sigma_R (SD of recruitment deviates) for all stocks")
    } else {
      map$log_sdr_s <- factor(rep(NA, ns))
      if (!silent) message_info("Fixed sigma_R (SD of recruitment deviates) for all stocks")
    }
  } else if (!silent && any(!is.na(map$log_sdr_s))) {
    message_info("Estimating sigma_R (SD of recruitment deviates) for stock ", paste(which(!is.na(map$log_sdr_s)), collapse = ", "))
  }

  if (nr == 1 || est_mov == "none") {
    map$mov_x_marrs <- factor(array(NA, dim(p$mov_x_marrs)))
    map$mov_g_ymars <- factor(array(NA, dim(p$mov_g_ymars)))
    map$mov_v_ymas <- factor(array(NA, dim(p$mov_v_ymas)))

    map$log_sdg_rs <- factor(array(NA, dim(p$log_sdg_rs)))
    map$t_corg_ps <- factor(array(NA, dim(p$t_corg_ps)))

    map$log_recdist_rs <- factor(matrix(NA, nr, ns))

    if (!silent) message_info("No stock movement or recruitment distribution parameters are estimated")

  } else {

    if (is.null(map$mov_x_marrs)) map$mov_x_marrs <- factor(array(NA, dim(p$mov_x_marrs)))

    if (is.null(map$mov_g_ymars)) {
      map$mov_g_ymars <- array(NA, c(ny, nm, na, nr, ns))

      # Group parameters based on data stratification
      if (length(MARSdata@Dtag@tag_yy) && length(MARSdata@Dtag@tag_aa)) {
        for (s in 1:ns) { # Estimate parameters for r = 2, ..., nr (softmax transformation)
          r_eff <- which(MARSdata@Dstock@presence_rs[, s])[-1]
          if (length(r_eff)) {
            for (i in 1:nrow(MARSdata@Dtag@tag_yy)) {
              yy <- which(MARSdata@Dtag@tag_yy[i, ] > 0)
              for (j in 1:nrow(MARSdata@Dtag@tag_aa)) {
                aa <- which(MARSdata@Dtag@tag_aa[j, ] > 0)

                gind_strat <- expand.grid(m = 1:nm, r = r_eff, s = s, yc = i, ac = j)
                gind_strat$par_no <- 1:nrow(gind_strat)

                gind <- expand.grid(y = yy, m = 1:nm, a = aa, r = r_eff, s = s, yc = i, ac = j) %>%
                  merge(gind_strat) %>%
                  as.matrix()

                if (all(is.na(map$mov_g_ymars))) {
                  parmin <- 0
                } else {
                  parmin <- max(map$mov_g_ymars, na.rm = TRUE)
                }
                map$mov_g_ymars[gind[, c("y", "m", "a", "r", "s")]] <- parmin + gind[, "par_no"]
              }
            }
          }
        }
      } else {
        for (s in 1:ns) { # Estimate parameters for r = 2, ..., nr (softmax transformation)
          r_eff <- which(MARSdata@Dstock@presence_rs[, s])[-1]
          if (length(r_eff)) {
            gind <- as.matrix(expand.grid(y = 1:ny, m = 1:nm, a = 1:na, r = r_eff, s = s))
            if (all(is.na(map$mov_g_ymars))) {
              parmin <- 0
            } else {
              parmin <- max(map$mov_g_ymars, na.rm = TRUE)
            }
            map$mov_g_ymars[gind] <- parmin + seq(1, nrow(gind))
          }
        }
      }
      map$mov_g_ymars <- factor(map$mov_g_ymars)
    }

    if (est_mov == "dist_random") {
      random <- c(random, "mov_g_ymars")
      map$mov_v_ymas <- factor(array(NA, dim(p$mov_v_ymas)))
      if (!silent) message_info("Movement estimated as random effects")
    } else {

      if (is.null(map$mov_v_ymas)) {
        map$mov_v_ymas <- array(NA, dim(p$mov_v_ymas))

        # Group parameters based on data stratification
        if (length(MARSdata@Dtag@tag_yy) && length(MARSdata@Dtag@tag_aa)) {
          for (s in 1:ns) { # Estimate parameters for r = 2, ..., nr (softmax transformation)
            nr_eff <- sum(MARSdata@Dstock@presence_rs[, s])
            if (nr_eff > 1) {
              for (i in 1:nrow(MARSdata@Dtag@tag_yy)) {
                yy <- which(MARSdata@Dtag@tag_yy[i, ] > 0)
                for (j in 1:nrow(MARSdata@Dtag@tag_aa)) {
                  aa <- which(MARSdata@Dtag@tag_aa[j, ] > 0)

                  vind_strat <- expand.grid(m = 1:nm, yc = i, ac = j)
                  vind_strat$par_no <- 1:nrow(vind_strat)

                  vind <- expand.grid(y = yy, m = 1:nm, a = aa, s = s, yc = i, ac = j) %>%
                    merge(vind_strat) %>%
                    as.matrix()

                  if (all(is.na(map$mov_v_ymas))) {
                    parmin <- 0
                  } else {
                    parmin <- max(map$mov_v_ymas, na.rm = TRUE)
                  }
                  map$mov_v_ymas[vind[, c("y", "m", "a", "s")]] <- parmin + vind[, "par_no"]
                }
              }
            }
          }
        } else {
          for (s in 1:ns) {  # Estimate parameters if nr_effective > 1
            nr_eff <- sum(MARSdata@Dstock@presence_rs[, s])
            if (nr_eff > 1) {
              vind <- as.matrix(expand.grid(y = 1:ny, m = 1:nm, a = 1:na, s = s))
              if (all(is.na(map$mov_v_ymas))) {
                parmin <- 0
              } else {
                parmin <- max(map$mov_v_ymas, na.rm = TRUE)
              }
              map$mov_v_ymas[vind] <- parmin + seq(1, nrow(vind))
            }
          }
        }

        map$mov_v_ymas <- factor(map$mov_v_ymas)
      }

      if (is.null(map$log_sdg_rs)) map$log_sdg_rs <- factor(array(NA, dim(p$log_sdg_rs)))
      if (is.null(map$t_corg_ps)) map$t_corg_ps <- factor(array(NA, dim(p$t_corg_ps)))

      if (!silent) message_info("Stock movement is estimated as fixed effects")
    }
  }

  if (is.null(map$log_recdist_rs)) {
    recdist_rs <- sapply(1:ns, function(s) {
      x <- ifelse(Dmodel@presence_rs[, s], NA, TRUE)
      x[which(x)[1]] <- NA
      return(x)
    })
    recdist_rs[!is.na(recdist_rs)] <- 1:sum(recdist_rs, na.rm = TRUE)
    map$log_recdist_rs <- factor(recdist_rs)
  }

  if (!silent && nr > 1 && any(!is.na(map$log_recdist_rs))) {
    message_info("Recruitment distribution will be estimated")
  }

  # Fleet parameters ----
  if (!est_qfs || ns == 1) {
    map$log_q_fs <- factor(matrix(NA, nf, ns))
    if (!silent && ns > 1) message_info("Fishery fleet catchability equal for all stocks")
  } else { # Fix q for first stock
    map$log_q_fs <- local({
      np <- (ns - 1) * nf
      m <- matrix(NA, nf, ns)
      m[, -1] <- 1:np
      factor(m)
    })
    if (!silent) message_info("Fishery catchability to be estimated for all fleets (relative to stock 1)")
  }
  if (condition == "F" && any(Dfishery@Cobs_ymfr < 1e-8)) {
    map$log_Fdev_ymfr <- local({
      m <- ifelse(Dfishery@Cobs_ymfr < 1e-8, NA, TRUE)
      m[!is.na(m)] <- 1:sum(m, na.rm = TRUE)
      factor(m)
    })
  } else if (condition == "catch") {
    map$log_Fdev_ymfr <- factor(array(NA, c(ny, nm, nf, nr)))
  }
  if (!silent && condition == "F") {
    message_info("F is an estimated parameter for all corresponding catches greater than 1e-8")
  }

  ## Fix dome parameter if selectivity is logistic or all parameters if mirrored to maturity
  if (any(!grepl("dome", Dfishery@sel_f))) {
    nsel <- max(Dfishery@sel_block_yf)
    sel_pf <- sapply(1:nsel, function(f) {
      sel_f <- Dfishery@sel_f[f]
      vec <- rep(TRUE, 3)
      if (sel_f %in% c("logistic_age", "logistic_length")) vec[3] <- NA
      if (sel_f %in% c("B", "SB")) vec[] <- NA
      return(vec)
    })
    sel_pf[!is.na(sel_pf)] <- 1:sum(sel_pf, na.rm = TRUE)
    map$sel_pf <- factor(sel_pf)
  }

  if (!silent) {
    message_info("Fishery selectivity setup:")

    fsel_start <- conv_selpar(p$sel_pf, type = Dfishery@sel_f, maxage = Dmodel@na, maxL = 0.95 * max(MARSdata@Dmodel@lmid))
    y <- if (length(Dlabel@year)) {
      Dlabel@year
    } else {
      1:ny
    }
    no_blocks <- apply(Dfishery@sel_block_yf, 2, function(x) length(unique(x)) == 1) %>% all()
    for (bb in unique(Dfishery@sel_block_yf)) {
      if (no_blocks) {
        if (length(Dfishery@CAAobs_ymafr) > 0) {
          nage <- sum(
            apply(Dfishery@CAAobs_ymafr[, , , bb, ], 1, function(x) sum(x, na.rm = TRUE)) > 0,
            na.rm = TRUE
          )
        } else {
          nage <- 0
        }
        if (length(Dfishery@CALobs_ymlfr)) {
          nlen <- sum(
            apply(Dfishery@CALobs_ymlfr[, , , bb, ], 1, function(x) sum(x, na.rm = TRUE)) > 0,
            na.rm = TRUE
          )
        } else {
          nlen <- 0
        }
        age <- paste(nage, "years age composition")
        len <- paste(nlen, "years length composition")
        if (length(Dlabel@fleet)) {
          fname <- paste0(" (", Dlabel@fleet[bb], ")")
        } else {
          fname <- ""
        }
        if (nage || nlen) {
          message_info("Fleet ", bb, fname, ": ", Dfishery@sel_f[bb], ", ", age, " and ", len)
        } else {
          message_info("Fleet ", bb, fname, ": ", Dfishery@sel_f[bb], ", no composition data")
        }
      } else {
        fleet <- lapply(1:nf, function(ff) {
          yy <- y[Dfishery@sel_block_yf[, ff] == bb]
          if (length(yy)) {
            if (all(diff(y) == 1)) {
              paste0(ff, " (", range(yy) %>% paste(collapse = "-"), ")")
            } else {
              paste0(ff, " (", range(yy) %>% paste(collapse = "-"), ", with gaps)")
            }
          } else {
            NULL
          }
        })
        message_info("Block ", bb, " (", Dfishery@sel_f[bb], ") assigned to fleet:\n", do.call(c, fleet) %>% paste(collapse = "\n"))
      }

      if (grepl("logistic", Dfishery@sel_f[bb])) {
        message_info("   Selectivity start values: full sel = ", round(fsel_start[1, bb], 2),
                     ", ascending limb SD = ", round(fsel_start[2, bb], 2))
      }
      if (grepl("dome", Dfishery@sel_f[bb])) {
        message_info("   Selectivity start values: full sel = ", round(fsel_start[1, bb], 2),
                     ", ascending limb SD = ", round(fsel_start[2, bb], 2),
                     ", descending limb SD = ", round(fsel_start[3, bb], 2))
      }
      message_info("\n")
    }
  }

  # Survey parameters ----
  ## Fix dome parameter or all parameters or all parameters if mirrored to fleet or maturity
  old_warn <- options()$warn
  options(warn = -1)
  on.exit(options(warn = old_warn))

  if (any(!grepl("dome", Dsurvey@sel_i))) {
    sel_pi <- sapply(1:Dsurvey@ni, function(i) {
      sel_i <- Dsurvey@sel_i[i]
      vec <- rep(TRUE, 3)
      if (sel_i %in% c("logistic_age", "logistic_length")) vec[3] <- NA
      if (sel_i %in% c("B", "SB") || !is.na(as.integer(sel_i))) vec[] <- NA
      return(vec)
    })
    sel_pi[!is.na(sel_pi)] <- 1:sum(sel_pi, na.rm = TRUE)

    map$sel_pi <- factor(sel_pi)
  }

  if (!silent) {
    message_info("Index selectivity setup:")
    isel_start <- conv_selpar(p$sel_pi, type = Dsurvey@sel_i, maxage = Dmodel@na, maxL = 0.95 * max(lmid))
    for (i in 1:Dsurvey@ni) {
      sel_i <- Dsurvey@sel_i[i]
      if (length(Dlabel@index)) {
        iname <- paste0(" (", Dlabel@index[i], ")")
      } else {
        iname <- ""
      }
      r <- which(apply(Dsurvey@samp_irs[i, , , drop = FALSE], 2, sum) > 0)
      rtext <- paste("region", paste(r, collapse = ", "))
      s <- which(apply(Dsurvey@samp_irs[i, , , drop = FALSE], 3, sum) > 0)
      stext <- paste("stock", paste(s, collapse = ", "))

      if (!is.na(as.integer(sel_i))) {
        message_info("Index ", i, iname, ": fleet ", sel_i, " in ", rtext, "; ", stext)
      } else {
        s <-
        message_info("Index ", i, iname, ": ", sel_i, " in ", rtext, "; ", stext)
      }

      if (grepl("logistic", sel_i)) {
        message_info("   Selectivity start values: full sel = ", round(isel_start[1, i], 2),
                     ", ascending limb SD = ", round(isel_start[2, i], 2))
      }
      if (grepl("dome", sel_i)) {
        message_info("   Selectivity start values: full sel = ", round(isel_start[1, i], 2),
                     ", ascending limb SD = ", round(isel_start[2, i], 2),
                     ", descending limb SD = ", round(isel_start[3, i], 2))
      }
    }
  }

  # Initial conditions ----
  if (any(Dfishery@Cinit_mfr < 1e-8)) {
    map$log_initF_mfr <- local({
      m <- ifelse(Dfishery@Cinit_mfr < 1e-8, NA, TRUE)
      m[!is.na(m)] <- 1:sum(m, na.rm = TRUE)
      factor(m)
    })
  }
  if (!silent && (is.null(map$log_initF_mfr) || any(!is.na(map$log_initF_mfr)))) {
    message_info("Initial equilibrium F will be estimated")
  }

  if (is.null(map$log_initrdev_as)) map$log_initrdev_as <- factor(matrix(NA, na, ns))
  if (!silent && (is.null(map$log_initrdev_as) || any(!is.na(map$log_initrdev_as)))) {
    message_info("Year 1 recruitment deviations will be estimated")
  }

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
        data.frame(Parameter = par_df[x], Number = nparam)
      } else {
        NULL
      }
    })
    output <- do.call(rbind, npar)

    message_info("Total estimated parameters: ", sum(output[, "Number"]))
    print(do.call(rbind, npar))
  }

  return(invisible(p))
}

par_df = c(
  "t_R0_s" = "Unfished recruitment",
  "t_h_s" = "Steepness",
  "mat_ps" = "Maturity ogive",
  "log_M_s" = "Natural mortality",
  "log_rdev_ys" = "Recruitment deviations",
  "log_sdr_s" = "Recruitment deviation standard deviation",
  "mov_x_marrs" = "Base movement",
  "mov_g_ymars" = "Movement, gravity to regions",
  "mov_v_ymas" = "Movement, region viscosity/retention",
  "log_sdg_rs" = "Region distribution standard deviation",
  "t_corg_ps" = "Region distribution correlation matrix",
  "log_recdist_rs" = "Region distribution of recruitment",
  "log_q_fs" = "Relative catchability of stocks by fleet",
  "log_Fdev_ymfr" = "F deviations",
  "sel_pf" = "Fishery selectivity",
  "sel_pi" = "Index selectivity",
  "log_initF_mfr" = "Equilibrium (year 1) F",
  "log_initrdev_as" = "Initial (year 1) recruitment deviations"
)
