

#' Fit MARS model
#'
#' Wrapper function that calls RTMB to create the model and perform the numerical optimization
#'
#' @param MARSdata Data object. Class \linkS4class{MARSdata}, validated by [check_data()]
#' @param parameters List of parameters, validated by [check_parameters]
#' @param map List of parameters indicated whether they are fixed and/or how they are shared. See [TMB::MakeADFun()].
#' @param random Character vector indicating the parameters that are random effects.
#' @param run_model Logical, indicates whether the model will be fitted through [stats::nlminb()].
#' @param do_sd Logical, indicates whether the standard deviations of parameters will be calculated with [TMB::sdreport()].
#' @param silent Logical, passed to `MakeADFun`
#' @param control Passed to [stats::nlminb()]
#' @param ... Other arguments to [TMB::MakeADFun()].
#' @importFrom methods new
#' @export
fit_MARS <- function(MARSdata, parameters, map = list(), random = NULL,
                 run_model = TRUE, do_sd = TRUE, silent = TRUE,
                 control = list(iter.max = 2e+05, eval.max = 4e+05), ...) {

  MARSdata@Misc$map <- map
  MARSdata@Misc$random <- random

  RTMB::TapeConfig(comparison = "tape")
  func <- function(p) .MARS(p, d = MARSdata)

  obj <- MakeADFun(
    func = func, parameters = parameters,
    map = map, random = random,
    silent = silent,
    ...
  )

  if (run_model) m <- optimize_RTMB_model(obj, do_sd = do_sd, control = control)

  M <- new(
    "MARSassess",
    obj = obj,
    report = obj$report(obj$env$last.par.best) %>% update_report(MARSdata)
  )

  if (run_model) {
    M@opt <- m$opt
    if (do_sd) M@SD <- m$SD
  }
  return(M)
}

update_report <- function(r, MARSdata) {

  if (is.null(r$F_yas)) {
    ny <- MARSdata@Dmodel@ny
    na <- MARSdata@Dmodel@na
    ns <- MARSdata@Dmodel@ns
    nf <- MARSdata@Dmodel@nf
    Fmax <- MARSdata@Dmodel@Fmax

    r$F_yas <- sapply2(1:ns, function(s) {
      sapply(1:na, function(a) {
        sapply(1:ny, function(y) {
          calc_summary_F(M = r$M_yas[y, a, s], N = apply(r$N_ymars[y, 1, a, , s, drop = FALSE], 3, sum),
                         CN = apply(r$CN_ymafrs[y, , a, , , s, drop = FALSE], 3, sum), Fmax = nf * Fmax)
        })
      })
    })
    r$Z_yas <- r$F_yas + r$M_yas
  }
  return(r)
}

.MARS <- function(p = list(), d) {

  # Assign data variables to environment, see OBS() for simulation ----
  getAllS4(d@Dmodel, d@Dstock, d@Dfishery, d@Dsurvey, d@DCKMR, d@Dtag)
  map <- d@Misc$map
  random <- d@Misc$random

  # Transform data ----
  delta_m <- 1/nm

  # Population arrays ----
  N_ymars <- array(0, c(ny + 1, nm, na, nr, ns))
  Nsp_yars <- array(NA_real_, c(ny, na, ns))
  SB_yrs <- array(NA_real_, c(ny, nr, ns))
  R_ys <-
    Rdev_ys <- array(NA_real_, c(ny, ns))
  FM_ymars <-
    Z_ymars <- array(NA_real_, c(ny, nm, na, nr, ns))
  if (length(HSP_s)) F_yas <- array(NA_real_, c(ny, na, ns))
  #LAK_ymals <- array(NA_real_, c(ny, nm, na, nl, ns))
  mov_ymarrs <- array(NA_real_, c(ny, nm, na, nr, nr, ns))

  # Fishery arrays ----
  sel_ymafs <- array(NA_real_, c(ny, nm, na, nf, ns))

  FM_ymafrs <-
    CN_ymafrs <- array(NA_real_, c(ny, nm, na, nf, nr, ns))
  if (any(CALobs_ymlfr > 0, na.rm = TRUE)) CN_ymlfrs <- array(NA_real_, c(ny, nm, nl, nf, nr, ns))
  CB_ymfrs <-
    VB_ymfrs <- array(NA_real_, c(ny, nm, nf, nr, ns))

  # Index of abundance arrays ----
  if (ni > 0) {
    IN_ymais <-
      sel_ymais <- array(NA_real_, c(ny, nm, na, ni, ns))
    if (any(IALobs_ymli > 0, na.rm = TRUE)) IN_ymlis <- array(NA_real_, c(ny, nm, na, nl, ns))
    VI_ymi <- I_ymi <- array(NA_real_, c(ny, nm, ni))
  }

  # Transform parameters ----
  ## Maturity at age ogive ----
  mat_yas <- sapply2(1:ns, function(s) {
    if (all(is.na(map$log_M_ps[, s]))) {
      matd_yas[, , s]
    } else {
      a50 <- na * plogis(p$mat_ps[1, s])
      a95 <- a50 + exp(p$mat_ps[2, s])

      a <- seq(1, na)
      m <- 1/(1 + exp(-log(19) * (a - a50)/(a95 - a50)))
      matrix(m, ny, na, byrow = TRUE)
    }
  })

  ## Natural mortality ----
  M_yas <- sapply2(1:ns, function(s) {
    if (is.na(map$log_M_s[s])) {
      Md_yas[, , s]
    } else {
      matrix(exp(p$log_M_s[s]), ny, na)
    }
  })

  ## Stock recruit parameters ----
  R0_s <- exp(p$t_R0_s + scale_s)
  h_s <- ifelse(SRR_s == "BH", 0.8 * plogis(p$t_h_s), exp(p$t_h_s)) + 0.2
  kappa_s <- sapply(1:ns, function(s) SRkconv(h_s[s], SRR = SRR_s[s]))

  phi_s <- sapply(1:ns, function(s) {
    calc_phi(M_yas[y_phi, , s], mat_yas[y_phi, , s] * fec_yas[y_phi, , s],
             delta = (m_spawn - 1 + delta_s[s]) * delta_m)
  })
  SB0_s <- R0_s * phi_s
  sralpha_s <- kappa_s/phi_s
  srbeta_s <- sapply(1:ns, function(s) SRbetaconv(h_s[s], R0_s[s], phi_s[s], SRR = SRR_s[s]))

  ## Recruitment deviates ----
  Rdev_ys[] <- exp(p$log_rdev_ys)

  ## Fishery selectivity ----
  selconv_pf <- conv_selpar(p$sel_pf, type = sel_f, maxage = na, maxL = 0.95 * max(lmid))
  sel_lf <- calc_sel_len(selconv_pf, lmid, type = sel_f)
  q_fs <- exp(p$log_q_fs)

  ## Index selectivity ----
  if (ni > 0) {
    selconv_pi <- conv_selpar(p$sel_pi, type = sel_i, maxage = na, maxL = 0.95 * max(lmid))
    sel_li <- calc_sel_len(selconv_pi, lmid, type = sel_i)
  }

  ## Stock distribution and movement parameters ----
  #dist_ymars[y, m, , , ] <- sapply2(1:ns, function(s) {
  #  sapply(1:na, function(a) softmax(p$mov_g_ymars[y, m, a, , s])) %>% t()
  #})
  for(y in 1:ny) {
    for(m in 1:nm) {
      mov_ymarrs[y, m, , , , ] <- sapply2(1:ns, function(s) {
        conv_mov(p$mov_x_ymarrs[y, m, , , , s], p$mov_g_ymars[y, m, , , s], p$mov_v_ymas[y, m, , s], na, nr)
      })
    }
  }

  # Miscellaneous penalty term, e.g., F > Fmax
  penalty <- 0

  # First year, first season initialization ----

  # Loops over years and seasons ----
  for(y in 1:ny) {
    for(m in 1:nm) {
      ## Calculate length-age key and fishery and index age selectivity ----
      #LAK_ymals[y, m, , , ] <- sapply2(1:ns, function(s) calc_LAK(len_ymas[y, m, , s], sdlen_ymas[y, m, , s], lbin))
      sel_ymafs[y, m, , , ] <- sapply2(1:ns, function(s) {
        calc_fsel_age(sel_lf, LAK_ymals[y, m, , , s], sel_f, selconv_pf, sel_block_yf[y, ],
                      mat = mat_yas[y, , s], a = seq(1, na))
      })
      if (ni > 0) {
        sel_ymais[y, m, , , ] <- sapply2(1:ns, function(s) {
          calc_isel_age(sel_li, LAK_ymals[y, m, , , s], sel_i, selconv_pi, sel_ymafs[y, m, , , s],
                        mat = mat_yas[y, , s], a = seq(1, na))
        })
      }

      ## This season's mortality ----
      Fsearch <- calc_F(
        Cobs = Cobs_ymfr[y, m, , ], N = N_ymars[y, m, , , ], sel = sel_ymafs[y, m, , , ],
        wt = fwt_ymafs[y, m, , , ], M = M_yas[y, , ], q_fs = q_fs, delta = delta_m,
        na = na, nr = nr, nf = nf, ns = ns, Fmax = Fmax, nitF = nitF, trans = "log"
      )
      penalty <- penalty + Fsearch[["penalty"]] # add penalty for exceeding Fmax

      FM_ymars[y, m, , , ] <- Fsearch[["F_ars"]]
      Z_ymars[y, m, , , ] <- Fsearch[["Z_ars"]]

      ## This season's fishery catch, vulnerable biomass, and total biomass ----
      CN_ymafrs[y, m, , , ] <- Fsearch[["CN_afrs"]]

      if (any(CALobs_ymlfr > 0, na.rm = TRUE)) { # If there's any length data
        CN_ymlfrs[y, m, , , , ] <- sapply2(1:ns, function(s) {
          sapply2(1:nr, function(r) {
            sapply(1:nf, function(f) CN_ymafrs[y, m, , f, r, s] %*% LAK_ymals[y, m, , , s])
          })
        })
      }

      CB_ymfrs[y, m, , , ] <- Fsearch[["CB_frs"]]
      VB_ymfrs[y, m, , ] <- apply(Fsearch[["VB_afrs"]], 2:4, sum)
      B_ymrs[y, m, , ] <- sapply(1:nr, function(r) rowSums(N_ymars[y, m, , r, ] * swt_ymas[y, m, , ])) %>% t()

      ## This year's spawning and recruitment ----
      if (m == m_spawn) {
        Nsp_yars[y, , , ] <- sapply2(1:ns, function(s) {
          sapply(1:nr, function(r) {
            natal_rs[r, s] * N_ymars[y, m, , r, s] * exp(-delta_s[s] * Z_ymars[y, m, , r, s]) * mat_yas[y, , s]
          })
        })
        SB_yrs[y, , ] <- sapply(1:ns, function(s) {
          sapply(1:nr, function(r) sum(Nsp_yars[y, , r, s] * fec_yas[y, , s]))
        })
        y_rec <- ifelse(m_rec > m_spawn, y, y+1) # When nm = 1, this only works if recruitment is age 1
        R_ys[y_rec, ] <- Rdev_ys[y_rec, ] * sapply(1:ns, function(s) {
          calc_recruitment(SB = sum(SB_yrs[y, , s]), SRR = SRR_s[s], a = sralpha_s[s], b = srbeta_s[s])
        })

        # Distribution of recruits is proportional to regional abundance of spawners
        mov_ymarrs[y_rec, m_rec, 1, , , ] <- sapply2(1:ns, function(s) {
          matrix(SB_yrs[y, , s]/sum(SB_yrs[y, , s]), nr, nr, byrow = TRUE)
        })
      }

      ## This year's index ----
      if (ni > 0) {
        IN_ymais[y, m, , , ] <- calc_index(
          N = N_ymars[y, m, , , ], Z = Z_ymars[y, m, , , ], sel = sel_ymais[y, m, , , ], samp = samp_irs,
          delta = delta_i, na = na, ns = ns, ni = ni
        )

        VI_ymi[y, m, ] <- sapply(1:ni, function(i) {
          I_s <- sapply(1:ns, function(s) {
            w <- if (unit_i[i] == "N") 1 else swt_ymas[y, m, , s]
            sum(IN_ymais[y, m, , i, s] * w)
          })
          sum(I_s)
        })

        if (any(IALobs_ymli > 0, na.rm = TRUE)) { # If there's any length data
          IN_ymlis <- sapply2(1:ns, function(s) {
            sapply(1:ni, function(i) IN_ymais[y, m, , i, s] %*% LAK_ymals[y, m, , , s])
          })
        }
      }

      ## Next season's abundance and total biomass ----
      # Have to advance the age classes in the season before spawning
      ynext <- ifelse(m == nm, y+1, y)
      mnext <- ifelse(m == nm, 1, m+1)
      N_ymars[ynext, mnext, , , ] <- calc_nextN(
        N = N_ymars[y, m, , , ], surv = exp(-Z_ymars[y, m, , , ]),
        na = na, nr = nr, ns = ns,
        advance_age = mnext == m_rec,
        R = R_ys[ynext, ],
        mov = mov_ymarrs[ynext, mnext, , , , ]
      )
    }

    ## Summary F and Z by year ----
    if (length(HSP_s)) {
      F_yas[y, , ] <- sapply(1:ns, function(s) {
        sapply(1:na, function(a) {
          calc_summary_F(M = M_yas[y, a, s], N = apply(N_ymars[y, 1, a, , s, drop = FALSE], 3, sum),
                         CN = apply(CN_ymafrs[y, , a, , , s, drop = FALSE], 3, sum), Fmax = nf * Fmax)
        })
      })
    }
  }

  # Likelihoods ----
  ## Marginal fishery age composition ----
  if (any(CAAobs_ymafr > 0, na.rm = TRUE)) {
    CN_ymafr <- apply(CN_ymafrs, 1:5, sum)

    loglike_CAA_ymfr <- sapply2(1:nr, function(r) {
      sapply2(1:nf, function(f) {
        sapply(1:nm, function(m) {
          sapply(1:ny, function(y) {
            pred <- CN_ymafr[y, m, , f, r]
            like_comp(obs = CAAobs_ymafr[y, m, , f, r], pred = pred, type = fcomp_like,
                      N = CAAN_ymfr[y, m, f, r], theta = CAAtheta_f[f],
                      stdev = sqrt(sum(pred)/pred))
          })
        })
      })
    })
  } else {
    loglike_CAA_ymfr <- 0
  }

  ## Marginal fishery length composition ----
  if (any(CALobs_ymlfr > 0, na.rm = TRUE)) { # If there's any length data
    CN_ymlfr <- apply(CN_ymlfrs, 1:5, sum)

    loglike_CAL_ymfr <- sapply2(1:nr, function(r) {
      sapply2(1:nf, function(f) {
        sapply(1:nm, function(m) {
          sapply(1:ny, function(y) {
            pred <- CN_ymlfr[y, m, , f, r]
            like_comp(obs = CALobs_ymlfr[y, m, , f, r], pred = pred, type = fcomp_like,
                      N = CALN_ymfr[y, m, f, r], theta = CALtheta_f[f],
                      stdev = sqrt(sum(pred)/pred))
          })
        })
      })
    })
  } else {
    loglike_CAL_ymfr <- 0
  }

  ## Index ----
  if (ni > 0) {
    q_i <- sapply(1:ni, function(i) calc_q(Iobs_ymi[, , i], B = VI_ymi[, , i]))
    I_ymi <- sapply2(1:ni, function(i) q_i[i] * VI_ymi[, , i])
    loglike_I_ymi <- ifelse(is.na(Iobs_ymi), 0, dnorm(log(Iobs_ymi/I_ymi), 0, Isd_ymi, log = TRUE))
  } else {
    loglike_I_ymi <- 0
  }

  if (ni > 0 && any(IAAobs_ymai > 0, na.rm = TRUE)) {
    IN_ymai <- apply(IN_ymais, 1:4, sum)

    loglike_IAA_ymi <- sapply2(1:ni, function(i) {
      sapply(1:nm, function(m) {
        sapply(1:ny, function(y) {
          pred <- IN_ymai[y, m, , i]
          like_comp(obs = IAAobs_ymai[y, m, , i], pred = pred, type = icomp_like,
                    N = IAAN_ymi[y, m, i], theta = IAAtheta_i[i],
                    stdev = sqrt(sum(pred)/pred))
        })
      })
    })
  } else {
    loglike_IAA_ymi <- 0
  }

  if (ni > 0 && any(IALobs_ymli > 0, na.rm = TRUE)) {
    IN_ymli <- apply(IN_ymlis, 1:4, sum)

    loglike_IAL_ymi <- sapply2(1:ni, function(i) {
      sapply(1:nm, function(m) {
        sapply(1:ny, function(y) {
          pred <- IN_ymli[y, m, , i]
          like_comp(obs = IALobs_ymli[y, m, , i], pred = pred, type = icomp_like,
                    N = IALN_ymi[y, m, i], theta = IALtheta_i[i],
                    stdev = sqrt(sum(pred)/pred))
        })
      })
    })
  } else {
    loglike_IAL_ymi <- 0
  }

  ## Stock composition ----
  if (ns > 1 && length(SC_ymafrs)) {

    loglike_SC_ymafr <- sapply2(1:nr, function(r) {
      sapply2(1:nrow(SC_ff), function(ff) { # Aggregate over fleets SC_ff
        sapply2(1:nrow(SC_aa), function(aa) { # Aggregate over age classes SC_aa
          sapply(1:nm, function(m) {
            sapply(1:ny, function(y) {
              pred <- apply(CN_ymafrs[y, m, SC_aa[aa, ], SC_ff[ff, ], r, , drop = FALSE], 6, sum)
              like_comp(obs = SC_ymafrs[y, m, aa, ff, r, ], pred = pred, type = SC_like,
                        N = SCN_ymafr[y, m, aa, ff, r], theta = SCtheta_f[ff],
                        stdev = SCstdev_f[ff])
            })
          })
        })
      })
    })

  } else {
    loglike_SC_ymafr <- 0
  }

  ## CKMR ----
  if (length(POP_s)) {
    pPOP_s <- lapply(1:ns, function(s) {
      calc_POP(t = POP_s[[s]]$t, a = POP_s[[s]]$a, y = POP_s[[s]]$y,
               N = apply(Nsp_yars[, , , s], 1:2, sum), fec = fec_yas[, , s])
    })
    loglike_POP_s <- lapply(1:ns, function(s) like_CKMR(n = POP_s[[s]]$n, m = POP_s[[s]]$m, p = pPOP_s[[s]], type = CKMR_like))
  } else {
    loglike_POP_s <- 0
  }

  if (length(HSP_s)) {
    Z_yas <- F_yas + M_yas
    pHSP_s <- lapply(1:ns, function(s) {
      calc_HSP(yi = HSP_s[[s]]$yi, yj = HSP_s[[s]]$yj,
               N = apply(Nsp_yars[, , , s], 1:2, sum), fec = fec_yas[, , s], Z = Z_yas[, , s])
    })
    loglike_HSP_s <- lapply(1:ns, function(s) like_CKMR(n = HSP_s[[s]]$n, m = HSP_s[[s]]$m, p = pHSP_s[[s]], type = CKMR_like))
  } else {
    loglike_HSP_s <- 0
  }

  loglike_tag <- 0

  loglike <- sum(loglike_I_ymi) + sum(loglike_IAA_ymi) + sum(loglike_IAL_ymi) +
    sum(loglike_CAA_ymfr) + sum(loglike_CAL_ymfr) + sum(loglike_SC_ymafr) +
    Reduce(sum, loglike_POP_s) + Reduce(sum, loglike_HSP_s) + sum(loglike_tag)

  # Priors ----
  sdr_s <- exp(p$log_sdr_s)
  bcr_s <- -0.5 * sdr_s * sdr_s

  logprior_rdev_ys <- sapply(1:ns, function(s) {
    ifelse(
      is.na(map$log_rdev_ys[, s]) | duplicated(map$log_rdev_ys[, s]),
      0,
      dnorm(p$log_rdev_ys[, s], bcr_s[s], sdr_s[s], log = TRUE)
    )
  })

  if (nr > 1 && "mov_g_ymars" %in% random) {
    sdg_s <- sapply2(1:ns, function(s) conv_Sigma(sigma = exp(p$log_sdg_rs[, s]), lower_diag = p$t_rhog_rs[, s]))

    logprior_dist_ymas <- sapply2(1:ns, function(s) {
      sapply2(1:na, function(a) {
        sapply(1:nm, function(m) {
          sapply(1:ny, function(y) dmvnorm(p$mov_g_ymars[y, m, a, , s], mu = 0, Sigma = sdg_s[, , s], log = TRUE))
        })
      })
    })
  } else {
    logprior_dist_ymas <- 0
  }

  logprior <- sum(logprior_rdev_ys) + sum(logprior_dist_ymas)

  # Objective function ----
  fn <- -1 * (logprior + loglike) + penalty

  # Report out variables ----
  ## Parameters ----
  ADREPORT(R0_s)
  ADREPORT(h_s)
  REPORT(kappa_s)
  REPORT(SB0_s)
  REPORT(sralpha_s)
  REPORT(srbeta_s)
  REPORT(sdr_s)

  REPORT(selconv_pf)
  REPORT(sel_lf)
  REPORT(q_fs)

  if (ni > 0) {
    REPORT(selconv_pi)
    REPORT(sel_li)

    ADREPORT(q_i)
  }

  REPORT(M_yas)
  REPORT(mat_yas)

  ## Population arrays ----
  REPORT(N_ymars)
  REPORT(SB_yrs)
  REPORT(R_ys)
  REPORT(Rdev_ys)
  REPORT(FM_ymars)
  REPORT(Z_ymars)
  REPORT(B_ymrs)
  if (nr > 1) REPORT(mov_ymarrs)

  ## Fishery arrays ----
  REPORT(sel_ymafs)
  REPORT(FM_ymafrs)
  REPORT(CN_ymafrs)
  if (any(CALobs_ymlfr > 0, na.rm = TRUE)) REPORT(CN_ymlfrs)
  REPORT(CB_ymfrs)
  #REPORT(VB_ymfrs)

  ## Index of abundance arrays ----
  if (ni > 0) {
    REPORT(sel_ymais)
    REPORT(I_ymi)
    if (any(IAAobs_ymai > 0, na.rm = TRUE)) REPORT(IN_ymais)
    if (any(IALobs_ymli > 0, na.rm = TRUE)) REPORT(IN_ymlis)
    REPORT(q_i)
  }

  ## CKMR ----
  if (length(POP_s)) REPORT(pPOP_s)

  if (length(HSP_s)) {
    REPORT(F_yas)
    REPORT(Z_yas)
    REPORT(pHSP_s)
  }

  ## Objective function values ----
  REPORT(loglike)
  REPORT(logprior)
  REPORT(penalty)
  REPORT(fn)

  REPORT(loglike_I_ymi)
  REPORT(loglike_IAA_ymi)
  REPORT(loglike_IAL_ymi)

  REPORT(loglike_CAA_ymfr)
  REPORT(loglike_CAL_ymfr)

  REPORT(loglike_SC_ymafr)

  REPORT(loglike_POP_s)
  REPORT(loglike_HSP_s)

  REPORT(loglike_tag)

  REPORT(logprior_rdev_ys)
  REPORT(logprior_dist_ymas)

  return(fn)
}

#' @importFrom methods slotNames
getAllS4 <- function (..., warn = TRUE) {
  fr <- parent.frame()
  dots <- list(...)

  for(i in 1:length(dots)) {
    nm <- slotNames(dots[[i]])
    for (j in nm) {
      if (warn && !is.null(fr[[j]])) warning("Object '", j, "' already defined")
      fr[[j]] <- slot(dots[[i]], j)
    }
  }
  invisible(NULL)
}
