

#' Fit MARS model
#'
#' Wrapper function that calls RTMB to create the model and perform the numerical optimization
#'
#' @param data List of data inputs, validated by \link{check_data}
#' @param parameters List of parameters, validated by \link{check_parameters}
#' @param map List of parameters indicated whether they are fixed and/or how they are shared. See \link[TMB]{MakeADFun}.
#' @param random Character vector indicating the parameters that are random effects.
#' @param run_model Logical, indicates whether the model will be fitted through \link[stats]{nlminb}.
#' @param do_sd Logical, indicates whether the standard deviations of parameters will be calculated with \link[TMB]{sdreport}.
#' @param ... Other arguments to \link[RTMB]{MakeADFun}.
#' @export
MARS <- function(data, parameters, map = list(), random = NULL, silent = TRUE, control = list(),
                 run_model = TRUE, do_sd = TRUE, ...) {

  data[["map"]] <- map
  data[["random"]] <- random

  RTMB::TapeConfig(comparison = "tape")
  func <- function(p) .MARS(p, d = data)

  obj <- MakeADFun(
    func = func, parameters = parameters,
    map = map, random = random,
    silent = silent,
    ...
  )

  if (run_model) {
    m <- optimize_RTMB_model(obj, do_sd = do_sd, control = control)
    m$report <- obj$report(obj$env$last.par.best)
  } else {
    m <- list()
  }
  m$obj <- obj

  return(m)
}

.MARS <- function(p = list(), d = list()) {

  # Assign data variables to environment, see OBS() for simulation ----
  getAll(d)

  # Transform data ----
  delta_m <- 1/nm

  phi_s <- sapply(1:ns, function(s) {
    calc_phi()
  })

  # Population arrays ----
  N_ymars <- array(NA_real_, c(ny + 1, nm, na, nr, ns))
  Nsp_yars <- array(NA_real_, c(ny, na, nr, ns))
  SB_ys <-
    R_ys <-
    Rdev_ys <- array(NA_real_, c(ny, ns))
  FM_ymars <-
    Z_ymars <-
    dist_ymars <- array(NA_real_, c(ny, nm, na, nr, ns))
  if (nrow(HSP)) F_yas <- array(NA_real_, c(ny, na, ns))
  LAK_ymals <- array(NA_real_, c(ny, nm, na, nl, ns))
  mov_ymarrs <- array(NA_real_, c(ny, nm, na, nr, nr, ns))

  # Fishery arrays ----
  fsel_ymafs <- array(NA_real_, c(ny, nm, na, nf, ns))

  FM_ymafrs <-
    CN_ymafrs <- array(NA_real_, c(ny, nm, na, nf, nr, ns))
  if (any(CALobs_ymlfs > 0, na.rm = TRUE)) CN_ymlfrs <- array(NA_real_, c(ny, nm, nl, nf, nr, ns))
  CB_ymfrs <-
    VB_ymfrs <- array(NA_real_, c(ny, nm, nf, nr, ns))

  # Index of abundance arrays ----
  if (ni > 0) {
    IN_ymais <-
      isel_ymais <- array(NA_real_, c(ny, nm, na, ni, ns))
    if (any(IALobs_ymli > 0, na.rm = TRUE)) IN_ymlis <- array(NA_real_, c(ny, nm, na, nl, ns))
    VI_ymi <- I_ymi <- array(NA_real_, c(ny, nm, ni))
  }

  # Transform parameters ----
  R0_s <- exp(p$R0x)/scale_s
  kappa_s <- exp(p$k_s) + 1

  Rdev_ys[] <- exp(p$log_rdev_ys)

  ## Fishery selectivity ----
  fsel_val <- conv_selpar(p$fsel, type = fsel_type, maxage = na - 1, Lmax = 0.95 * max(len_ymas))
  fsel_len <- calc_sel_len(fsel_val, lmid, type = fsel_type)

  ## Index selectivity ----
  if (ni > 0) {
    isel_val <- conv_selpar(p$isel, type = isel_type, maxage = na - 1, Lmax = 0.95 * max(len_ymas))
    isel_len <- calc_sel_len(isel_val, lmid, type = isel_type)
  }

  # Miscellaneous penalty term, e.g., F > Fmax
  penalty <- 0

  # Stock recruit parameters ----
  h_s <- sapply(1:ns, function(s) SRhconv(kappa_s[s], SRR = SRR_s[s]))
  B0_s <- R0_s * phi_s
  sralpha_s <- kappa_s/phi_s
  srbeta_s <- sapply(1:ns, function(s) SRbetaconv(h_s[s], R0_s[s], phi_s[s], SRR = SRR_s[s]))

  # First year, first season initialization ----

  # Loops over years and seasons ----
  for(y in 1:ny) {
    for(m in 1:nm) {
      ## Stock distribution and movement parameters ----
      dist_ymars[y, m, , , ] <- sapply2(1:ns, function(s) {
        sapply(1:na, function(a) softmax(p$dist[y, m, a, , s])) %>% t()
      })
      mov_ymarrs[y, m, , , , ] <- sapply2(1:ns, function(s) {
        conv_mov(p$mov_x[y, m, , , , s], p$mov_g[y, m, , s], p$mov_v[y, m, , s], nr, na)
      })

      ## Calculate length-age key and fishery and index age selectivity ----
      LAK_ymals[y, m, , , ] <- sapply2(1:ns, function(s) calc_LAK(len_ymas[y, m, , s], sdlen_ymas[y, m, , s], lbin))
      fsel_ymafs[y, m, , , ] <- sapply2(1:ns, function(s) {
        calc_fsel_age(fsel_len, LAK_ymals[y, m, , , s], fsel_type, fsel_val, sel_block[y, ], na - 1)
      })
      if (ni > 0) {
        isel_ymais[y, m, , , ] <- sapply2(1:ns, function(s) {
          calc_isel_age(isel_len, LAK_ymals[y, m, , , s], isel_type, isel_val, fsel_ymafs[y, m, , , s], na - 1)
        })
      }

      ## This season's mortality ----
      Fsearch <- calc_F(
        Cobs = Cobs_ymfr[y, m, , ], N = N_ymars[y, m, , , ], sel = fsel_ymafs[y, m, , , ],
        wt = fwt_yafs[y, , , ], M = M_yas[y, , ], q_fs = q_fs, delta = delta_m,
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
        Nsp_yars[y, , , ] <- sapply(1:ns, function(s) {
          Nsp <- N_ymars[y, m, , r_natal[, s], s, drop = FALSE] *
            exp(-spawn_time_frac * Z_ymars[y, m, , r_natal[, s], s, drop = FALSE])
          apply(Nsp, 3:4, sum)
        })
        SB_ys[y, ] <- sapply(1:ns, function(s) sum(Nsp_yars[y, , , s] * fec_yas[y, , , s]))
        y_rec <- ifelse(m_rec > m_spawn, y, y+1) # When nm = 1, this only works if recruitment is age 1
        R_ys[y_rec, ] <- Rdev_ys[y_rec, ] * sapply(1:ns, function(s) {
          calc_recruitment(SB = SB[y, s], SRR = SRR_s[s], a = sralpha_s[s], b = srbeta_s[s])
        })
      }

      ## This year's index ----
      if (ni > 0) {
        IN_ymais[y, m, , , ] <- calc_index(
          N = N_ymars[y, m, , , ], Z = Z_ymars[y, m, , , ], sel = isel_ymais[y, m, , , ], samp = isamp_irs,
          delta = delta_i, na = na, ns = ns, ni = ni
        )

        VI_ymi[y, m, ] <- sapply(1:ni, function(i) {
          I_s <- sapply(1:ns, function(s) {
            w <- if (iunit == "N") 1 else swt_ymas[y, m, , s]
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
      if (m == nm) {
        N_ymars[y+1, 1, , , ] <- calc_nextN(
          N = N_ymars[y, nm, , , ], surv = exp(-Z_ymars[y, nm, , , ]),
          na = na, nr = nr, ns = ns,
          advance_age = m == m_rec, type = dist_type, R = R_ys[y+1, ],
          dist = dist_ymars[y+1, 1, , , ], mov = mov_ymarrs[y+1, 1, , , , ]
        )
      } else {
        N_ymars[y, m+1, , , ] <- calc_nextN(
          N = N_ymars[y, m, , , ], surv = exp(-Z_ymars[y, m, , , ]),
          na = na, nr = nr, ns = ns,
          advance_age = m == m_rec, type = dist_type, R = R_ys[y, ],
          dist = dist_ymars[y, m+1, , , ], mov = mov_ymarrs[y, m+1, , , , ]
        )
      }
    }

    ## Summary F and Z by year ----
    if (length(HSP)) {
      F_yas[y, , ] <- sapply(1:ns, function(s) {
        sapply(1:na, function(a) {
          calc_summary_F(M = M[y, a, s], N = apply(N_ymars[y, 1, a, , s, drop = FALSE], 3, sum),
                         CN = apply(CN_ymafrs[y, , a, , , s, drop = FALSE], 3, sum), Fmax = nf * Fmax)
        })
      })
    }
  }

  # Likelihoods ----
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
          like_comp(obs = IAAobs_ymai[y, m, , i], pred = IN_ymai[y, m, , i], type = comp_like, N = , theta = )
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
          like_comp(obs = IALobs_ymli[y, m, , i], pred = IN_ymli[y, m, , i], type = comp_like, N = , theta = )
        })
      })
    })
  } else {
    loglike_IAL_ymi <- 0
  }

  ## Marginal fishery age composition ----
  if (any(CAAobs_ymafr > 0, na.rm = TRUE)) {
    CN_ymafr <- apply(CN_ymafrs, 1:5, sum)

    loglike_CAA_ymfr <- sapply2(1:nr, function(r) {
      sapply2(1:nf, function(f) {
        sapply(1:nm, function(m) {
          sapply(1:ny, function(y) {
            like_comp(obs = CAAobs_ymafr[y, m, , f, r], pred = CN_ymafr[y, m, , f, r], type = comp_like, N = , theta = )
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
            like_comp(obs = CALobs_ymlfr[y, m, , f, r], pred = CN_ymlfr[y, m, , f, r], type = comp_like, N = , theta = )
          })
        })
      })
    })
  } else {
    loglike_CAL_ymfr <- 0
  }

  ## Stock composition ----
  if (ns > 1 && length(d@SC_ymafrs)) {

    loglike_SC_ymafr <- sapply2(1:nr, function(r) {
      sapply2(1:nf, function(f) {
        sapply2(1:length(SC_a), function(aa) { # Aggegrate over age classes SC_a
          sapply(1:nm, function(m) {
            sapply(1:ny, function(y) {
              like_comp(obs = SC_ymafrs[y, m, aa, f, r, ], pred = colSums(CN_ymafrs[y, m, SC_a[[aa]], f, r, ]),
                        type = SC_like, stdev = )
            })
          })
        })
      })
    })

  } else {
    loglike_SC_ymafr <- 0
  }

  ## CKMR ----
  if (length(POP)) {
    p_POP <- lapply(1:ns, function(s) {
      calc_POP(t = POP[[s]]$t, a = POP[[s]]$a, y = POP[[s]]$y,
               N = apply(Nsp_yars[, , , s], 1:2, sum),
               fec = fec_yas[, , s])
    })
    loglike_CKMR_POP <- lapply(1:ns, function(s) like_CKMR(n = POP[[s]]$n, m = POP[[s]]$m, p = p_POP[[s]], type = d@CKMR@like))
  } else {
    loglike_CKMR_POP <- 0
  }

  if (length(HSP)) {
    Z_yas <- F_yas + M_yas
    p_HSP <- lapply(1:ns, function(s) {
      calc_HSP(yi = HSP[[s]]$yi, yj = HSP[[s]]$yj,
               N = apply(Nsp_yars[, , , s], 1:2, sum),
               fec = fec_yas[, , s], Z = Z_yas[, , s])
    })
    loglike_CKMR_HSP <- lapply(1:ns, function(s) like_CKMR(n = HSP[[s]]$n, m = HSP[[s]]$m, p = p_HSP[[s]], type = d@CKMR@like))
  } else {
    loglike_CKMR_HSP <- 0
  }

  #loglike_tag <- 0

  loglike <- sum(loglike_I_ymi) + sum(loglike_IAA_ymi) + sum(loglike_IAL_ymi) +
    sum(loglike_CAA_ymfr) + sum(loglike_CAL_ymfr) + sum(loglike_SC_ymafr) +
    Reduce(sum, loglike_CKMR_POP) + Reduce(sum, loglike_CKMR_HSP) #+ sum(loglike_tag)

  # Priors ----
  sdr_s <- exp(p$log_sdr_s)
  rbc_s <- -0.5 * sdr_s * sdr_s

  logprior_rdev_ys <- sapply(1:ns, function(s) { # Check for duplicated rec devs
    ifelse(is.na(map$log_rdev_ys[, s]), 0, dnorm(p$log_rdev_ys[, s], rbc_s[s], sdr_s[s], log = TRUE))
  })

  logprior_dist <- 0
  #logprior_dist <- sapply(1:ns, function(s))

  logprior <- sum(logprior_rdev_ys) + sum(logprior_dist)

  # Objective function ----
  fn <- -1 * (logprior + loglike) + penalty

  # Report out variables ----
  REPORT(CN_ymafrs)
  REPORT(CN_ymlfrs)

  REPORT(loglike)

  REPORT(logprior_rdev_ys)
  REPORT(logprior_dist)
  REPORT(fn)

  return(fn)
}
