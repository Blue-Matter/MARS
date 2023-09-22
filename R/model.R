

model <- function(data, parameters, map = list(), random = NULL, silent = TRUE, control = list(),
                  run_model = TRUE, do_sd = TRUE, ...) {

  RTMB::TapeConfig(comparison = "tape")
  func <- function(p) .model(p, d = data)

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

.model <- function(p = list(), d = list()) {

  # Assign data variables to environment, see OBS() for simulation ----
  getAll(d)

  # Transform data ----
  delta <- 1/nm

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
  LAK_ymals <- array(NA_real_, c(ny, nm, na, nl, ns))

  # Fishery arrays ----
  FM_ymafs <-
    sel_ymafs <- array(NA_real_, c(ny, nm, na, nf, ns))

  CN_ymafs <- array(NA_real_, c(ny, nm, na, nf, ns))
  CN_ymlfs <- array(NA_real_, c(ny, nm, nl, nf, ns))
  CB_ymfs <- array(NA_real_, c(ny, nm, nf, ns))

  VB_ymfs <- array(NA_real_, c(ny, nm, nf, ns))

  # Transform parameters ----
  R0_s <- exp(p$R0x)/scale_s
  kappa_s <- exp(p$k_s) + 1

  B0_s <- R0_s * phi_s

  Rdev_ys[] <- exp(p$log_Rdev)

  # Fishery selectivity
  fsel_val <- conv_selpar(p$fsel, type = fsel_type, na = na, Lmax = 0.9 * max(len_ymas))
  fsel_len <- calc_sel_len()

  # Miscellaneous penalty term, e.g., F > Fmax
  penalty <- 0

  # Stock recruit parameters ----
  h_s <- sapply(1:ns, function(s) SRhconv(kappa_s[s], SRR = SRR_s[s]))
  sralpha_s <- kappa_s/phi_s
  srbeta_s <- sapply(1:ns, function(s) SRbetaconv(h_s[s], R0_s[s], phi_s[s], SRR = SRR_s[s]))

  # First year, first season initialization ----

  # Loop over years and seasons ----
  for(y in 1:ny) {
    for(m in 1:nm) {
      # Stock distribution parameters
      #dist_ymars[y, m, , , ] <- conv_dist()

      # Calculate length-age key and fishery age selectivity ----
      LAK_ymals[y, m, , , ] <- sapply2(1:ns, function(s) calc_LAK(len_ymas[y, m, , s], sdlen_ymas[y, m, , s], lbin))
      sel_ymafs[y, m, , , ] <- calc_sel_age()

      ## This season's mortality ----
      Fsearch <- Newton_F(
        Cobs = Cobs_ymf[y, m, ], N = N_ymars[y, m, , , ], sel = sel_ymafs[y, m, , , ],
        wt = wt_yafs[y, , , ], M = M_yas[y, , s],
        fleet_area = fleet_area, q_fs = q_fs, delta = delta,
        na = na, nr = nr, nf = nf, ns = ns, Fmax = Fmax, nitF = nitF, trans = trans
      )
      penalty <- penalty + Fsearch[["penalty"]] # add penalty for exceeding Fmax

      FM_ymars[y, m, , , ] <- Fsearch[["F_ars"]]
      Z_ymars[y, m, , , ] <- sapply2(1:nr, function(r) { # a s r before aperm
        FM_ymars[y, m, , r, ] + delta * M_yas[y, , ]
      }) %>%
        aperm(c(1, 3, 2))

      ## This season's fishery catch ----
      CN_ymafs[y, m, , , ] <- Fsearch[["CN_afs"]]
      CN_ymlfs[y, m, , , ] <- sapply2(1:ns, function(s) {
        sapply(1:nf, function(f) CN_ymafs[y, m, , f, s] %*% LAK_ymals[y, m, , , s])
      })
      CB_ymfs[y, m, , ] <- Fsearch[["CB_fs"]]
      VB_ymfs[y, m, , ] <- apply(Fsearch[["VB_afs"]], 2:3, sum)

      ## This year's spawning and recruitment ----
      if (m == m_spawn) {
        Nsp_yars[y, , , ] <- N_ymars[y, m, , , ] * exp(-spawn_time_frac * Z_ymars[y, m, , , ])
        SB_ys[y, ] <- sapply(1:ns, function(s) {
          sapply(1:nr, function(r) Nsp_yars[y, , r, s] * fec_yas[y, , s]) %>% sum()
        })
        R_ys[y, ] <- Rdev_ys[y, ] * sapply(1:ns, function(s) {
          predict_recruitment(SB = SB[y, s], SRR = SRR_s[s], a = sralpha_s[s], b = srbeta_s[s])
        })
      }

      ## Next season's abundance ----
      # Need to advance ages, assign recruitment, and merge plus group
      Nnext_as <- array(N_ymars[y, m, , , ] * exp(-Z_ymars[y, m, , , ]), c(na, nr, ns)) %>% apply(c(1, 3), sum)
      if (m == nm) {
        N_ymars[y+1, 1, , , ] <- sapply2(1:ns, function(s) Nnext_as[, s] * dist_ymars[y+1, 1, , , s])
      } else {
        N_ymars[y, m+1, , , ] <- sapply2(1:ns, function(s) Nnext_as[, s] * dist_ymars[y, m+1, , , s])
      }
    }
  }

  # Index arrays ----

  # Index selectivity ----

  # Index ----

  # Likelihoods ----
  #log_like_index <- 0
  #log_like_IAA <- 0
  #log_like_IAL <- 0
  #log_like_CAA <- 0
  #log_like_CAL <- 0
  #log_like_SOO <- 0
  #log_like_tag <- 0
  #log_like <- sum(log_like_index) + sum(log_like_IAA) + sum(log_like_IAL) +
  #  sum(log_like_CAA) + sum(log_like_CAL) + sum(log_like_SOO) + sum(log_like_tag)

  # Priors ----
  log_prior_rdev <- sapply(1:ns, function(s) {
    ifelse(est_rdev[, s], dnorm(p$log_rdev[, s], -0.5*sdr*sdr, sdr, log = TRUE), 0)
  })
  #log_prior_dist <- sapply(1:ns, function(s))
  log_prior <- sum(log_prior_rdev) + sum(log_prior_dist)

  # Objective function ----
  fn <- -1 * (log_prior + log_like) + penalty

  # Report out variables ----


  return(fn)
}
