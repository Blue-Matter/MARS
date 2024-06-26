
#' Multi-fleet, multi-area, multi-stock population dynamics model
#'
#' Project age-structured populations forward in time. Frequently used to calculate
#' equilibrium abundance and biomass for which there is no analytic solution
#' due to seasonal movement.
#'
#' @param ny Integer, number of years for the projection
#' @param nm Integer, number of seasons
#' @param initN_ars Abundance in the first year, first season. Array `[a, r, s]`
#' @param mov_ymarrs Movement array `[y, m, a, r, r, s]`. If missing, uses a diagonal matrix (no movement among areas).
#' @param M_yas Natural mortality (per year). Array `[y, a, s]`
#' @param SRR_s Character vector by `s` for the stock recruit relationship. See [calc_recruitment()] for options
#' @param sralpha_s Numeric vector by `s` for the stock recruit alpha parameter
#' @param srbeta_s Numeric vector by `s` for the stock recruit beta parameter
#' @param mat_yas Maturity ogive. Array `[y, a, s]`
#' @param fec_yas Fecundity schedule (spawning output of mature individuals). Array `[y, a, s]`
#' @param Rdev_ys Recruitment deviations. Matrix `[y, s]`
#' @param m_spawn Integer, season of spawning
#' @param m_rec Integer, season of recruitment
#' @param delta_s Numeric vector by `s`. Fraction of season that elapses when spawning occurs, e.g., midseason spawning when `delta_s = 0.5`.
#' @param natal_rs Matrix `[r, s]`. The fraction of the mature stock `s` in region `r` that spawns at
#' time of spawning. See example in [Dstock-class].
#' @param recdist_rs Matrix `[r, s]`. The fraction of the incoming recruitment of stock `s` that settles in region `r`.
#' @param fwt_ymafs Fishery weight at age. Array `[y, m, a, f, s]`
#' @param sel_ymafs Fishery selectivity. Array `[y, m, a, f, s]`
#' @param condition Whether the fishing mortality is conditioned on the catch or specified F argument.
#' @param F_ymfr Fishing mortality (per season). Array `[y, m, f, r]`. Only used if `condition = "F"`.
#' @param Cobs_ymfr Fishery catch (weight). Array `[y, m, f, r]`. Only used if `condition = "catch"` to solve for F (see [calc_F()]).
#'
#' @examples
#' unfished_pop <- calc_population()
#' @inheritParams calc_F
#' @return
#' A named list containing:
#'
#' * `N_ymars` Stock abundance
#' * `F_ymars` Fishing mortality (summed across fleets)
#' * `F_ymfr` Fishing mortality (by fleet and region)
#' * `Z_ymars` Total mortality
#' * `F_ymafrs` Fishing mortality (disaggregated by fleet)
#' * `CN_ymafrs` Catch at age (abundance)
#' * `CB_ymfrs` Fishery catch (weight)
#' * `VB_ymfrs` Vulnerable biomass available to the fishing fleets
#' * `Nsp_yars` Spawning abundance (in the spawning season)
#' * `Npsp_yars` Mature abundance (that does not if outside natal regions)
#' * `S_yrs` Spawning output
#' * `R_ys` Recruitment
#' * `penalty` Numeric quadratic penalty if apical fishing mortality (by fleet) exceeds `Fmax`. See [calc_F()].
#'
#' @export
calc_population <- function(ny = 10, nm = 4, na = 20, nf = 1, nr = 4, ns = 2,
                            initN_ars = array(1, c(na, nr, ns)),
                            mov_ymarrs,
                            M_yas = array(0.3, c(ny, na, ns)),
                            SRR_s = rep("BH", ns),
                            sralpha_s = rep(1e16, ns),
                            srbeta_s = rep(1e16, ns),
                            mat_yas = array(1, c(ny, na, ns)), fec_yas = array(1, c(ny, na, ns)),
                            Rdev_ys = matrix(1, ny, ns),
                            m_spawn = 1, m_rec = 1, delta_s = rep(0, ns), natal_rs = matrix(1, nr, ns),
                            recdist_rs = matrix(1/nr, nr, ns),
                            fwt_ymafs = array(1, c(ny, nm, na, nf, ns)),
                            q_fs = matrix(1, nf, ns), sel_ymafs = array(1, c(ny, nm, na, nf, ns)),
                            condition = c("F", "catch"),
                            F_ymfr = array(0, c(ny, nm, nf, nr)),
                            Cobs_ymfr = matrix(1e-8, c(ny, nm, nf, nr)),
                            Fmax = 2, nitF = 5L) {

  # Dispatch method for AD variables ----
  `[<-` <- RTMB::ADoverload("[<-")

  condition <- match.arg(condition)

  # Population arrays ----
  delta_m <- 1/nm
  N_ymars <- array(0, c(ny + 1, nm, na, nr, ns))
  Npsp_yars <-
    Nsp_yars <- array(NA_real_, c(ny, na, nr, ns))
  S_yrs <- array(NA_real_, c(ny, nr, ns))
  R_ys <- array(NA_real_, c(ny, ns))
  F_ymars <-
    Z_ymars <- array(NA_real_, c(ny, nm, na, nr, ns))
  if (missing(mov_ymarrs)) {
    mov_ymarrs <- diag(nr) %>% array(c(nr, nr, ny, nm, na, ns)) %>% aperm(c(3:5, 1:2, 6))
  }

  # Fishery arrays ----
  if (condition == "catch") {
    F_ymfr <- array(NA_real_, c(ny, nm, nf, nr))
  }
  F_ymafrs <-
    CN_ymafrs <- array(NA_real_, c(ny, nm, na, nf, nr, ns))
  CB_ymfrs <-
    VB_ymfrs <- array(NA_real_, c(ny, nm, nf, nr, ns))

  # Penalty for exceeding Fmax (if conditioning on catch)
  penalty <- 0

  N_ymars[1, 1, , , ] <- initN_ars

  # Loops over years and seasons ----
  for(y in 1:ny) {
    for(m in 1:nm) {
      ## This season's mortality ----
      if (condition == "catch") {
        Fsearch <- calc_F(
          Cobs = Cobs_ymfr[y, m, , ], N = N_ymars[y, m, , , ], sel = sel_ymafs[y, m, , , ],
          wt = fwt_ymafs[y, m, , , ], M = M_yas[y, , ], q_fs = q_fs, delta = delta_m,
          na = na, nr = nr, nf = nf, ns = ns, Fmax = Fmax, nitF = nitF, trans = "log"
        )
        penalty <- penalty + Fsearch[["penalty"]] # Report penalty for exceeding Fmax

        F_ymafrs[y, m, , , , ] <- Fsearch[["F_afrs"]]
        F_ymars[y, m, , , ] <- Fsearch[["F_ars"]]
        Z_ymars[y, m, , , ] <- Fsearch[["Z_ars"]]
        F_ymfr[y, m, , ] <- Fsearch[["F_index"]]

        ## This season's fishery catch, vulnerable biomass, and total biomass ----
        CN_ymafrs[y, m, , , , ] <- Fsearch[["CN_afrs"]]
        CB_ymfrs[y, m, , , ] <- Fsearch[["CB_frs"]]
        VB_ymfrs[y, m, , , ] <- apply(Fsearch[["VB_afrs"]], 2:4, sum)
      } else {

        ind_afrs <- as.matrix(expand.grid(y = y, m = m, a = 1:na, f = 1:nf, r = 1:nr, s = 1:ns))
        fs_afrs <- ind_afrs[, c("f", "s")]
        ymfr_afrs <- ind_afrs[, c("y", "m", "f", "r")]
        ymafs_afrs <- ind_afrs[, c("y", "m", "a", "f", "s")]
        ymars_afrs <- ind_afrs[, c("y", "m", "a", "r", "s")]

        F_ymafrs[y, m, , , , ] <- q_fs[fs_afrs] * F_ymfr[ymfr_afrs] * sel_ymafs[ymafs_afrs]
        F_ymars[y, m, , , ] <- apply(F_ymafrs[y, m, , , , , drop = FALSE], c(3, 5:6), sum)

        ind_ars <- as.matrix(expand.grid(y = y, m = m, a = 1:na, r = 1:nr, s = 1:ns))
        yas_ars <- ind_ars[, c("y", "a", "s")]

        Z_ymars[y, m, , , ] <- F_ymars[y, m, , , ] + delta_m * M_yas[yas_ars]
        CN_ymafrs[y, m, , , , ] <- F_ymafrs[y, m, , , , ] * (1 - exp(-Z_ymars[ymars_afrs])) * N_ymars[ymars_afrs] / Z_ymars[ymars_afrs]

        ymafs_afrs <- ind_afrs[, c("y", "m", "a", "f", "s")]

        CB_ymfrs[y, m, , , ] <- array(CN_ymafrs[ind_afrs] * fwt_ymafs[ymafs_afrs], c(na, nf, nr, ns)) %>%
          apply(2:4, sum)

        VB_ymfrs[y, m, , , ] <- array(
          sel_ymafs[ymafs_afrs] * fwt_ymafs[ymafs_afrs] * N_ymars[ymars_afrs],
          c(na, nf, nr, ns)
        ) %>%
          apply(2:4, sum)
      }

      ## This year's spawning and recruitment ----
      if (m == m_spawn) {
        Npsp_yars[y, , , ] <- sapply2(1:ns, function(s) {
          sapply(1:nr, function(r) {
            N_ymars[y, m, , r, s] * exp(-delta_s[s] * Z_ymars[y, m, , r, s]) * mat_yas[y, , s]
          })
        })
        Nsp_yars[y, , , ] <- sapply2(1:ns, function(s) {
          sapply(1:nr, function(r) natal_rs[r, s] * Npsp_yars[y, , r, s])
        })
        S_yrs[y, , ] <- sapply(1:ns, function(s) {
          sapply(1:nr, function(r) sum(Nsp_yars[y, , r, s] * fec_yas[y, , s]))
        })
        R_ys[y, ] <- Rdev_ys[y, ] * sapply(1:ns, function(s) {
          calc_recruitment(sum(S_yrs[y, , s]), SRR = SRR_s[s], a = sralpha_s[s], b = srbeta_s[s])
        })
      }

      ## Next season's abundance and total biomass ----
      # Have to advance the age classes in the season before spawning
      ynext <- ifelse(m == nm, y+1, y)
      mnext <- ifelse(m == nm, 1, m+1)
      ylast <- min(ynext, ny)

      y_spawn <- ifelse(m_rec > m_spawn, ynext, ynext-1) # If m_rec <= m_spawn in year 1, multiseason?

      Rnext <- if (y_spawn > 0) { # For year 1, no recruitment has been calculated yet
        if (y_spawn <= ny) {
          R_ys[y_spawn, ]
        } else {
          0
        }
      } else {
        apply(initN_ars[1, , , drop = FALSE], 3, sum)
      }

      N_ymars[ynext, mnext, , , ] <- calc_nextN(
        N = N_ymars[y, m, , , ], surv = exp(-Z_ymars[y, m, , , ]),
        na = na, nr = nr, ns = ns,
        advance_age = mnext == m_rec,
        R = Rnext,
        recdist = recdist_rs,
        mov = mov_ymarrs[ylast, mnext, , , , ]
      )
    }
  }

  if (condition == "F") penalty <- penalty + sum(posfun(Fmax, F_ymfr))

  out <- list(
    N_ymars = N_ymars,
    F_ymars = F_ymars,
    F_ymfr = F_ymfr,
    Z_ymars = Z_ymars,
    F_ymafrs = F_ymafrs,
    CN_ymafrs = CN_ymafrs,
    CB_ymfrs = CB_ymfrs,
    VB_ymfrs = VB_ymfrs,
    Nsp_yars = Nsp_yars,
    Npsp_yars = Npsp_yars,
    S_yrs = S_yrs,
    R_ys = R_ys,
    penalty = penalty
  )

  return(out)
}

#' Equilibrium spawners per recruit by projection
#'
#' Project a population forward in time using [calc_population()] with constant recruitment and
#' seasonal dynamics (growth, movement-by-season) to obtain per recruit parameters. Note that the fishing
#' mortality among fleets and stocks remain linked by matrix `q_fs`.
#'
#' @inheritParams calc_population
#'
#' @param F_mfr Equilibrium fishing mortality (per season). Matrix `[m, f, r]`
#' @param sel_mafs Selectivity by season, age, fleet, stock. Array `[m, a, f, s]`
#' @param fwt_mafs Fishery weight array by season, age, fleet, stock. Array `[m, a, r, r]`. Can be used
#' calculate yield per recruit.
#' @param mov_marrs Movement array `[m, a, r, r, s]`. If missing, uses a diagonal matrix (no movement among areas).
#' @param M_as Natural mortality. Matrix `[a, s]`
#' @param mat_as Maturity at age. Matrix `[a, s]`
#' @param fec_as Fecundity at age. Matrix `[a, s]`
#' @details
#' The initial population vector will be the survival at age evenly by the number of regions `nr`.
#' @return A named list returned by [calc_population()].
#' @seealso [calc_phi_simple()]
#' @export
calc_phi_project <- function(ny, nm, na, nf = 1, nr, ns = 1,
                             F_mfr = array(0, c(nm, nf, nr)),
                             sel_mafs = array(1, c(nm, na, nf, ns)),
                             fwt_mafs = array(1, c(nm, na, nf, ns)),
                             q_fs = matrix(1, nf, ns),
                             M_as, mov_marrs,
                             mat_as, fec_as,
                             m_spawn = 1, m_rec = 1,
                             delta_s = rep(0, ns),
                             natal_rs = matrix(1, nr, ns),
                             recdist_rs = matrix(1/nr, nr, ns)) {
  delta_m <- 1/nm

  if (missing(mov_marrs)) {
    mov_ymarrs <- diag(nr) %>% array(c(nr, nr, ny, nm, na, ns)) %>% aperm(c(3:5, 1:2, 6))
  } else {
    mov_marrs <- array(mov_marrs, c(nm, na, nr, nr, ns))
    mov_ymarrs <- array(mov_marrs, c(nm, na, nr, nr, ns, ny)) %>% aperm(c(6, 1:5))
  }
  M_yas <- array(M_as, c(na, ns, ny)) %>% aperm(c(3, 1, 2))

  SRR <- rep("BH", ns)
  sralpha <- srbeta <- rep(1e16, ns)

  mat_yas <- array(mat_as, c(na, ns, ny)) %>% aperm(c(3, 1, 2))
  fec_yas <- array(fec_as, c(na, ns, ny)) %>% aperm(c(3, 1, 2))
  sel_ymafs <- array(sel_mafs, c(nm, na, nf, ns, ny)) %>% aperm(c(5, 1:4))
  fwt_ymafs <- array(fwt_mafs, c(nm, na, nf, ns, ny)) %>% aperm(c(5, 1:4))
  F_ymfr <- array(F_mfr, c(nm, nf, nr, ny)) %>% aperm(c(4, 1:3))

  initNPR_ars <- sapply2(1:ns, function(s) {
    NPR_ar <- sapply(1:nr, function(r) {
      F_af <- lapply(1:nf, function(f) sel_ymafs[1, 1, , f, s] * q_fs[f, s] * F_ymfr[1, 1, f, r])
      Z_a <- M_yas[1, , s] + Reduce("+", F_af)
      calc_NPR(Z_a)
    })
    return(NPR_ar/nr)
  })

  pop_phi <- calc_population(
    ny, nm, na, nf, nr, ns, initN_ars = initNPR_ars,
    mov_ymarrs, M_yas,
    SRR_s = SRR, sralpha_s = sralpha, srbeta_s = srbeta,
    mat_yas, fec_yas, Rdev_ys = matrix(1, ny, ns),
    m_spawn, m_rec, delta_s, natal_rs,
    fwt_ymafs = fwt_ymafs, q_fs,
    sel_ymafs = sel_ymafs,
    condition = "F",
    F_ymfr = F_ymfr
  )
  return(pop_phi)
}

#' Simple spawners per recruit calculation
#'
#' Calculate spawners per recruit by individual stock. Appropriate for a population model with a single
#' spatial area and an annual time steps, i.e. single season.
#'
#' @inheritParams calc_phi_project
#' @param Z Total mortality at age
#' @param mat_a Maturity at age. Vector
#' @param fec_a Fecundity at age. Vector
#' @param delta Fraction of season that elapses when spawning occurs, e.g., midseason spawning when `delta = 0.5`.
#' @return [calc_phi_simple()] returns the spawners per recruit.
#' @seealso [calc_phi_project()]
#' @export
calc_phi_simple <- function(Z, fec_a, mat_a, delta = 0) {
  NPR <- calc_NPR(Z)
  sum(NPR * exp(-Z * delta) * fec_a * mat_a)
}

#' @rdname calc_phi_simple
#' @return [calc_NPR()] returns a vector of the numbers per recruit at age, i.e., survival.
#' @param plusgroup Logical, whether the largest age class is a plusgroup accumulator age
#' @export
calc_NPR <- function(Z, na = length(Z), plusgroup = TRUE) {
  surv <- exp(-Z)
  is_ad <- inherits(Z, "advector")
  if (is_ad) {
    `[<-` <- RTMB::ADoverload("[<-")
    NPR <- advector(numeric(na))
  } else {
    NPR <- numeric(na)
  }
  NPR[1] <- 1
  for(a in 2:na) NPR[a] <- NPR[a-1] * surv[a-1]
  if (plusgroup) NPR[na] <- NPR[na]/(1 - surv[na])
  return(NPR)
}
