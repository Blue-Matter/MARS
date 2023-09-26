
#' Newton-Raphson search for fishing mortality
#'
#' Performs a root finding routine to find the index of F that minimizes the difference between
#' observed catch and the value predicted by the Baranov equation.
#'
#' @param Cobs Observed catch. Vector by `f`
#' @param N Stock abundance at the beginning of the time step. Array `[a, r, s]`
#' @param sel Selectivity. Array `[a, f, s]`
#' @param wt Fishery weight at age. Array `[a, f, s]`
#' @param M Instantaneous natural mortality. Units of per year `[a, s]`
#' @param fleet_area Integers that specify the region `r` in which fleet `f` operates. Vector by `f`
#' @param q_fs Relative catchability of stock `s` for fleet `f`. Matrix `[f, s]`
#' @param delta Numeric, the duration of time in years corresponding to the observed catch, e.g., 0.25 is a quarterly time step.
#' @param na Integer, number of age classes
#' @param nr Integer, number of regions
#' @param nf Integer, number of fleets
#' @param ns Integer, number of stocks
#' @param Fmax Numeric, the maximum Findex value
#' @param nitF Integer, number of iterations for the Newton-Raphson routine
#' @param trans Whether to perform the search in log or logit space
#' @return
#' A list containing:
#'
#' - `F_afs` Fishing mortality array
#' - `F_ars` Fishing mortality array
#' - `F_index` Index of fishing mortality. Vector by `f`
#' - `CB_fs` Catch (biomass) matrix
#' - `CN_afs` Catch (abundance) array
#' - `VB_afs` Vulnerable biomass at the beginning of the time step.
#' - `penalty` Penalty term returned by \link{posfun} when `F_index` exceeds `Fmax`
#' - `fn` Difference between predicted and observed catch at the last iteration. Vector by `f`
#' - `gr` Gradient of `fn` with respect to `F_index` in either log or logit space at the last iteration. Vector by `f`
#'
#' @details
#' The predicted catch for fleet `f` is
#' \deqn{
#' C^{\textrm{pred}}_f = \sum_s \sum_a v_{a,f,s} q_{f,s} F_f \dfrac{1 - \exp(-Z_{a,s})}{Z_{a,s}} N_{a,f,s} w_{a,f,s}
#' }
#'
#' The Newton-Raphson routine minimizes \eqn{f(\vec{x}) = \vec{C}^{\textrm{pred}} - \vec{C}^{\textrm{obs}}} where the vector arrow indexes fleet.
#'
#' If `trans = "log"`, \eqn{\vec{F} = \exp(\vec{x})}.
#'
#' If `trans = "logit"`, \eqn{\vec{F} = F_{\textrm{max}}/(1 + \exp(-\vec{x}))}.
#'
#' The gradient with respect to \eqn{\vec{x}} is
#' \deqn{
#' f'(\vec{x}) = \sum_s \sum_a v_{a,f,s} q_{f,s} N_{a,f,s} w_{a,f,s} \left(\dfrac{\alpha\gamma}{\beta}\right)'
#' }
#'
#' \deqn{
#' \left(\dfrac{\alpha\gamma}{\beta}\right)' = \dfrac{(\alpha\gamma' + \alpha'\gamma)\beta - \alpha\gamma\beta'}{\beta^2}
#' }
#'
#' where
#'
#' \tabular{l}{
#' \eqn{\alpha_f = F_f} \cr
#' \eqn{\beta_{a,s} = Z_{a,s} = M_{a,s} + \sum_f v_{a,f,s} q_{f,s} F_f} \cr
#' \eqn{\gamma_{a,s} = 1 - \exp(-Z_{a,s})} \cr
#' \eqn{\beta'_{a,f,s} = v_{a,f} q_{f,s} \alpha'_f} \cr
#' \eqn{\gamma'_{a,f,s} = \exp(-Z_{a,s})\beta'_{a,f,s}}
#' }
#'
#' If `trans = "log"`, \eqn{\alpha'_f = \alpha_f}.
#'
#' If `trans = "logit"`, \eqn{\alpha'_f = F_{\textrm{max}}\exp(-x_f)/(1 + \exp(-x_f))^2}.
#'
#' This function solves for \eqn{\vec{x}} by iterating until \eqn{f(\vec{x})} approaches zero. In iteration \eqn{i+1}:
#' \deqn{\vec{x}_{i+1} = \vec{x}_i - \dfrac{f(\vec{x}_i)}{f'(\vec{x}_i)}}.
#' @author Q. Huynh
#' @export
Newton_F <- function(Cobs, N, sel, wt, M, fleet_area, q_fs, delta = 1,
                     na = dim(N)[1], nr = dim(N)[2], ns = dim(N)[3], nf = length(Cobs),
                     Fmax = 2, nitF = 5L, trans = c("log", "logit")) {

  trans <- match.arg(trans)
  if (missing(fleet_area) && ns == 1 && nr == 1 && nf == 1) fleet_area <- 1
  if (missing(q_fs) && ns == 1 && nr == 1 && nf == 1) q_fs <- 1
  if (any(fleet_area > nr)) stop("Values in fleet_area vector can't exceed ", nr)

  N <- array(N, c(na, nr, ns))
  sel <- array(sel, c(na, nf, ns))
  wt <- array(wt, c(na, nf, ns))
  M <- matrix(M, na, ns)
  q_fs <- matrix(q_fs, nf, ns)

  # Initialize search ----
  fn <- gr <- x_loop <- matrix(NA_real_, nitF + 1, nf)

  VB_afs <- sapply2(1:ns, function(s) {
    sapply(1:nf, function(f) N[, fleet_area[f], s] * sel[, f, s] * wt[, f, s])
  })
  VB_init <- apply(VB_afs, 2, sum)
  F_init <- Cobs/(Cobs + VB_init)

  if (trans == "log") {
    x_loop[1, ] <- log(F_init)
  } else {
    x_loop[1, ] <- qlogis(F_init/Fmax)
  }

  # Run search for Findex ----
  penalty <- 0 # posfun penalty if F_index > Fmax
  ln_Fmax <- log(Fmax)
  for(i in seq(1, nitF + 1)) {
    if (trans == "log") {
      if (i < nitF + 1) {
        F_loop <- exp(x_loop[i, ])
      } else {  # Last iteration just calculates f and g
        F_loop <- CondExpGt(x_loop[i, ], ln_Fmax, Fmax, exp(x_loop[i, ]))
        penalty <- penalty + sum(posfun(ln_Fmax, x_loop[i, ]))
      }
    } else {
      F_loop <- Fmax * plogis(x_loop[i, ])
    }

    F_fs <- sapply(1:ns, function(s) q_fs[, s] * F_loop)
    F_afs <- sapply2(1:ns, function(s) {
      sapply(1:nf, function(f) F_fs[f, s], sel[, f, s])
    })

    Z_as <- apply(F_afs, c(1, 3), sum) + delta * M
    .gamma_as <- 1 - exp(-Z_as)

    CN_afs <- sapply2(1:ns, function(s) {
      sapply(1:nf, function(f) F_afs[, f, s] * N[, fleet_area[f], s] * .gamma_as[, s] / Z_as[, s])
    })
    CB_f <- apply(CN_afs * wt, 2, sum)

    fn[i, ] <- CB_f - Cobs

    if (trans == "log") {
      deriv_F <- F_loop # f
    } else {
      deriv_F <- Fmax * exp(-x_loop[i, ])/(1 + exp(-x_loop[i, ]))/(1 + exp(-x_loop[i, ]))
    }
    deriv_Z_afs <- sapply2(1:ns, function(s) {
      sapply(1:nf, function(f) q_fs[f, s] * deriv_F[f] * sel[, f, s])
    })
    deriv_gamma_afs <- sapply2(1:ns, function(s) {
      sapply(1:nf, function(f) exp(-Z_as[, s]) * deriv_Z_afs[, f, s])
    })

    constants <- sapply2(1:ns, function(s) { # a s f after aperm
      sapply(1:nf, function(f) q_fs[f, s] * sel[, f, s] * N[, fleet_area[f], s] * wt[, f, s])
    }) %>%
      aperm(c(1, 3, 2))

    deriv1 <- sapply2(1:ns, function(s) outer(deriv_F, .gamma_as[, s]) + F_loop * t(deriv_gamma_afs[, , s])) # f a s
    deriv2 <- sapply2(1:nf, function(f) deriv1[f, , ] * Z_as) # a s f
    deriv3 <- sapply2(1:ns, function(s) deriv2[, s, ] -  outer(.gamma_as[, s], F_loop) * deriv_Z_afs[, , s]) # a f s
    deriv4 <- sapply2(1:nf, function(f) deriv3[, f, ]/Z_as/Z_as) # a s f

    gr[i, ] <- apply(constants * deriv4, 3, sum)

    if (i < nitF + 1) {
      x_loop[i+1, ] <- x_loop[i, ] - fn[i, ]/gr[i, ]
    }
  }

  F_ars <- sapply2(1:nr, function(r) { # a s r before aperm
    i <- fleet_area == r
    apply(F_afs[, i, , drop = FALSE], c(1, 3), sum)
  }) %>%
    aperm(c(1, 3, 2))

  CB_fs <- apply(CN_afs * wt, 2:3, sum)

  list(F_afs = F_afs, F_ars = F_ars, F_index = F_loop,
       CB_fs = CB_fs, CN_afs = CN_afs, VB_afs = VB_afs,
       penalty = penalty, fr = fn[i+1, ], gr = gr[i+1, ])
}



conv_selpar <- function(x, type, nf = length(type), na, Lmax) {

  sel_par <- sapply(1:nf, function(f) {
    sd_asc <- exp(x[2, f])
    sd_desc <- exp(x[3, f])
    if (grepl("age", type[f])) {
      Aapical <- na * plogis(x[1, f])
      v <- c(Aapical, sd_asc, sd_desc)
    } else {
      Lapical <- Lmax * plogis(x[1, f])
      v <- c(Lapical, sd_asc, sd_desc)
    }
    return(v)
  })

  return(sel_par)
}

calc_sel_len <- function(sel_par, lmid, type, nf = length(type)) {

  sel_len <- sapply(1:nf, function(f) {

    if (grepl("length", type[f])) {
      ex_asc <- (lmid - sel_par[1, f])/sel_par[2, f]
      ex2_asc <- -1 * ex_asc^2
      asc <- 2^ex2_asc

      if (grepl("logistic", type[f])) {
        desc <- 1
      } else {
        ex_desc <- (lmid - sel_par[1, f])/sel_par[3, f]
        ex2_desc <- -1 * ex_desc^2
        desc <- 2^ex2_desc
      }
      v <- CondExpLt(lmid, sel_par[1, f], asc, desc)
    } else {
      v <- rep(NA_real_, length(lmid))
    }
    return(v)
  })

  return(sel_len)
}

calc_sel_age <- function(sel_par, sel_len, nf, type) {

}

calc_phi <- function(Z, fec, spawn_time_frac = 0) {
  NPR <- calc_NPR(exp(-Z))
  sum(NPR * fec * exp(-Z * spawn_time_frac))
}

calc_NPR <- function(surv, na = length(surv), plusgroup = TRUE) {
  NPR <- numeric(na)
  NPR[1] <- 1
  for(a in 2:na) NPR[a] <- NPR[a-1] * surv[a-1]
  if (plusgroup) NPR[na] <- NPR[na]/(1 - surv[na])
  return(NPR)
}


predict_recruitment <- function(x, SRR = c("BH", "Ricker"), eq = FALSE, ...) {
  SRR <- match.arg(SRR)
  dots <- list(...)

  if (all(names(dots) %in% c("h", "R0", "phi0"))) {
    h <- dots$h
    R0 <- dots$R0
    phi0 <- dots$phi0

    if (SRR == "BH") {
      a <- 4*h/(1-h)/phi0
      b <- (5*h)/(1-h)/phi0/R0
    } else {
      a <- (5*h)^1.25/phi0
      b <- log((5*h)^1.25)/phi0/R0
    }

  } else if (all(names(dots) %in% c("a", "b"))) {
    a <- dots$a
    b <- dots$b
  } else {
    stop("No stock recruit parameters found")
  }

  if (SRR == "BH") {
    val <- if (eq) a*x/(1 + b*x) else (a*x - 1)/b/x
  } else {
    val <- if (eq) a*x*exp(-b*x) else log(a*x)/b/x
  }
  return(val)
}

SRalphaconv <- function(h, phi, SRR = c("BH", "Ricker")) {
  SRR <- match.arg(SRR)
  switch(
    SRR,
    "BH" = 4*h/(1-h)/phi,
    "Ricker" = (5*h)^1.25/phi
  )
}

SRbetaconv <- function(h, R0, phi, SRR = c("BH", "Ricker")) {
  SRR <- match.arg(SRR)
  switch(
    SRR,
    "BH" = (5*h-1)/(1-h)/phi/R0,
    "Ricker" = log((5*h)^1.25)/phi/R0
  )
}

SRhconv <- function(k, SRR = c("BH", "Ricker")) {
  SRR <- match.arg(SRR)
  switch(
    SRR,
    "BH" = k/(4+k),
    "Ricker" = 0.2*k^0.8
  )
}

calc_LAK <- function(len_a, sd_la, lbin, nl = length(lbin) - 1) {
  stopifnot(length(sd_la) == length(len_a))

  lak_al <- sapply(1:nl, function(j) {
    if (j == nl) {
      1 - pnorm(lbin[j], len_a, sd_la)
    } else if (j == 1) {
      pnorm(lbin[j+1], len_a, sd_la)
    } else {
      pnorm(lbin[j+1], len_a, sd_la) - pnorm(lbin[j], len_a, sd_la)
    }
  })

  return(lak_al)
}


#' Project stock abundance to the next time step
#'
#' This function generates the abundance array by calculating survival from current mortality, then advances
#' age classes, re-distributes the stock, and adds recruitment.
#'
#' @param N Abundance at current time step. Array `[a, r, s]`
#' @param surv Survival during the current time step. Array `[a, r, s]`
#' @param na Integer, number of age classes
#' @param nr Integer, number of regions
#' @param ns Integer, number of stocks
#' @param advance_age Logical, whether the animals advance to their next age class
#' @param R Incoming recruitment. Only assigned if `advance_age = TRUE`. Vector of `s`
#' @param type Character that indicates whether the stock distribution is assigned based on a distribution array or a movement array.
#' @param dist Distribution of the stock in the next time step. Array `[a, r, s]`
#' @param mov Movement array in the next time step. Array `[a, r, r, s]`
#' @return Abundance at the next time step. Array `[a, r, s]`
#' @export
calc_nextN <- function(N, surv, na = dim(N)[1], nr = dim(N)[2], ns = dim(N)[3],
                       advance_age = TRUE, R = numeric(ns),
                       type = c("dist", "mov"), dist = array(1/nr, c(na, nr, ns)), mov) {

  type <- match.arg(type)

  N <- array(N, c(na, nr, ns))
  surv <- array(surv, c(na, nr, ns))
  dist <- array(dist, c(na, nr, ns))
  mov <- array(mov, c(na, nr, nr, ns))

  # Apply survival and advance age class ----
  if (advance_age) {
    Nsurv_ars <- array(NA_real_, c(na, nr, ns))
    Nsurv_ars[2:na, , ] <- N[2:na - 1, , ] * surv[2:na - 1, , ]
    Nsurv_ars[na, , ] <- Nsurv_ars[na, , ] + N[na, , ] * surv[na, , ]
  } else {
    Nsurv_ars <- array(N * surv, c(na, nr, ns))
  }

  # Distribute stock ----
  Nnext_ars <- sapply2(1:ns, function(s) {

    if (type == "dist") { # Distribution vectors
      Ntotal <- apply(Nsurv_ars[, , s], 1, sum)
      if (advance_age) Ntotal[1] <- R_s
      Nout_ar <- Ntotal * dist[, , s]
    } else { # Movement matrix
      Nout_ra <- sapply(2:na, function(a) Nsurv_ars[a, , s] %*% mov[a, , , s])
      if (advance_age) R_r <- R[s] * dist[1, , s]
      Nout_ar <- rbind(R_r, t(Nout_ra))
    }
    return(Nout_ar)
  })

  return(Nnext_ars)
}
