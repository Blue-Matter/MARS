
#' Newton-Raphson search for fishing mortality
#'
#' Performs a root finding routine to find the index of F that minimizes the difference between
#' observed catch and the value predicted by the Baranov equation.
#'
#' @param Cobs Observed catch. Matrix `[f, r]`
#' @param N Stock abundance at the beginning of the time step. Array `[a, r, s]`
#' @param sel Selectivity. Array `[a, f, s]`
#' @param wt Fishery weight at age. Array `[a, f, s]`
#' @param M Instantaneous natural mortality. Units of per year `[a, s]`
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
#' - `F_afrs` Fishing mortality array
#' - `F_ars` Fishing mortality array (summed across fleets)
#' - `Z_ars` Total mortality array
#' - `F_index` Index of fishing mortality. Matrix `[f, r]`
#' - `CB_frs` Catch (biomass) array
#' - `CN_afrs` Catch (abundance) array
#' - `VB_afrs` Vulnerable biomass at the beginning of the time step. Array
#' - `penalty` Penalty term returned by [posfun()] when `F_index` exceeds `Fmax`
#' - `fn` Difference between predicted and observed catch at the last iteration. Matrix `[f, r]`
#' - `gr` Gradient of `fn` with respect to `F_index` in either log or logit space at the last iteration. Vector by `[f, r]`
#'
#' @details
#' The predicted catch for fleet `f` in region `r` is
#' \deqn{
#' C^{\textrm{pred}}_{f,r} = \sum_s \sum_a v_{a,f,s} q_{f,s} F_{f,r} \dfrac{1 - \exp(-Z_{a,r,s})}{Z_{a,r,s}} N_{a,r,s} w_{a,f,s}
#' }
#'
#' The Newton-Raphson routine minimizes \eqn{f(x_{f,r}) = C_{f,r}^{\textrm{pred}} - C_{f,r}^{\textrm{obs}}}.
#'
#' If `trans = "log"`, \eqn{F_{f,r} = \exp(x_{f,r})}.
#'
#' If `trans = "logit"`, \eqn{F_{f,r} = F_{\textrm{max}}/(1 + \exp(x_{f,r}))}.
#'
#' The gradient with respect to \eqn{\vec{x}} is
#' \deqn{
#' f'(x_{f,r}) = \sum_s \sum_a v_{a,f,s} q_{f,s} N_{a,r,s} w_{a,f,s} \left(\dfrac{\alpha\gamma}{\beta}\right)'
#' }
#'
#' \deqn{
#' \left(\dfrac{\alpha\gamma}{\beta}\right)' = \dfrac{(\alpha\gamma' + \alpha'\gamma)\beta - \alpha\gamma\beta'}{\beta^2}
#' }
#'
#' where
#'
#' \tabular{l}{
#' \eqn{\alpha_{f,r} = F_{f,r}} \cr
#' \eqn{\beta_{a,r,s} = Z_{a,r,s} = M_{a,s} + \sum_f v_{a,f,s} q_{f,s} F_{f,r}} \cr
#' \eqn{\gamma_{a,r,s} = 1 - \exp(-Z_{a,r,s})} \cr
#' \eqn{\beta'_{a,f,r,s} = v_{a,f} q_{f,s} \alpha'_{f,r}} \cr
#' \eqn{\gamma'_{a,f,r,s} = \exp(-Z_{a,r,s})\beta'_{a,f,r,s}}
#' }
#'
#' If `trans = "log"`, \eqn{\alpha'_{f,r} = \alpha_{f,r}}.
#'
#' If `trans = "logit"`, \eqn{\alpha'_{f,r} = F_{\textrm{max}}\exp(-x_{f,r})/(1 + \exp(-x_{f,r}))^2}.
#'
#' This function solves for \eqn{\vec{x}} by iterating until \eqn{f(\vec{x})} approaches zero, where the vector arrow
#' indexes over fleet and region. In iteration \eqn{i+1}:
#' \deqn{\vec{x}_{i+1} = \vec{x}_i - \dfrac{f(\vec{x}_i)}{f'(\vec{x}_i)}}.
#' @author Q. Huynh
#' @export
calc_F <- function(Cobs, N, sel, wt, M, q_fs, delta = 1,
                   na = dim(N)[1], nr = dim(N)[2], ns = dim(N)[3], nf = length(Cobs),
                   Fmax = 2, nitF = 5L, trans = c("log", "logit")) {

  trans <- match.arg(trans)
  if (missing(q_fs) && ns == 1 && nr == 1 && nf == 1) q_fs <- 1
  if (is.null(dim(Cobs)) && nr == 1 && nf == 1) Cobs <- matrix(Cobs, nf, nr)

  N <- array(N, c(na, nr, ns))
  sel <- array(sel, c(na, nf, ns))
  wt <- array(wt, c(na, nf, ns))
  M <- matrix(M, na, ns)
  q_fs <- matrix(q_fs, nf, ns)

  # Initialize search ----
  fn <- gr <- x_loop <- list()
  #matrix(NA_real_, nitF + 1, nf)

  VB_afrs <- sapply2(1:ns, function(s) {
    sapply2(1:nf, function(r) {
      sapply(1:nf, function(f) N[, r, s] * sel[, f, s] * wt[, f, s])
    })
  })
  VB_fr <- apply(VB_afrs, c(2, 3), sum)
  F_init <- Cobs/(Cobs + VB_fr)

  if (trans == "log") {
    x_loop[[1]] <- log(F_init)
  } else {
    x_loop[[1]] <- qlogis(F_init/Fmax)
  }

  # Run search for Findex ----
  penalty <- 0 # posfun penalty if F_index > Fmax
  ln_Fmax <- log(Fmax)
  for(i in seq(1, nitF + 1)) {
    if (trans == "log") {
      if (i < nitF + 1) {
        F_loop <- exp(x_loop[[1]])
      } else {  # Last iteration just calculates f and g
        F_loop <- CondExpGt(x_loop[[i]], ln_Fmax, Fmax, exp(x_loop[[i]]))
        F_loop <- CondExpLe(Cobs, 1e-8, 0, F_loop)
        penalty <- penalty + sum(posfun(ln_Fmax, x_loop[[i]]))
      }
    } else {
      F_loop <- Fmax * plogis(x_loop[[i]])
    }

    F_frs <- sapply2(1:ns, function(s) q_fs[, s] * F_loop)
    F_afrs <- sapply2(1:ns, function(s) {
      sapply2(1:nr, function(r) {
        sapply(1:nf, function(f) F_frs[f, r, s] * sel[, f, s])
      })
    })
    F_ars <- apply(F_afrs, c(1, 3, 4), sum)
    Z_ars <- sapply2(1:nr, function(r) F_ars[, r, ] + delta * M) %>% # a s r
      aperm(c(1, 3, 2))
    .gamma_ars <- 1 - exp(-Z_ars)

    CN_afrs <- sapply2(1:ns, function(s) {
      sapply2(1:nr, function(r) {
        sapply(1:nf, function(f) calc_Baranov(F_afrs[, f, r, s], Z_ars[, r, s], N[, r, s]))
      })
    })
    CB_afsr <- sapply2(1:nr, function(r) CN_afrs[, , r, ] * wt)
    CB_fr <- apply(CB_afsr, c(2, 4), sum)

    fn[[i]] <- CB_fr - Cobs

    if (trans == "log") {
      deriv_F <- F_loop # f r
    } else {
      deriv_F <- Fmax * exp(-x_loop[[i]])/(1 + exp(-x_loop[[i]]))/(1 + exp(-x_loop[[i]]))
    }
    deriv_Z_fars <- sapply2(1:ns, function(s) {
      sapply2(1:nr, function(r) q_fs[, s] * deriv_F[, r] * t(sel[, , s]))
    })
    deriv_gamma_afrs <- sapply2(1:ns, function(s) {
      sapply2(1:nr, function(r) exp(-Z_ars[, r, s]) * t(deriv_Z_fars[, , r, s]))
    })

    constants_afrs <- sapply2(1:ns, function(s) {
      sapply2(1:nr, function(r) {
        sapply(1:nf, function(f) q_fs[f, s] * sel[, f, s] * N[, r, s] * wt[, f, s])
      })
    })

    deriv1_fars <- sapply(1:ns, function(s) {
      sapply2(1:nr, function(r) {
        outer(deriv_F[, r], .gamma_ars[, r, s]) + F_loop[, r] * t(deriv_gamma_afrs[, , r, s])
      })
    })
    deriv2_arsf <- sapply2(1:nf, function(f) deriv1_fars[f, , , ] * Z_ars)
    deriv3_afrs <- sapply2(1:ns, function(s) {
      sapply2(1:nr, function(r) {
        deriv2_arsf[, r, s, ] -  outer(.gamma_ars[, r, s], F_loop[, r]) * t(deriv_Z_fars[, , r, s])
      })
    })
    deriv4_arsf <- sapply2(1:nf, function(f) deriv3_afrs[, f, , ]/Z_ars/Z_ars)

    gr[[i]] <- apply(constants_afrs * aperm(deriv4_arsf, c(1, 4, 2, 3)), c(2, 33), sum)

    if (i <= nitF) x_loop[[i+1]] <- x_loop[[i]] - fn[[i]]/gr[[i]]
  }

  CB_frs <- apply(CB_afsr, c(2, 4, 3), sum)

  list(F_afrs = F_afrs, F_ars = F_ars, F_index = F_loop, Z_ars = Z_ars,
       CB_frs = CB_frs, CN_afrs = CN_afrs, VB_afrs = VB_afrs,
       penalty = penalty, fr = fn[[i+1]], gr = gr[[i+1]])
}

calc_Baranov <- function(FM, Z, N) FM/Z * (1 - exp(-Z)) * N


#' @importFrom stats uniroot
.calc_summary_F <- function(FM, M, N, CN) calc_Baranov(FM, FM + M, N) - CN
calc_summary_F <- function(M, N, CN, Fmax) {
  out <- uniroot(.calc_summary_F, interval = c(0, Fmax), M = M, N = N, CN = CN)
  out$root
}

calc_phi <- function(Z, fec, delta = 0) {
  NPR <- calc_NPR(exp(-Z))
  sum(NPR * fec * exp(-Z * delta))
}

calc_NPR <- function(surv, na = length(surv), plusgroup = TRUE) {
  NPR <- numeric(na)
  NPR[1] <- 1
  for(a in 2:na) NPR[a] <- NPR[a-1] * surv[a-1]
  if (plusgroup) NPR[na] <- NPR[na]/(1 - surv[na])
  return(NPR)
}


#' Calculate recruitment from stock-recruit function
#'
#' @param x Numeric, either the spawning output or the equilibrium spawners per recruit, from which
#' the recruitment will be calculated. See argument `eq`.
#' @param SRR Character to indicate the functional form of the stock recruit function
#' @param eq Logical, indicates whether `x` is the spawning output (`FALSE`) or equilibrium spawners per recruit (`TRUE`)
#' @param ... Parameters of the SRR function. Provide one of two sets of variables:
#' \itemize{
#' \item `h`, `R0` and `phi0`
#' \item `a` and `b` (alpha, beta values)
#' }
#' @examples
#' calc_recruitment(10, SRR = "Ricker", a = 2, b = 0.5)
#' calc_recruitment(10, SRR = "Ricker", h = 0.9, R0 = 1, phi0 = 1)
#' @export
calc_recruitment <- function(x, SRR = c("BH", "Ricker"), eq = FALSE, ...) {
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
  k <- SRkconv(h, SRR)
  k/phi
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

SRkconv <- function(h, SRR = c("BH", "Ricker")) {
  SRR <- match.arg(SRR)
  switch(
    SRR,
    "BH" = 4*h/(1-h),
    "Ricker" = (5*h)^1.25
  )
}

#' Length-at-age key
#'
#' Calculates the probability distribution of length-at-age using the normal probability density function
#'
#' @param len_a Vector of length-at-age
#' @param sd_la Vector of standard deviation in length-at-age
#' @param lbin Vector of the lower boundary of the length bins
#' @param nl Integer, number of length bins (default is `length(lbin) - 1`)
#' @return Matrix by age (rows) and length (columns)
#' @export
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
# #' @param type Character that indicates whether the stock distribution is assigned based on a distribution array or a movement array.
# #' @param dist Distribution of the stock in the next time step. Array `[a, r, s]`
#' @param mov Movement array in the next time step. Array `[a, r, r, s]`. Rows denote region of origin and columns denote region of destination.
#' @param plusgroup Logical, whether the last age class is an accumulator plus group.
#' @return Abundance at the next time step. Array `[a, r, s]`
#' @export
calc_nextN <- function(N, surv, na = dim(N)[1], nr = dim(N)[2], ns = dim(N)[3],
                       advance_age = TRUE, R = numeric(ns),
                       #type = c("dist", "mov"),
                       #dist = array(1/nr, c(na, nr, ns)),
                       mov = array(1/nr, c(na, nr, nr, ns)),
                       plusgroup = TRUE) {

  type <- match.arg(type)

  N <- array(N, c(na, nr, ns))
  surv <- array(surv, c(na, nr, ns))
  #dist <- array(dist, c(na, nr, ns))
  mov <- array(mov, c(na, nr, nr, ns))

  # Apply survival and advance age class ----
  if (advance_age) {
    Nsurv_ars <- array(NA_real_, c(na, nr, ns))
    Nsurv_ars[2:na, , ] <- N[2:na - 1, , ] * surv[2:na - 1, , ]
    if (plusgroup) Nsurv_ars[na, , ] <- Nsurv_ars[na, , ] + N[na, , ] * surv[na, , ]
  } else {
    Nsurv_ars <- N * surv
  }

  # Distribute stock ----
  Nnext_ars <- sapply2(1:ns, function(s) {
    Nout_ra <- sapply(2:na, function(a) Nsurv_ars[a, , s] %*% mov[a, , , s])
    if (advance_age) {
      R_ar <- R[s]/rep(nr, nr)
      R_r <- R_ar %*% mov[1, , , s]
    }
    Nout_ar <- rbind(R_r, t(Nout_ra))

    #if (type == "dist") { # Distribution vectors
    #  Ntotal <- apply(Nsurv_ars[, , s], 1, sum)
    #  if (advance_age) Ntotal[1] <- R
    #  Nout_ar <- Ntotal * dist[, , s]
    #} else { # Movement matrix
    #  Nout_ra <- sapply(2:na, function(a) Nsurv_ars[a, , s] %*% mov[a, , , s])
    #  if (advance_age) R_r <- R[s] * dist[1, , s]
    #  Nout_ar <- rbind(R_r, t(Nout_ra))
    #}
    return(Nout_ar)
  })

  return(Nnext_ars)
}


#' Calculate index at age
#'
#' For indices of abundance, the function calculates the numbers vulnerable to the survey.
#'
#' @param N Stock abundance at the beginning of the time step. Array `[a, r, s]`
#' @param Z Instantaneous total mortality. Array `[a, r, s]`
#' @param sel Index selectivity. Array `[a, i, s]`
#' @param na Integer, number of age classes
#' @param ni Integer, number of indices
#' @param ns Integer, number of stocks
#' @param samp Boolean indicates which regions and stocks are sampled by the index. Array `[i, r, s]`
#' @param delta Fraction of time step when the index samples the population. Vector by `i`
#' @return Index at age. Array `[a, i, s]`
#' @details
#' The index is calculated as
#' \deqn{
#' I_{a,i,s} = v_{a,i,s} \sum_r N_{a,r,s} \exp(-\delta_i Z_{a,r,s}) \times \mathbb{1}(r \in R_i) \mathbb{1}(s \in S_i)
#' }
#'
#' where \eqn{R_i} and \eqn{S_i} denote the regions and stocks, respectively, sampled by index \eqn{i}. For example,
#' \eqn{R_2 = 1} denotes that the second index of abundance only samples region 1. These are informed by array `samp` where
#' `samp[i, r, s] = 1` indicates that stock `s` in region `r` is sampled by index `i`.
#' @export
calc_index <- function(N, Z, sel, na = dim(N)[1], ni = dim(sel)[2], ns = dim(N)[3],
                       samp = array(1, c(ni, na, ns)), delta = rep(0, ni)) {

  N_ais <- sapply(1:ns, function(s) {
    sapply(1:ni, function(i) {
      r_i <- samp[i, , s]
      if (sum(r_i)) {
        N_ars <- N[, r_i, s, drop = FALSE] * exp(-delta[i] * Z[, r_i, s, drop = FALSE])
        N_a <- apply(N_ars, 1, sum)
      } else {
        N_a <- numeric(na)
      }
      return(N_a)
    })
  })
  IN_ais <- N_ais * sel
  return(IN_ais)
}

calc_q <- function(Iobs, B) {
  i <- !is.na(Iobs) & Iobs > 0
  n <- sum(i)
  num <- log(Iobs[i]/B[i]) %>% sum()
  q <- exp(num/n)
  return(q)
}

#' Calculate movement matrix for all age classes
#'
#' Movement matrices are calculated for all age classes from a base matrix and a gravity model formulation
#' (Carruthers et al. 2016).
#'
#' @param x Base log-movement parameters. See details. Array `[a, r, r]`
#' @param g Gravity model attractivity term. Tendency to move to region `r`. Matrix `[a, r]`
#' @param v Gravity model viscosity term. Tendency to stay in same region. Vector by `a`
#' @param na Integer, number of ages
#' @param nr Integer, number of regions
#' @param aref Integer, reference age class
#' @details
#' Rows index region of origin and columns denote region of destination.
#'
#' In log space, the movement matrix \eqn{m_a} for age class \eqn{a} from region \eqn{r} to \eqn{r'} is the sum of base matrix \eqn{x} and
#' gravity matrix \eqn{G}:
#' \deqn{m_{a,r,r'} = x_{a,r,r'} + G_{a,r,r'}}
#'
#' To essentially exclude movement from \eqn{r} to \eqn{r'}, set \eqn{x_{a,r,r'} = -1000}.
#'
#' Gravity matrix \eqn{G} includes an attractivity term \eqn{g} and viscosity term \eqn{v}:
#'
#' \deqn{G_{a,r,r'} =
#' \begin{cases}
#' g'_{a,r'} + v_a \quad & r = r'\\
#' g'_{a,r'} \quad & \textrm{otherwise}
#' \end{cases}
#' }
#'
#' Vector \eqn{g'} are offset terms relative to the value for the reference age class:
#' \deqn{g'_{a,r'} =
#' \begin{cases}
#' g_{a,r} \quad & a = a_{ref}\\
#' g_{a,r} + g_{a=aref,r} \quad & \textrm{otherwise}
#' \end{cases}
#' }
#'
#' The movement matrix in normal space is obtained by the softmax transformation:
#' \deqn{M_{a,r,r'} = \dfrac{\exp(m_{a,r,r'})}{\sum_{r'}\exp(m_{a,r,r'})}}
#'
#' If \eqn{x} and \eqn{v} are zero, then the movement matrix simply distributes the total stock
#' abundance into the various regions as specified in \eqn{g'}.
#' @references
#' Carruthers, T.R., et al. 2015. Modelling age-dependent movement: an application to red and
#' gag groupers in the Gulf of Mexico. CJFAS 72: 1159-1176. \doi{10.1139/cjfas-2014-0471}
#' @return Movement array `[a, r, r]`
#' @export
conv_mov <- function(x, g, v, na = dim(x)[1], nr = dim(x)[2], aref = ceiling(0.5 * na)) {
  stopifnot(length(g) == na)
  stopifnot(length(g) == length(v))
  x <- array(x, c(na, nr, nr))
  g <- matrix(g, na, nr)

  mov_rra <- sapply2(1:na, function(a) {
    gg <- g[a, ] + (a != aref) * g[aref, ]
    ln_mov <- x[a, , ] + v[a] * diag(nr) + matrix(gg, nr, nr, byrow = TRUE)
    apply(ln_mov, 1, softmax) %>% t()
  }) %>%
    array(c(nr, nr, na))

  return(aperm(mov_rra, c(3, 1, 2)))
}


#' Predict the probability of CKMR kinship pairs
#'
#' Calculate the probability of observing a parent-offspring pair (`calc_POP`) and
#' half-sibling pair (`calc_HSP`) for closed-kin mark recapture (CKMR) for an age-structured
#' model.
#'
#' @param t Vector, capture year of parent `i`
#' @param a Vector, age at capture of parent `i`
#' @param y Vector, birth year of offspring `j`
#' @param N Abundance of mature spawners. Matrix by `[y, a]`
#' @param fec Fecundity schedule of mature spawners. Matrix by `[y, a]`
#' @return A vector of probabilities.
#' @seealso [like_CKMR()]
#' @section Parent-offspring pairs:
#' The parent-offspring probability is calculated from Bravington et al. 2016, eq 3.4:
#'
#' \deqn{p_{\textrm{POP}} = 2 \times \dfrac{f(y_j,y_j - (t_i - a_i))}{\sum_a f(y_j,a) N(y_j,a)}}
#'
#' where \eqn{y_j - (t_i - a_i)} is the parental age in year \eqn{y_j}. Scalar 2 accounts for the fact
#' that the parent could be either a mother or a father.
#' `calc_POP` is vectorized with respect to `t`, `a`, and `y`.
#' @references
#' Bravington, M.V. et al. 2016. Close-Kin Mark-Recapture. Stat. Sci. 31: 259-274.
#' \doi{10.1214/16-STS552}
#' @author Q. Huynh with contribution from Y. Tsukahara (Fisheries Research Institute, Japan)
#' @export
calc_POP <- function(t, a, y, N, fec) {
  if (is.null(t)) return(0)
  a_yj <- y - t + a # Age of parent in birth year of offspring, Bravington 2016, eq 3.4
  if (!length(a_yj) == length(y)) stop("Vectors t, a, y need to be the same length")

  rel_RO <- sapply(1:length(a_yj), function(j) {
    RO <- sum(N[y[j], ] * fec[y[j], ])
    fec[y[j], a_yj[j]]/RO
  })
  p <- 2 * rel_RO
  return(p)
}

#' @name calc_POP
#' @section Half-sibling pairs:
#' The half-sibling probability is calculated from Bravington et al. 2016, eq 3.10, and expanded
#' by Hillary et al. 2018, Supplement S2.8.1 for
#' age-specific survival and fecundity of the parent:
#'
#' \deqn{p_{\textrm{HSP}} = \sum_a\left(
#' \dfrac{N(y_i, a)f(y_i, a)}{\sum_{a'} N(y_i, a')f(y_i,a')}
#' \exp(-\sum_{t = 0}^{y_j - y_i - 1} Z(y_i + t,a' + t))
#' \dfrac{f(y_j,a+y_j-y_i)}{\sum_{a'} N(y_j,a')f(y_i,a')}
#' \right)
#' }
#'
#' - The first ratio is the probability that a fish at age \eqn{a} in year \eqn{y_i} is the parent of \eqn{i}.
#' - The exponential term is that fish's survival from year \eqn{y_i} to \eqn{y_j}.
#' - The second ratio is the probability that the parent of \eqn{i}, age \eqn{a+y_j-y_i} in year \eqn{y_j}, is the parent of \eqn{j}.
#'
#' The parent is not observed in the HSP, so we sum the probabilities over all potential ages in year \eqn{y_i}.
#' `calc_HSP` is vectorized with respect to `yi` and `yj`.
#' @references
#' Hillary, R.M. et al. 2018. Genetic relatedness reveals total population size of white sharks in eastern Australia and
#' New Zealand. Sci. Rep. 8: 2661. \doi{10.1038/s41598-018-20593-w}
#' @param yi Vector, birth year of sibling `i`. Must be older than sibling `j`.
#' @param yj Vector, birth year of sibling `j`.
#' @param Z Instantaneous total mortality rate. Matrix by `[y, a]`
#' @export
calc_HSP <- function(yi, yj, N, fec, Z) {
  if (is.null(yi)) return(0)
  stopifnot(length(yi) == length(yj))
  sapply(1:length(yi), function(x) {
    .calc_HSP(yi[x], yj[x], N = N, fec = fec, Z = Z)
  })
}

.calc_HSP <- function(yi, yj, N, fec, Z, na = ncol(N)) {
  delta_ij <- yj - yi

  RO_yi <- sum(N[yi, ] * fec[yi, ])
  RO_yj <- sum(N[yj, ] * fec[yj, ])

  p_ai <- sapply(1:na, function(ai) {
    relRO_di <- fec[yi, ai]/RO_yi
    relRO_dj <- fec[yj, min(ai + delta_ij, na)]/RO_yj

    Z_ij <- sapply(seq(0, delta_ij - 1), function(t) Z[yi + t, min(ai + t, na)])
    surv_dij <- exp(-sum(Z_ij))

    N[yi, ai] * relRO_di * relRO_dj * surv_dij
  })
  p <- sum(p_ai)
  return(p)
}
