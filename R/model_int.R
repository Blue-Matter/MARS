
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
#' @param q_fs Relative catchability of stock `s` for fleet `f`. Defaults to 1 if missing. Matrix `[f, s]`
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
  if (missing(q_fs)) q_fs <- 1

  Cobs <- matrix(Cobs, nf, nr)
  N <- array(N, c(na, nr, ns))
  sel <- array(sel, c(na, nf, ns))
  wt <- array(wt, c(na, nf, ns))
  M <- matrix(M, na, ns)
  q_fs <- matrix(q_fs, nf, ns)

  # Initialize search ----
  fn <- gr <- x_loop <- list()

  ind_afrs <- as.matrix(expand.grid(a = 1:na, f = 1:nf, r = 1:nr, s = 1:ns))
  ind_ars <- as.matrix(expand.grid(a = 1:na, r = 1:nr, s = 1:ns))

  frs_afrs <- ind_afrs[, c("f", "r", "s")]
  afs_afrs <- ind_afrs[, c("a", "f", "s")]
  ars_afrs <- ind_afrs[, c("a", "r", "s")]
  fr_afrs <- ind_afrs[, c("f", "r")]
  fs_afrs <- ind_afrs[, c("f", "s")]

  as_ars <- ind_ars[, c("a", "s")]

  VB_afrs <- array(N[ars_afrs] * sel[afs_afrs] * wt[afs_afrs], c(na, nf, nr, ns))
  VB_fr <- apply(VB_afrs, c(2, 3), sum)
  if (inherits(Cobs, "advector")) {
    Cobs_loop <- CondExpLt(Cobs, 1e-8, 1e-9, Cobs)
  } else {
    Cobs_loop <- ifelse(Cobs < 1e-8, 1e-9, Cobs)
  }
  F_init <- Cobs_loop/(Cobs_loop + VB_fr)

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
      F_loop <- CondExpGt(x_loop[[i]], ln_Fmax, Fmax, exp(x_loop[[i]]))
      if (i == nitF + 1) { # Last iteration used to calculate corresponding f and g with F and penalty
        F_loop <- CondExpLt(Cobs, 1e-8, 0, F_loop)
        penalty <- penalty + sum(posfun(Fmax, F_loop))
      }
    } else {
      F_loop <- Fmax * plogis(x_loop[[i]])
    }

    F_frs <- sapply2(1:ns, function(s) q_fs[, s] * F_loop) %>% array(c(nf, nr, ns))
    F_afrs <- array(F_frs[frs_afrs] * sel[afs_afrs], c(na, nf, nr, ns))
    F_ars <- apply(F_afrs, c(1, 3, 4), sum)
    Z_ars <- array(F_ars + delta * M[as_ars], c(na, nr, ns))
    .gamma_ars <- 1 - exp(-Z_ars)

    CN_afrs <- array(F_afrs * .gamma_ars[ars_afrs] * N[ars_afrs]/ Z_ars[ars_afrs], c(na, nf, nr, ns))
    CB_afrs <- array(CN_afrs * wt[afs_afrs], c(na, nf, nr, ns))
    CB_fr <- apply(CB_afrs, 2:3, sum)

    fn[[i]] <- CB_fr - Cobs

    if (trans == "log") {
      deriv_F <- F_loop # f r
    } else {
      deriv_F <- Fmax * exp(-x_loop[[i]])/(1 + exp(-x_loop[[i]]))/(1 + exp(-x_loop[[i]]))
    }

    deriv_Z_afrs <- array(q_fs[fs_afrs] * deriv_F[fr_afrs] * sel[afs_afrs], c(na, nf, nr, ns))
    deriv_gamma_afrs <- array(exp(-Z_ars[ars_afrs]) * deriv_Z_afrs, c(na, nf, nr, ns))

    constants_afrs <- array(q_fs[fs_afrs] * sel[afs_afrs] * N[ars_afrs] * wt[afs_afrs], c(na, nf, nr, ns))

    deriv1_afrs <- array(
      deriv_F[fr_afrs] * .gamma_ars[ars_afrs] + F_loop[fr_afrs] * deriv_gamma_afrs,
      c(na, nf, nr, ns)
    )
    deriv2_afrs <- array(deriv1_afrs * Z_ars[ars_afrs], c(na, nf, nr, ns))
    deriv3_afrs <- array(
      deriv2_afrs - .gamma_ars[ars_afrs] * F_loop[fr_afrs] * deriv_Z_afrs,
      c(na, nf, nr, ns)
    )
    deriv_afrs <- array(deriv3_afrs/Z_ars[ars_afrs]/Z_ars[ars_afrs], c(na, nf, nr, ns))

    gr[[i]] <- apply(constants_afrs * deriv_afrs, c(2, 3), sum)

    if (i <= nitF) x_loop[[i+1]] <- x_loop[[i]] - fn[[i]]/gr[[i]]
  }
  #browser(expr = any(Cobs >= 1e-8))
  CB_frs <- apply(CB_afrs, 2:4, sum)
  list(F_afrs = F_afrs, F_ars = F_ars, F_index = F_loop, Z_ars = Z_ars,
       CB_frs = CB_frs, CN_afrs = CN_afrs, VB_afrs = VB_afrs,
       penalty = penalty, fn = fn[[nitF + 1]], gr = gr[[nitF + 1]])
}

calc_Baranov <- function(FM, Z, N) FM/Z * (1 - exp(-Z)) * N


#' @importFrom stats uniroot
.calc_summary_F <- function(FM, M, N, CN) calc_Baranov(FM, FM + M, N) - CN
calc_summary_F <- function(M, N, CN, Fmax) {

  if (is.na(CN) || is.na(N)) {
    out <- NA_real_
  } else if (CN > N) {
    out <- Inf
  } else {
    out <- try(uniroot(.calc_summary_F, interval = c(0, Fmax), M = M, N = N, CN = CN)$root, silent = TRUE)
    if (is.character(out)) out <- Fmax
  }
  return(out)
}

#' Calculate recruitment from stock-recruit function
#'
#' @param x Numeric, either the spawning output or the equilibrium spawners per recruit, from which
#' the recruitment will be calculated. See argument `eq`.
#' @param SRR Character to indicate the functional form of the stock recruit function
#' @param eq Logical, indicates whether `x` is the spawning output (`FALSE`) or equilibrium spawners per recruit (`TRUE`)
#' @param ... Parameters of the SRR function. Provide one of two sets of variables:
#' 1. `h`, `R0` and `phi0`, or
#' 2. `a` and `b` (alpha, beta values)
#'
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
    val <- if (eq) (a*x - 1)/b/x else a*x/(1 + b*x)
  } else {
    val <- if (eq) log(a*x)/b/x else a*x*exp(-b*x)
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

#' Calculate von Bertalanffy length-at-age
#'
#' Returns an array of length-at-age with seasonal dimension.
#' Useful for [MARSdata-class] inputs.
#'
#' @param Linf_s Vector by stock `s` of asymptotic length
#' @param K_s Vector by stock `s` of the growth coefficient
#' @param t0_s Vector by `s` of the age at length zero.
#' @param ns Integer, number of stocks
#' @param nm Integer, number of seasons
#' @param ny Integer, number of years
#' @param a Integer vector of ages
#' @return Array `[y, m, a, s]`
#' @examples
#' len_ymas <- calc_growth(c(30, 40), c(0.4, 0.2), c(-1, -1))
#'
#' # Calculate stock weight at age
#' a_s <- rep(1e-6, 2)
#' b_s <- c(3, 3.1)
#'
#' ns <- length(a_s)
#' swt_ymas <- sapply(1:ns, function(s) {
#'   a_s[s]*len_ymas[, , , s]^b_s[s]
#' }, simplify = "array")
#' @export
calc_growth <- function(Linf_s, K_s, t0_s, ns = length(Linf_s), nm = 4, ny = 20, a = seq(1, 10)) {
  len_ymas <- sapply2(1:ns, function(s) {
    sapply2(a, function(aa) {
      sapply(1:nm, function(m) {
        tt <- aa + (m - 1)/nm
        len_y <- Linf_s[s] * (1 - exp(-K_s[s] * (tt - t0_s[s])))
        rep(len_y, ny)
      })
    })
  })
  return(len_ymas)
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
#' @param R Incoming total recruitment. Vector length `s`. Only assigned if `advance_age = TRUE`.
#' @param recdist Distribution of incoming recruitment. Matrix `[r, s]`. Only assigned if `advance_age = TRUE`.
#' @param mov Movement array in the next time step. Array `[a, r, r, s]`. Rows denote region of origin and columns denote region of destination.
#' @param plusgroup Logical, whether the last age class is an accumulator plus group.
#' @return Abundance at the next time step. Array `[a, r, s]`
#' @export
calc_nextN <- function(N, surv, na = dim(N)[1], nr = dim(N)[2], ns = dim(N)[3],
                       advance_age = TRUE, R = numeric(ns),
                       mov = array(1/nr, c(na, nr, nr, ns)),
                       recdist = matrix(1/nr, nr, ns),
                       plusgroup = TRUE) {

  N <- array(N, c(na, nr, ns))
  surv <- array(surv, c(na, nr, ns))
  mov <- array(mov, c(na, nr, nr, ns))

  # Apply survival and advance age class ----
  if (advance_age) {
    is_ad <- inherits(N, "advector") || inherits(surv, "advector") || inherits(R, "advector")
    if (is_ad) {
      `[<-` <- RTMB::ADoverload("[<-")
    }
    Nsurv_ars <- array(0, c(na, nr, ns))
    Nsurv_ars[2:na, , ] <- N[2:na - 1, , ] * surv[2:na - 1, , ]
    if (plusgroup) Nsurv_ars[na, , ] <- Nsurv_ars[na, , ] + N[na, , ] * surv[na, , ]
    Nsurv_ars[1, , ] <- sapply(1:ns, function(s) recdist[, s] * R[s])
  } else {
    Nsurv_ars <- N * surv
  }

  # Distribute stock ----
  if (nr > 1) {
    ind_arrs <- as.matrix(expand.grid(a = 1:na, rf = 1:nr, rt = 1:nr, s = 1:ns))
    arfs_arrs <- ind_arrs[, c("a", "rf", "s")]

    Nnext_arrs <- array(Nsurv_ars[arfs_arrs] * mov, c(na, nr, nr, ns))
    Nnext_ars <- apply(Nnext_arrs, c(1, 3, 4), sum)
  } else {
    Nnext_ars <- Nsurv_ars
  }

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
#' @param nr Integer, number of regions
#' @param ns Integer, number of stocks
#' @param ni Integer, number of indices
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
calc_index <- function(N, Z, sel, na = dim(N)[1], nr = dim(N)[2], ns = dim(N)[3], ni = dim(sel)[2],
                       samp = array(1, c(ni, nr, ns)), delta = rep(0, ni)) {

  N <- array(N, c(na, nr, ns))
  Z <- array(Z, c(na, nr, ns))
  sel <- array(sel, c(na, ni, ns))

  ind_airs <- as.matrix(expand.grid(a = 1:na, i = 1:ni, r = 1:nr, s = 1:ns))
  irs_airs <- ind_airs[, c("i", "r", "s")]
  ars_airs <- ind_airs[, c("a", "r", "s")]
  ais_airs <- ind_airs[, c("a", "i", "s")]
  i_airs <- ind_airs[, "i"]

  IN_airs <- array(
    N[ars_airs] * samp[irs_airs] * sel[ais_airs] * exp(-delta[i_airs] * Z[ars_airs]),
    c(na, ni, nr, ns)
  )
  IN_ais <- apply(IN_airs, c(1, 2, 4), sum)

  #N_ais <- sapply2(1:ns, function(s) {
  #  sapply(1:ni, function(i) {
  #    r_i <- samp[i, , s]
  #    if (sum(r_i)) {
  #      N_ars <- N[, r_i, s, drop = FALSE] * exp(-delta[i] * Z[, r_i, s, drop = FALSE])
  #      N_a <- apply(N_ars, 1, sum)
  #    } else {
  #      N_a <- numeric(na)
  #    }
  #    return(N_a)
  #  })
  #})
  #IN_ais2 <- N_ais * sel
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
  x <- array(x, c(na, nr, nr))
  g <- matrix(g, na, nr)

  mov_rra <- sapply2(1:na, function(a) {
    gg <- g[a, ] + (a != aref) * g[aref, ]
    ln_mov <- x[a, , ] + v[a] * diag(nr) + matrix(gg, nr, nr, byrow = TRUE)
    matrix(ln_mov, nr, nr) %>%
      apply(1, softmax) # This operation flips matrix orientation. FROM = col, TO = row
  }) %>%
    array(c(nr, nr, na))

  return(aperm(mov_rra, 3:1))
}

#' Equilibrium distribution from movement matrix
#'
#' Applies the movement matrix several times in order to obtain the equilibrium spatial distribution of a movement matrix.
#' Not used in the model but useful for reporting.
#' @return Numeric vector of length `nr`
#' @param x Movement matrix, a square matrix with rows corresponding to origin (sum to 1), and columns corresponding to destination
#' @param nr Number of regions
#' @param start The initial distribution. Vector of length `nr`
#' @param nit Integer, the number of times the movement matrix will be applied
#' @export
calc_eqdist <- function(x, nr = dim(x)[2], start = rep(1/nr, nr), nit = 20) {
  if (inherits(start, "advector") || inherits(x, "advector")) {
    `[<-` <- RTMB::ADoverload("[<-")
  }

  N <- array(NA_real_, c(nit, nr))
  N[1, ] <- start
  for (i in 2:nit - 1) N[i+1, ] <- colSums(N[i, ] * x)
  return(N[nit-1, ])
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
#' by Hillary et al. 2018, Supplement S2.8.1 for age-specific survival and fecundity of the parent:
#'
#' \deqn{p_{\textrm{HSP}} = 4 \times
#' \sum_a\left(
#' \dfrac{N(y_i, a)f(y_i, a)}{\sum_{a'} N(y_i, a')f(y_i,a')}\times
#' \exp(-\sum_{t = 0}^{y_j - y_i - 1} Z(y_i + t,a + t))\times
#' \dfrac{f(y_j,a+y_j-y_i)}{\sum_{a'} N(y_j,a')f(y_j,a')}
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
  p <- 4 * sum(p_ai)
  return(p)
}
