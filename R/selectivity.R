
#' Selectivity at age and length
#'
#' @description
#' A series of functions to calculate selectivity at age and length from a matrix of parameters.
#'
#' \itemize{
#' \item [conv_selpar()] converts parameters from log or logit space to real units.
#' \item [calc_sel_len()] calculates selectivity at length.
#' \item [calc_fsel_age()] calculates selectivity at age for fisheries, and coordinates dummy fleets.
#' \item [calc_isel_age()] calculates selectivity at age for indices, and can map selectivity from fisheries
#' or population parameters (e.g, mature biomass).
#' }
#' @param x Estimated parameters. Matrix `[3, f]`
#' @param type Character string to indicate the functional form of selectivity. Options include: "logistic_length", "dome_length",
#' "logistic_age", "dome_age", an integer (`f`) to map index selectivity to the corresponding fleet `f` (will be coerced to integer),
#' or "SB" to fix to maturity at age schedule.
#' @param maxage Maximum value of the age of full selectivity
#' @param maxL Maximum value of the length of full selectivity
#' @section Converting selectivity parameters (conv_selpar):
#' The first row of `x` corresponds to the length or age of full selectivity: \eqn{\mu_f = L_{max}/(1 + \exp(-x_{1,f}))}
#'
#' The second row of `x` corresponds to the width of the ascending limb for selectivity:
#' \eqn{\sigma_f^{asc} = \exp(x_{2,f})}
#'
#' The third row of `x` corresponds to the width of the descending limb for selectivity (if dome-shaped):
#' \eqn{\sigma_f^{des} = \exp(x_{3,f})}
#' @return
#' [conv_selpar()] returns a matrix of the same dimensions as `x`.
#' @export
conv_selpar <- function(x, type, maxage, maxL) {
  nf <- length(type)
  stopifnot(ncol(x) == nf)
  stopifnot(nrow(x) >= 3)

  sel_par <- sapply(1:nf, function(f) {
    sd_asc <- exp(x[2, f])
    sd_desc <- exp(x[3, f])
    if (grepl("age", type[f])) {
      Aapical <- maxage * plogis(x[1, f])
      v <- c(Aapical, sd_asc, sd_desc)
    } else if (grepl("length", type[f])) {
      Lapical <- maxL * plogis(x[1, f])
      v <- c(Lapical, sd_asc, sd_desc)
    } else {
      v <- rep(NA_real_, 3)
    }
    return(v)
  })

  return(sel_par)
}

#' @rdname conv_selpar
#' @aliases calc_sel_len
#' @param sel_par Matrix of parameters returned by [conv_selpar()]
#' @param lmid Midpoints of length bins for calculating selectivity at length
#' @section Length selectivity (calc_sel_len):
#' The selectivity at length is
#' \deqn{
#' s_{\ell} =
#' \begin{cases}
#' 2^\alpha & L_{\ell} < \mu_f\\
#' 2^\beta & L_{\ell} \ge \mu_f\\
#' \end{cases}
#' }
#' where
#' \eqn{
#' \alpha = -(L_\ell - \mu_f)^2/(\sigma_f^{asc})^2
#' }
#' and
#' \eqn{
#' \beta = -(L_\ell - \mu_f)^2/(\sigma_f^{des})^2
#' }
#' @return
#' [calc_sel_len()] returns a matrix `[l, f]`, i.e., `[length(lmid), length(type)]`.
#' @export
calc_sel_len <- function(sel_par, lmid, type) {
  nf <- length(type)

  sel_lf <- sapply(1:nf, function(f) {
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

  return(sel_lf)
}

#' @rdname conv_selpar
#' @param sel_len Selectivity at length matrix returned by [calc_sel_len()]
#' @param LAK Length-at-age probability matrix. Matrix `[a, length(lmid)]`
#' @param sel_block Integer vector. Length `length(type)`. See details below.
#' @param mat Maturity at age vector
#' @param a Integer vector of ages corresponding to the rows of `LAK` (as well as `mat`)
#' @section Age selectivity (calc_fsel_age):
#' The equivalent selectivity at age is converted from the length values (`sel_len`) as
#' \deqn{
#' s_a = \sum_\ell s_\ell \times \textrm{Prob}(L_{\ell}|a)
#' }
#'
#' If selectivity is explicitly in age units, then the values are directly calculated
#' from parameters in `sel_par`.
#'
#' Vector `sel_block` assigns the output selectivity from a different column of the input matrix
#' and facilitates time-varying selectivity in time blocks. For example, `sel_block[1] <- 2` means
#' that selectivity values in the first column of the output is based on the second column of the
#' input matrices (`sel_len[, 2]` or `sel_par[, 2]`).
#' @return
#' [calc_fsel_age()] returns a matrix `[a, f]`, i.e., `[a, length(sel_block)]`
#' @export
calc_fsel_age <- function(sel_len, LAK, type, sel_par, sel_block = seq(1, length(type)), mat, a = seq(1, nrow(LAK))) {
  nf <- length(sel_block)

  sel_af <- sapply(1:nf, function(ff) {
    f <- sel_block[ff]
    if (grepl("length", type[f])) {
      v <- sel_len[, f] %*% t(LAK)
    } else if (grepl("age", type[f])){

      ex_asc <- (a - sel_par[1, f])/sel_par[2, f]
      ex2_asc <- -1 * ex_asc^2
      asc <- 2^ex2_asc

      if (grepl("logistic", type[f])) {
        desc <- 1
      } else {
        ex_desc <- (a - sel_par[1, f])/sel_par[3, f]
        ex2_desc <- -1 * ex_desc^2
        desc <- 2^ex2_desc
      }
      v <- CondExpLt(a, sel_par[1, f], asc, desc)
    } else if (type[f] == "SB") {
      v <- mat
    } else if (type[f] == "free") {
      v <- plogis(sel_par[, f])
    }
    return(v)
  })

  return(sel_af)
}

#' @rdname conv_selpar
#' @param fsel_age Matrix returned by [calc_fsel_age()]
#' @return
#' [calc_isel_age()] returns a matrix `[a, i]`, i.e., `[a, length(type)]`
#' @export
calc_isel_age <- function(sel_len, LAK, type, sel_par, fsel_age, maxage, mat, a = seq(1, nrow(LAK))) {

  old_warn <- options()$warn
  options(warn = -1)
  on.exit(options(warn = old_warn))

  ni <- length(type)

  sel_ai <- sapply(1:ni, function(i) {
    ti <- type[i]
    tii <- as.integer(type[i])

    if (is.na(tii)) {

      if (grepl("length", ti)) {
        v <- sel_len[, i] %*% t(LAK)
      } else if (grepl("age", ti)) {

        ex_asc <- (a - sel_par[1, i])/sel_par[2, i]
        ex2_asc <- -1 * ex_asc^2
        asc <- 2^ex2_asc

        if (grepl("logistic", ti)) {
          desc <- 1
        } else {
          ex_desc <- (a - sel_par[1, i])/sel_par[3, i]
          ex2_desc <- -1 * ex_desc^2
          desc <- 2^ex2_desc
        }
        v <- CondExpLt(a, sel_par[1, i], asc, desc)
      } else if (ti == "SB") {
        v <- mat
      } else if (ti == "free") {
        v <- plogis(sel_par[, i])
      }

    } else {
      v <- fsel_age[, tii]
    }
    return(v)
  })
  return(sel_ai)
}
