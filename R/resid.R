

#' Calculate model residuals
#'
#' Extract residuals from fitted model
#'
#' @param object [MARSassess-class] object returned by `fit_MARS()`
#' @param vars Character vector to indicate which residuals will be calculated.
#' Available choices from [MARSdata-class] object are:
#' "Cinit_mfr", "Cobs_ymfr", "CAAobs_ymafr", "CALobs_ymlfr",
#' "Iobs_ymi", "IAAobs_ymai", "IALobs_ymli", "SC_ymafrs"
#' @param type Character, 'response' for the `log(observed/predicted)` values or 'pearson' for calculating Z-scores.
#' Composition data always use 'pearson'.
#' @param ... Not used
#'
#' @return A named list based on `vars` argument.
#' @aliases residuals
#' @importFrom stats residuals
#' @export
residuals.MARSassess <- function(object, vars, type = c("response", "pearson"), ...) {

  vars_choices <- c("Cinit_mfr", "Cobs_ymfr", "CAAobs_ymafr", "CALobs_ymlfr",
                    "Iobs_ymi", "IAAobs_ymai", "IALobs_ymli", "SC_ymafrs", "tag_ymarrs")
  if (missing(vars)) vars <- vars_choices
  vars <- match.arg(vars, choices = vars_choices, several.ok = TRUE)
  type <- match.arg(type)
  dat <- get_MARSdata(object)

  res <- list()

  # Init catch ----
  if (any(vars == "Cinit_mfr") && !is.null(object@report$initCB_mfrs)) {
    Cinitres <- log(dat@Dfishery@Cinit_mfr/apply(object@report$initCB_mfrs, 1:3, sum))
    if (type == "pearson") {
      Cinitres <- Cinitres/0.01
    }
    res$Cinit_mfr <- Cinitres
  }

  # Catch ----
  if (any(vars == "Cobs_ymfr")) {
    Cres <- log(dat@Dfishery@Cobs_ymfr/apply(object@report$CB_ymfrs, 1:4, sum))
    if (type == "pearson") {
      if (length(dat@Dfishery@Csd_ymfr)) {
        Csd <- dat@Dfishery@Csd_ymfr
      } else {
        Csd <- 0.01
      }
      Cres <- Cres/Csd
    }
    res$Cobs_ymfr <- Cres
  }


  # Catch at age ----
  if (any(vars == "CAAobs_ymafr") && length(dat@Dfishery@CAAobs_ymafr)) {
    res$CAAobs_ymafr <- array(NA, dim(dat@Dfishery@CAAobs_ymafr))
    pred <- apply(object@report$CN_ymafrs, 1:5, sum)
    for(y in 1:dat@Dmodel@ny) {
      for(m in 1:dat@Dmodel@nm) {
        for(f in 1:dat@Dfishery@nf) {
          for(r in 1:dat@Dmodel@nr) {
            pred_i <- pred[y, m, , f, r]
            res$CAAobs_ymafr[y, m, , f, r] <- resid_comp(
              obs = dat@Dfishery@CAAobs_ymafr[y, m, , f, r],
              pred = pred_i,
              like = dat@Dfishery@fcomp_like,
              N = dat@Dfishery@CAAN_ymfr[y, m, f, r],
              theta = dat@Dfishery@CAAtheta_f[f],
              stdev = sum(pred_i)/pred_i
            )
          }
        }
      }
    }
  }

  # Catch at length ----
  if (any(vars == "CALobs_ymlfr") && length(dat@Dfishery@CALobs_ymlfr)) {
    res$CALobs_ymlfr <- array(NA, dim(dat@Dfishery@CALobs_ymlfr))
    pred <- apply(object@report$CN_ymlfrs, 1:5, sum)
    for(y in 1:dat@Dmodel@ny) {
      for(m in 1:dat@Dmodel@nm) {
        for(f in 1:dat@Dfishery@nf) {
          for(r in 1:dat@Dmodel@nr) {
            pred_i <- pred[y, m, , f, r]
            res$CALobs_ymlfr[y, m, , f, r] <- resid_comp(
              obs = dat@Dfishery@CALobs_ymlfr[y, m, , f, r],
              pred = pred_i,
              like = dat@Dfishery@fcomp_like,
              N = dat@Dfishery@CALN_ymfr[y, m, f, r],
              theta = dat@Dfishery@CALtheta_f[f],
              stdev = sum(pred_i)/pred_i
            )
          }
        }
      }
    }
  }

  # Indices ----
  if (any(vars == "Iobs_ymi") && length(dat@Dsurvey@Iobs_ymi)) {
    res$Iobs_ymi <- log(dat@Dsurvey@Iobs_ymi/object@report$I_ymi)
    if (type == "pearson") {
      res$Iobs_ymi <- res$Iobs_ymi/dat@Dsurvey@Isd_ymi
    }
  }

  # Indices at age ----
  if (any(vars == "IAAobs_ymai") && length(dat@Dsurvey@IAAobs_ymai)) {
    res$IAAobs_ymai <- array(NA, dim(dat@Dsurvey@IAAobs_ymai))
    pred <- apply(object@report$IN_ymais, 1:4, sum)
    for(y in 1:dat@Dmodel@ny) {
      for(m in 1:dat@Dmodel@nm) {
        for(i in 1:dat@Dmodel@ni) {
          pred_i <- pred[y, m, , i]
          res$IAAobs_ymai[y, m, , i] <- resid_comp(
            obs = dat@Dsurvey@IAAobs_ymai[y, m, , i],
            pred = pred_i,
            like = dat@Dsurvey@icomp_like,
            N = dat@Dsurvey@IAAN_ymi[y, m, i],
            theta = dat@Dsurvey@IAAtheta_i[i],
            stdev = sum(pred_i)/pred_i
          )
        }
      }
    }
  }

  # Indices at length ----
  if (any(vars == "IALobs_ymli") && length(dat@Dsurvey@IALobs_ymli)) {
    res$IALobs_ymli <- array(NA, dim(dat@Dsurvey@IALobs_ymli))
    pred <- apply(object@report$IN_ymlis, 1:4, sum)
    for(y in 1:dat@Dmodel@ny) {
      for(m in 1:dat@Dmodel@nm) {
        for(i in 1:dat@Dmodel@ni) {
          pred_i <- pred[y, m, , i]
          res$IALobs_ymli[y, m, , i] <- resid_comp(
            obs = dat@Dsurvey@IALobs_ymli[y, m, , i],
            pred = pred_i,
            like = dat@Dsurvey@icomp_like,
            N = dat@Dsurvey@IALN_ymi[y, m, i],
            theta = dat@Dsurvey@IALtheta_i[i],
            stdev = sum(pred_i)/pred_i
          )
        }
      }
    }
  }

  # Stock composition (at age) ----
  if (any(vars == "SC_ymafrs") && length(dat@Dfishery@SC_ymafrs)) {
    SC_ymafrs <- dat@Dfishery@SC_ymafrs
    res$SC_ymafrs <- array(NA, dim(dat@Dfishery@SC_ymafrs))
    pred <- object@report$SCpred_ymafrs
    for(y in 1:dat@Dmodel@ny) {
      for(m in 1:dat@Dmodel@nm) {
        for(a in 1:dim(SC_ymafrs)[3]) {
          for(f in 1:dim(SC_ymafrs)[4]) {
            for(r in 1:dat@Dmodel@nr) {
              res$SC_ymafrs[y, m, a, f, r, ] <- resid_comp(
                obs = SC_ymafrs[y, m, a, f, r, ],
                pred = pred[y, m, a, f, r, ],
                like = dat@Dfishery@SC_like,
                N = dat@Dfishery@SCN_ymafr[y, m, a, f, r],
                theta = dat@Dfishery@SCtheta_f[f],
                stdev = dat@Dfishery@SCstdev_ymafrs[y, m, a, f, r, ]
              )
            }
          }
        }
      }
    }
  }


  # Stock composition (at age) ----
  if (any(vars == "SC_ymafrs") && length(dat@Dfishery@SC_ymafrs)) {
    SC_ymafrs <- dat@Dfishery@SC_ymafrs
    res$SC_ymafrs <- array(NA, dim(SC_ymafrs))
    pred <- object@report$SCpred_ymafrs
    for(y in 1:dat@Dmodel@ny) {
      for(m in 1:dat@Dmodel@nm) {
        for(a in 1:dim(SC_ymafrs)[3]) {
          for(f in 1:dim(SC_ymafrs)[4]) {
            for(r in 1:dat@Dmodel@nr) {
              res$SC_ymafrs[y, m, a, f, r, ] <- resid_comp(
                obs = SC_ymafrs[y, m, a, f, r, ],
                pred = pred[y, m, a, f, r, ],
                like = dat@Dfishery@SC_like,
                N = dat@Dfishery@SCN_ymafr[y, m, a, f, r],
                theta = dat@Dfishery@SCtheta_f[f],
                stdev = dat@Dfishery@SCstdev_ymafrs[y, m, a, f, r, ]
              )
            }
          }
        }
      }
    }
  }

  # Tag movement ----
  if (any(vars == "tag_ymarrs") && length(dat@Dtag@tag_ymarrs)) {
    tag_ymarrs <- dat@Dtag@tag_ymarrs
    res$tag_ymarrs <- array(NA, dim(tag_ymarrs))
    pred <- object@report$tagpred_ymarrs
    for(y in 1:dim(tag_ymarrs)[1]) {
      for(m in 1:dat@Dmodel@nm) {
        for(a in 1:dim(tag_ymarrs)[3]) {
          for(rf in 1:dat@Dmodel@nr) {
            for(s in 1:dat@Dmodel@ns) {
              res$tag_ymarrs[y, m, a, rf, , s] <- resid_comp(
                obs = tag_ymarrs[y, m, a, rf, , s],
                pred = pred[y, m, a, rf, , s],
                like = dat@Dtag@tag_like,
                N = dat@Dtag@tagN_ymars[y, m, a, rf, s],
                theta = dat@Dtag@tagtheta_s[s],
                stdev = dat@Dtag@tagstdev_s[s]
              )
            }
          }
        }
      }
    }
  }

  return(res)
}

# Vectors
resid_comp <- function(obs, pred, like, ...) {

  dots <- list(...)

  obs_prob <- obs/sum(obs, na.rm = TRUE)
  pred_prob <- pred/sum(pred, na.rm = TRUE)

  # Observed minus predicted
  num <- switch(
    like,
    "multinomial" = dots$N * (obs_prob - pred_prob),
    "dirmult1" = dots$N * (obs_prob - pred_prob),
    "dirmult2" = dots$N * (obs_prob - pred_prob),
    "lognormal" = ifelse(obs > 0, log(obs_prob/pred_prob), NA),
    NA
  )

  # Variance
  denom <- switch(
    like,
    "multinomial" = dots$N * pred_prob * (1 - pred_prob),
    "dirmult1" = dots$N * pred_prob * (1 - pred_prob) * dots$N * (1 + dots$theta) / (1 + dots$theta * dots$N),
    "dirmult2" = dots$N * pred_prob * (1 - pred_prob) * (dots$N + dots$theta) / (1 + dots$theta),
    "lognormal" = dots$stdev * dots$stdev,
    NA
  )

  res <- num/sqrt(denom)
  return(res)
}

plot_resid_Cobs <- function(fit, f = 1, ...) {
  vars <- "Cobs_ymfr"

  dat <- get_MARSdata(fit)

  res <- residuals(fit, vars = vars, ...)
  x <- apply(res[[vars]][, , f, , drop = FALSE], c(1, 2, 4), identity)

  name <- dat@Dlabel@region

  year <- make_yearseason(dat@Dlabel@year, dat@Dmodel@nm)
  x <- collapse_yearseason(x)

  color <- make_color(ncol(x), type = "region")
  fname <- dat@Dlabel@fleet[f]

  matplot(year, x, xlab = "Year", ylab = paste(fname, "catch residual"), typ = "o", col = color, pch = 16)
  abline(h = 0, lty = 2, col = "grey60")
  if (ncol(x) > 1) legend("topleft", legend = name, col = color, lwd = 1, pch = 16, horiz = TRUE)

  invisible()
}

plot_resid_Iobs <- function(fit, i = 1, ...) {
  vars <- "Iobs_ymi"

  dat <- get_MARSdata(fit)

  res <- residuals(fit, vars = vars, ...)[[vars]]
  if (is.null(res)) return(invisible())

  x <- apply(res[, , i, drop = FALSE], 1:2, identity)

  name <- dat@Dlabel@index[i]

  year <- make_yearseason(dat@Dlabel@year, dat@Dmodel@nm)
  x <- collapse_yearseason(x)

  plot(year, x, typ = "n", xlab = "Year", ylab = paste(name, "residual"))
  lines(year[!is.na(x)], x[!is.na(x)], col = "grey70", lty = 2)
  points(year, x, typ = "o", pch = 16)
  abline(h = 0, lty = 1, col = "grey60")

  invisible()
}

plot_resid_CAA <- function(fit, f = 1, r = 1, do_hist = FALSE, ...) {
  vars <- "CAAobs_ymafr"

  dat <- get_MARSdata(fit)

  res <- residuals(fit, vars = vars, ...)[[vars]]
  if (is.null(res)) return(invisible())

  x <- apply(res[, , , f, r, drop = FALSE], 1:3, identity)

  year <- make_yearseason(dat@Dlabel@year, dat@Dmodel@nm)
  x <- collapse_yearseason(x)

  ind <- rowSums(x, na.rm = TRUE) != 0
  yeardiff <- c(diff(year), diff(year)[1])

  .plot_resid_comp(
    year[ind], dat@Dlabel@age, x[ind, , drop = FALSE],
    xlab = "Year", ylab = "Age", do_hist = do_hist,
    xdiff = yeardiff[ind]
  )
}

plot_resid_CAL <- function(fit, f = 1, r = 1, do_hist = FALSE, ...) {
  vars <- "CALobs_ymlfr"

  dat <- get_MARSdata(fit)

  res <- residuals(fit, vars = vars, ...)[[vars]]
  if (is.null(res)) return(invisible())

  x <- apply(res[, , , f, r, drop = FALSE], 1:3, identity)

  year <- make_yearseason(dat@Dlabel@year, dat@Dmodel@nm)
  x <- collapse_yearseason(x)

  ind <- rowSums(x, na.rm = TRUE) != 0
  yeardiff <- c(diff(year), diff(year)[1])

  .plot_resid_comp(
    year[ind], dat@Dmodel@lmid, x[ind, , drop = FALSE],
    xlab = "Year", ylab = "Length", do_hist = do_hist,
    xdiff = yeardiff[ind]
  )
}

plot_resid_IAA <- function(fit, i = 1, do_hist = FALSE, ...) {
  vars <- "IAAobs_ymai"

  dat <- get_MARSdata(fit)

  res <- residuals(fit, vars = vars, ...)[[vars]]
  if (is.null(res)) return(invisible())

  x <- apply(res[, , , i, drop = FALSE], 1:3, identity)

  year <- make_yearseason(dat@Dlabel@year, dat@Dmodel@nm)
  x <- collapse_yearseason(x)

  ind <- rowSums(x, na.rm = TRUE) != 0
  yeardiff <- c(diff(year), diff(year)[1])

  .plot_resid_comp(
    year[ind], dat@Dlabel@age, x[ind, , drop = FALSE],
    xlab = "Year", ylab = "Age", do_hist = do_hist,
    xdiff = yeardiff[ind]
  )
}

plot_resid_IAL <- function(fit, i = 1, do_hist = FALSE, ...) {
  vars <- "IALobs_ymli"

  dat <- get_MARSdata(fit)

  res <- residuals(fit, vars = vars, ...)[[vars]]
  if (is.null(res)) return(invisible())

  x <- apply(res[, , , i, drop = FALSE], 1:3, identity)

  year <- make_yearseason(dat@Dlabel@year, dat@Dmodel@nm)
  x <- collapse_yearseason(x)

  ind <- rowSums(x, na.rm = TRUE) != 0
  yeardiff <- c(diff(year), diff(year)[1])

  .plot_resid_comp(
    year[ind], dat@Dmodel@lmid, x[ind, , drop = FALSE],
    xlab = "Year", ylab = "Length", do_hist = do_hist,
    xdiff = yeardiff[ind]
  )
}

plot_resid_SC <- function(fit, a = 1, f = 1, r = 1, do_hist = FALSE, ...) {
  vars <- "SC_ymafrs"
  dat <- get_MARSdata(fit)

  res <- residuals(fit, vars = vars, ...)[[vars]]
  if (is.null(res)) return(invisible())

  x <- apply(res[, , a, f, r, , drop = FALSE], c(1, 2, 6), identity)

  year <- make_yearseason(dat@Dlabel@year, dat@Dmodel@nm)
  x <- collapse_yearseason(x)

  ind <- rowSums(x, na.rm = TRUE) != 0
  yeardiff <- c(diff(year), diff(year)[1])

  .plot_resid_comp(
    year[ind], dat@Dlabel@stock, x[ind, , drop = FALSE],
    xlab = "Year", ylab = "Stock", do_hist = do_hist,
    xdiff = yeardiff[ind]
  )
}

plot_resid_tagmov <- function(fit, yy = 1, aa = 1, s = 1, ...) {
  vars <- "tag_ymarrs"
  dat <- get_MARSdata(fit)

  res <- residuals(fit, vars = vars, ...)[[vars]]
  if (is.null(res)) return(invisible())

  x <- apply(res[yy, , aa, , , s, drop = FALSE], c(2, 4, 5), identity)

  .plot_resid_comp(z = x, do_hist = TRUE)
}

#' @importFrom graphics hist box
.plot_resid_comp <- function(x = 1:nrow(z), y = 1:ncol(z), z, xlab = "Year", ylab = "Age", zmax = 2,
                             do_hist = FALSE, ydiff, xdiff) {

  if (all(is.na(z))) return(invisible())
  if (do_hist) return(hist(z, xlab = "Residuals", main = ""))

  if (is.character(y)) {
    yval <- y
    y <- 1:ncol(z)
    yaxt <- "n"
  } else {
    yaxt <- "s"
  }

  zz <- pmin(z, zmax) %>% pmax(-zmax) %>% round(2)
  zlegend <- seq(-zmax, zmax, 0.01)
  cols <- hcl.colors(length(zlegend), palette = "Blue-Red 2", alpha = 1) %>%
    structure(names = zlegend)

  if (missing(ydiff)) {
    ydiff <- diff(y)
    ydiff <- c(ydiff, ydiff[1])
  }

  if (missing(xdiff)) {
    xdiff <- diff(x)
    xdiff <- c(xdiff, xdiff[1])
  }

  border <- ifelse(any(xdiff < 0.5), NA, "grey60")
  rect_diff <- ifelse(any(xdiff < 0.5), 0.475, 0.5)

  plot(NULL, typ = "n", xlab = xlab, ylab = ylab, xlim = range(x), ylim = range(y), yaxt = yaxt)
  sapply(1:nrow(z), function(i) {
    if (any(!is.na(zz[i, ]))) {
      rect(xleft = x[i] - rect_diff * xdiff[i], xright = x[i] + rect_diff * xdiff[i],
           ybottom = y - rect_diff * ydiff, ytop = y + rect_diff * ydiff,
           col = cols[as.character(zz[i, ])],
           border = border)
    }
  })
  legend("topleft", legend = c(-zmax, -zmax/2, 0, zmax/2, zmax),
         col = "grey60", pt.cex = 1, pt.bg = cols[zlegend %in% c(-zmax, -zmax/2, 0, zmax/2, zmax)],
         pch = 22, horiz = TRUE)
  if (yaxt == "n") {
    axis(2, at = y, labels = yval)
  }
  box()

  invisible()
}

