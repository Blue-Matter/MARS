
# Data plots ----
#' @name plot-MARS-data
#' @title Plotting functions for data in MARS model
#' @description A set of functions to plot data variables and predicted values (catch, age composition, etc.)
#' @return Various base graphics plots
NULL

#' @rdname plot-MARS-data
#' @aliases plot_catch
#'
#' @param fit [MARSassess-class] object returned by [fit_MARS()]
#' @param by Character to indicate dimension for multivariate plots
#' @param f Integer for the corresponding fleet
#' @param prop Logical, whether to plot proportions (TRUE) or absolute numbers
#' @param annual Logical, whether to plot annual values (summed over seasons)
#' @details
#' - `plot_catch` plots the fishery catch by stock or region (either whole numbers or proportions)
#'
#' @export
plot_catch <- function(fit, f = 1, by = c("region", "stock"), prop = FALSE, annual = FALSE) {
  by <- match.arg(by)
  var <- "CB_ymfrs"

  Dlabel <- get_MARSdata(fit)@Dlabel
  year <- Dlabel@year
  nm <- max(length(Dlabel@season), 1)

  if (by == "stock") {
    x <- apply(fit@report[[var]][, , f, , , drop = FALSE], c(1, 2, 5), sum) # C_yms
    name <- Dlabel@stock
  } else {
    x <- apply(fit@report[[var]][, , f, , , drop = FALSE], c(1, 2, 4), sum)  # C_ymr
    name <- Dlabel@region
  }

  if (annual) {
    x <- apply(x, c(1, 3), sum)
  } else {
    year <- make_yearseason(year, nm)
    x <- collapse_yearseason(x)
  }

  color <- make_color(ncol(x), type = by)
  fname <- Dlabel@fleet[f]

  if (prop) {
    propplot(x, cols = color, leg.names = name, xval = year, ylab = paste(fname, "catch proportion"))
  } else {
    matplot(year, x, xlab = "Year", ylab = paste(fname, "catch"), typ = "o", col = color, pch = 16,
            ylim = c(0, 1.1) * range(x, na.rm = TRUE), zero_line = TRUE)
    if (ncol(x) > 1) legend("topleft", legend = name, col = color, lwd = 1, pch = 16, horiz = TRUE)
  }

  invisible()
}

#' @rdname plot-MARS-data
#' @aliases plot_index
#' @param i Integer, indexes the survey
#' @param zoom Logical, for `plot_index()`. If \code{TRUE}, plots a subset of years with observed data points. Otherwise, plots
#' predicted values over all model years.
#' @details
#' - `plot_index` plots indices of abundance
#' @export
plot_index <- function(fit, i = 1, zoom = FALSE) {
  dat <- get_MARSdata(fit)
  Dlabel <- dat@Dlabel
  year <- Dlabel@year
  nm <- max(length(Dlabel@season), 1)

  iname <- Dlabel@index[i]

  ipred <- apply(fit@report$I_ymi[, , i, drop = FALSE], 1:2, identity)
  iobs <- apply(dat@Dsurvey@Iobs_ymi[, , i, drop = FALSE], 1:2, identity)
  isd <- apply(dat@Dsurvey@Isd_ymi[, , i, drop = FALSE], 1:2, identity)

  year <- make_yearseason(year, nm)
  ipred <- collapse_yearseason(ipred)
  iobs <- collapse_yearseason(iobs)
  isd <- collapse_yearseason(isd)

  ind <- !is.na(iobs)

  if (sum(ind)) {
    if (zoom) {
      year <- year[ind]
      ipred <- ipred[ind]
      iobs <- iobs[ind]
      isd <- isd[ind]
    }

    iupper <- exp(log(iobs) + 1.96 * isd)
    ilower <- exp(log(iobs) - 1.96 * isd)

    plot(year, iobs, xlab = "Year", ylab = iname, typ = "p", pch = 16,
         ylim = c(0, 1.1) * range(ipred, iupper, na.rm = TRUE), zero_line = TRUE)
    arrows(year, y0 = ilower, y1 = iupper, length = 0)
    lines(year, ipred, lwd = 2, col = 2, typ = ifelse(length(year) > 10, "l", "o"))
  }

  invisible()
}


#' @rdname plot-MARS-data
#' @aliases plot_CAA
#' @param f Integer, indexes the fleet
#' @param r Integer, indexes the region
#' @details
#' - `plot_CAA` plots the fishery catch at age
#' @export
plot_CAA <- function(fit, f = 1, r = 1) {
  dat <- get_MARSdata(fit)

  if (sum(dat@Dfishery@CAAN_ymfr, na.rm = TRUE)) {
    N <- apply(dat@Dfishery@CAAN_ymfr[, , f, r, drop = FALSE], 1:2, identity) %>% t() %>% as.numeric()
    N[is.na(N)] <- 0

    if (any(N > 0)) {
      Dlabel <- dat@Dlabel
      nm <- max(length(Dlabel@season), 1)

      year <- Dlabel@year
      year <- make_yearseason(year, nm)

      pred <- apply(fit@report$CN_ymafrs[, , , f, r, , drop = FALSE], 1:3, sum)
      pred <- collapse_yearseason(pred) %>% apply(1, function(x) x/sum(x, na.rm = TRUE)) %>% t()
      pred[is.na(pred)] <- 0

      obs <- apply(dat@Dfishery@CAAobs_ymafr[, , , f, r, drop = FALSE], 1:3, identity)
      obs <- collapse_yearseason(obs) %>% apply(1, function(x) x/sum(x, na.rm = TRUE)) %>% t()

      include <- rowSums(obs, na.rm = TRUE) > 0

      plot_composition(obs[include, , drop = FALSE], pred[include, , drop = FALSE], xval = Dlabel@age, ylab = "Proportion",
                       zval = year[include], N = N[include])
    }
  }

  invisible()
}

#' @rdname plot-MARS-data
#' @aliases plot_CAL
#' @details
#' - `plot_CAL` plots the catch at length
#' @export
plot_CAL <- function(fit, f = 1, r = 1) {
  dat <- get_MARSdata(fit)

  if (sum(dat@Dfishery@CALN_ymfr, na.rm = TRUE)) {
    N <- apply(dat@Dfishery@CALN_ymfr[, , f, r, drop = FALSE], 1:2, identity) %>% t() %>% as.numeric()
    N[is.na(N)] <- 0

    if (any(N > 0)) {
      Dlabel <- dat@Dlabel
      nm <- max(length(Dlabel@season), 1)

      year <- Dlabel@year
      year <- make_yearseason(year, nm)

      pred <- apply(fit@report$CN_ymlfrs[, , , f, r, , drop = FALSE], 1:3, sum)
      pred <- collapse_yearseason(pred) %>% apply(1, function(x) x/sum(x, na.rm = TRUE)) %>% t()
      pred[is.na(pred)] <- 0

      obs <- apply(dat@Dfishery@CALobs_ymlfr[, , , f, r, drop = FALSE], 1:3, identity)
      obs <- collapse_yearseason(obs) %>% apply(1, function(x) x/sum(x, na.rm = TRUE)) %>% t()

      include <- rowSums(obs, na.rm = TRUE) > 0

      plot_composition(obs[include, , drop = FALSE], pred[include, , drop = FALSE], xval = dat@Dmodel@lmid,
                       xlab = "Length", ylab = "Proportion",
                       zval = year[include], N = N[include])
    }
  }

  invisible()
}

#' @rdname plot-MARS-data
#' @aliases plot_IAA
#' @details
#' - `plot_IAA` plots the index age composition
#' @export
plot_IAA <- function(fit, i = 1) {
  dat <- get_MARSdata(fit)

  if (sum(dat@Dsurvey@IAAN_ymi, na.rm = TRUE)) {
    N <- apply(dat@Dsurvey@IAAN_ymi[, , i, drop = FALSE], 1:2, identity) %>% t() %>% as.numeric()
    N[is.na(N)] <- 0

    if (any(N > 0)) {
      Dlabel <- dat@Dlabel
      nm <- max(length(Dlabel@season), 1)

      year <- Dlabel@year
      year <- make_yearseason(year, nm)

      pred <- apply(fit@report$IN_ymais[, , , i, , drop = FALSE], 1:3, sum)
      pred <- collapse_yearseason(pred) %>% apply(1, function(x) x/sum(x, na.rm = TRUE)) %>% t()
      pred[is.na(pred)] <- 0

      obs <- apply(dat@Dsurvey@IAAobs_ymai[, , , i, drop = FALSE], 1:3, identity)
      obs <- collapse_yearseason(obs) %>% apply(1, function(x) x/sum(x, na.rm = TRUE)) %>% t()

      include <- rowSums(obs, na.rm = TRUE) > 0

      plot_composition(obs[include, , drop = FALSE], pred[include, , drop = FALSE], xval = Dlabel@age, ylab = "Proportion",
                       zval = year[include], N = N[include])
    }
  }

  invisible()
}

#' @rdname plot-MARS-data
#' @aliases plot_IAL
#' @details
#' - `plot_IAL` plots the index length composition
#' @export
plot_IAL <- function(fit, i = 1) {
  dat <- get_MARSdata(fit)

  if (sum(dat@Dsurvey@IALN_ymi, na.rm = TRUE)) {
    N <- apply(dat@Dsurvey@IALN_ymi[, , i, drop = FALSE], 1:2, identity) %>% t() %>% as.numeric()
    N[is.na(N)] <- 0

    if (any(N > 0)) {
      Dlabel <- dat@Dlabel
      nm <- max(length(Dlabel@season), 1)

      year <- Dlabel@year
      year <- make_yearseason(year, nm)

      pred <- apply(fit@report$IN_ymlis[, , , i, , drop = FALSE], 1:3, sum)
      pred <- collapse_yearseason(pred) %>% apply(1, function(x) x/sum(x, na.rm = TRUE)) %>% t()
      pred[is.na(pred)] <- 0

      obs <- apply(dat@Dsurvey@IALobs_ymli[, , , i, drop = FALSE], 1:3, identity)
      obs <- collapse_yearseason(obs) %>% apply(1, function(x) x/sum(x, na.rm = TRUE)) %>% t()

      include <- rowSums(obs, na.rm = TRUE) > 0

      plot_composition(obs[include > 0, ], pred[include > 0, ], xval = dat@Dmodel@lmid,
                       xlab = "Length", ylab = "Proportion",
                       zval = year[include > 0], N = N[include > 0])
    }
  }

  invisible()
}

#' @rdname plot-MARS-data
#' @param ff Integer, indexes the aggregate fleet (for stock composition data)
#' @param aa Integer, indexes the aggregate age class (for stock composition and tag data)
#' @aliases plot_SC
#' @details
#' - `plot_SC` plots the stock composition
#' @importFrom graphics matlines
#' @export
plot_SC <- function(fit, ff = 1, aa = 1, r = 1, prop = FALSE) {
  dat <- get_MARSdata(fit)

  if (dat@Dmodel@ns == 1) stop("Stock composition figure not needed.")

  if (sum(dat@Dfishery@SC_ymafrs, na.rm = TRUE)) {
    N <- apply(dat@Dfishery@SC_ymafrs[, , aa, ff, r, , drop = FALSE], 1:2, sum)
    N[is.na(N)] <- 0

    if (any(N > 0)) {
      Dlabel <- dat@Dlabel
      nm <- max(length(Dlabel@season), 1)

      year <- Dlabel@year
      year <- make_yearseason(year, nm)

      pred <- apply(fit@report$SCpred_ymafrs[, , aa, ff, r, , drop = FALSE], c(1, 2, 6), identity)
      pred <- collapse_yearseason(pred) %>% apply(1, function(x) x/sum(x, na.rm = TRUE)) %>% t()
      pred[is.na(pred)] <- 0

      N <- collapse_yearseason(N)

      obs <- apply(dat@Dfishery@SC_ymafrs[, , aa, ff, r, , drop = FALSE], c(1, 2, 6), identity)
      obs <- collapse_yearseason(obs) %>% apply(1, function(x) x/sum(x, na.rm = TRUE)) %>% t()

      if (prop) {
        color <- make_color(ncol(pred), type = "stock")
        propplot(pred, cols = color, leg.names = Dlabel@stock, xval = year)

        if (dat@Dmodel@ns == 2) {
          obs_cumsum <- apply(obs, 1, cumsum) %>% t()
          matlines(obs_cumsum[, -dat@Dmodel@ns, drop = FALSE], type = "o", col = 2, pch = 16, ylim = c(0, 1))
        }

      } else {
        include <- rowSums(obs, na.rm = TRUE) > 0

        plot_composition(obs[include, , drop = FALSE], pred[include, , drop = FALSE], xval = 1:dat@Dmodel@ns,
                         xlab = "Stock", ylab = "Proportion",
                         zval = year[include], N = N[include],
                         xaxislab = Dlabel@stock)
      }
    }
  }

  invisible()
}

#' @importFrom graphics lines mtext
plot_composition <- function(obs, pred = NULL, xval = 1:ncol(obs), xlab = "Age",
                             ylab = "Value", zval = 1:nrow(obs), N = rowSums(obs),
                             xaxislab = 1:ncol(obs), ncol = 4, nrow = 4) {

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(mfcol = c(nrow, ncol), mar = rep(0, 4), oma = c(5.1, 5.1, 2.1, 2.1))

  ylim <- c(0, 1.1) * range(pred, obs, na.rm = TRUE)

  yaxp <- c(0, max(pretty(ylim, n = 4)), 4)
  if (max(obs, pred, na.rm = TRUE) == 1) yaxp <- c(0, 1, 4)
  las <- 1
  nplot <- nrow(pred)
  nset <- nrow * ncol
  xaxt_manual <- any(xaxislab != xval)
  N <- round(N, 2)

  for(i in 1:nplot) {
    yaxt <- ifelse(i %% nset %in% c(1:nrow), "s", "n") # TRUE = first column
    xaxt <- ifelse(i < nplot && i %% nrow %in% c(1:(nrow-1)), "n", "s") # TRUE = first three rows

    plot(xval, obs[i, ], typ = "o", ylim = ylim, yaxp = yaxp, xaxt = ifelse(xaxt_manual, "n", xaxt), yaxt = yaxt, las = las)
    if (xaxt == "s") axis(1, at = xval, labels = xaxislab)
    if (!is.null(pred)) lines(xval, pred[i, ], lwd = 2, col = 2)
    legend("topright", legend = c(zval[i], ifelse(is.null(N), "", paste0("N = ", N[i]))), bty = "n", xjust = 1)

    if (i %% nset == 1) {
      mtext(xlab, side = 1, line = 3, outer = TRUE)
      mtext(ylab, side = 2, line = 3.5, outer = TRUE)
    }
  }

  invisible()
}

#' @rdname plot-MARS-data
#' @param yy Integer, indexes the aggregate years (for the tag data)
#' @param s Integer, indexes the stock
#' @aliases plot_tagmov
#' @details
#' - `plot_tagmov` plots the tag movements
#' @export
plot_tagmov <- function(fit, s = 1, yy = 1, aa = 1) {
  dat <- get_MARSdata(fit)

  if (dat@Dmodel@nr == 1) stop("Tag movement figure not needed.")

  if (sum(dat@Dtag@tag_ymarrs, na.rm = TRUE)) {
    N <- apply(dat@Dtag@tag_ymarrs[yy, , aa, , , s, drop = FALSE], c(2, 4), sum)
    N[is.na(N)] <- 0

    if (any(N > 0)) {
      Dlabel <- dat@Dlabel

      pred <- apply(fit@report$tagpred_ymarrs[yy, , aa, , , s, drop = FALSE], c(2, 4, 5), identity) %>%
        collapse_yearseason()

      obs <- apply(dat@Dtag@tag_ymarrs[yy, , aa, , , s, drop = FALSE], c(2, 4, 5), identity)
      obs2 <- apply(obs, c(1, 2), function(x) x/sum(x, na.rm = TRUE)) %>% aperm(c(2, 3, 1)) %>%
        collapse_yearseason() # Season + region of origin
      obs2[is.na(obs2)] <- 0

      N <- collapse_yearseason(N)

      nrow <- dat@Dmodel@nr
      ncol <- dat@Dmodel@nm

      z <- paste(
        paste("Origin:", rep(Dlabel@region, dat@Dmodel@nm)),
        paste(rep(paste("\nSeason", 1:dat@Dmodel@nm), each = dat@Dmodel@nr))
      )

      plot_composition(obs2, pred, xval = 1:dat@Dmodel@nr,
                       xlab = "Destination", ylab = "Proportion",
                       zval = z, N = N,
                       xaxislab = Dlabel@region, ncol = ncol, nrow = nrow)
    }
  }

  invisible()
}
