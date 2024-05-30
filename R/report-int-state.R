

#' @importFrom grDevices hcl.colors
make_color <- function(n, type = c("fleet", "region", "stock"), alpha = 1) {
  if (n > 1) {
    type <- match.arg(type)
    pal <- switch(
      type,
      "fleet" = "Set2",
      "region" = "Sunset",
      "stock" = "Plasma"
    )
    grDevices::hcl.colors(n, palette = pal, alpha = alpha)
  } else {
    "black"
  }
}

#' @importFrom graphics barplot box legend axis
propplot <- function(x, cols, leg.names, xval, ylab = "Proportion", border = ifelse(nrow(x) > 20, NA, "grey60")) {
  p <- apply(x, 1, function(x) x/sum(x, na.rm = TRUE))
  colnames(p) <- NULL
  if (is.null(dim(p))) p <- matrix(p, length(xval), 1)

  if (missing(cols)) {
    ncat <- nrow(p)
    cols <- make_color(ncat)
  }

  plot(NULL, xlab = "Year", ylab = ylab, xlim = c(0, ncol(p)), xaxt = "n", ylim = c(0, 1), yaxs = "i", xaxs = "i")
  barplot(p, add = TRUE, col = cols, width = 1, space = 0, border = border)
  if (!missing(leg.names) && length(leg.names) > 1) {
    legend("topleft", legend = leg.names, fill = cols, border = border, horiz = TRUE)
  }
  box()

  xt <- pretty(xval)
  xt <- xt[xt %in% xval]
  xp <- match(xt, xval)

  axis(1, at = xp - 0.5, labels = xt)

  invisible()
}


# State variable plots ----
#' @name plot-MARS-state
#' @title Plotting functions for fitted MARS model
#' @description A set of functions to plot state variables (biomass, recruitment time series, etc.)
#' @return Various base graphics plots
NULL

#' @rdname plot-MARS-state
#' @aliases plot_S
#'
#' @param fit [MARSassess-class] object returned by [fit_MARS()]
#' @param by Character to indicate dimension for multivariate plots
#' @param s Integer for the corresponding stock
#' @param prop Logical, whether to plot proportions (TRUE) or absolute numbers
#' @details
#' - `plot_S` plots spawning output by stock or region (either whole numbers or proportions for the latter)
#'
#' @export
plot_S <- function(fit, by = c("total", "stock", "region"), r, s, prop = FALSE) {
  by <- match.arg(by)
  var <- "S_yrs"

  Dlabel <- get_MARSdata(fit)@Dlabel
  year <- Dlabel@year
  ny <- length(year)

  if (by == "total") {
    x <- apply(fit@report[[var]], c(1, 3), sum) # S_ys
    name <- Dlabel@stock

    if (!missing(s)) {
      x <- x[, s, drop = FALSE]
      name <- name[s]
    }

    x <- structure(x, dimnames = list(year = year, stock = name))
  } else if (by == "stock") {

    if (missing(r)) stop("Need region integer r to plot spawning output by stock")

    name <- Dlabel@stock
    ns <- length(name)
    x <- matrix(fit@report[[var]][, r, ], ny, ns) %>% # S_ys
      structure(dimnames = list(year = year, stock = name))

  } else if (by == "region") {
    if (missing(s)) stop("Need stock integer s to plot spawning output by region")
    name <- Dlabel@region
    nr <- length(name)

    x <- matrix(fit@report[[var]][, , s], ny, nr) %>% # S_yr
      structure(dimnames = list(year = year, region = name))
  }
  color <- make_color(ncol(x), type = ifelse(by == "total", "stock", by))

  if (prop) {
    propplot(x, cols = color, leg.names = name, xval = year, ylab = "Spawning fraction")
  } else {

    if (by == "total" && !missing(s)) {
      ylab <- paste(name, "spawning output")
    } else {
      ylab <- "Spawning output"
    }
    matplot(year, x, xlab = "Year", ylab = ylab, typ = "o", col = color, pch = 16,
            ylim = c(0, 1.1) * range(x, na.rm = TRUE), zero_line = TRUE)
    if (ncol(x) > 1) legend("topleft", legend = name, col = color, lwd = 1, pch = 16, horiz = TRUE)
  }

  out <- array2DF(x, responseName = ifelse(prop, "p", "S"))
  invisible(out)
}


#' @rdname plot-MARS-state
#' @aliases plot_B
#' @details
#' - `plot_B` plots total biomass by stock or region (either whole numbers or proportions for the latter)
#'
#' @export
plot_B <- function(fit, by = c("total", "stock", "region"), r, s, prop = FALSE) {
  by <- match.arg(by)
  var <- "B_ymrs"

  Dlabel <- get_MARSdata(fit)@Dlabel
  year <- Dlabel@year
  ny <- length(year)
  nm <- max(length(Dlabel@season), 1)

  if (by == "total") {
    x <- apply(fit@report[[var]], c(1, 2, 4), sum) # B_yms
    name <- Dlabel@stock

    if (!missing(s)) {
      x <- x[, , s, drop = FALSE]
      name <- name[s]
    }
  } else if (by == "stock") {

    if (missing(r)) stop("Need region integer r to plot spawning output by stock")

    name <- Dlabel@stock
    ns <- length(name)
    x <- array(fit@report[[var]][, , r, ], c(ny, nm, ns)) # B_yms

  } else {

    if (missing(s)) stop("Need stock integer s to plot spawning output by region")
    name <- Dlabel@region
    nr <- length(name)

    x <- array(fit@report[[var]][, , , s], c(ny, nm, nr)) # B_ymr

  }

  year <- make_yearseason(year, nm)
  x <- collapse_yearseason(x) # B_tr or B_ts

  color <- make_color(ncol(x), type = ifelse(by == "total", "stock", by))

  if (prop) {
    propplot(x, cols = color, leg.names = name, xval = year, ylab = "Biomass fraction")
  } else {

    if (by == "total" && !missing(s)) {
      ylab <- paste(name, "spawning output")
    } else {
      ylab <- "Spawning output"
    }
    matplot(year, x, xlab = "Year", ylab = "Total biomass", typ = "o", col = color, pch = 16,
            ylim = c(0, 1.1) * range(x, na.rm = TRUE), zero_line = TRUE)
    if (ncol(x) > 1) legend("topleft", legend = name, col = color, lwd = 1, pch = 16, horiz = TRUE)
  }

  out <- array2DF(x, responseName = ifelse(prop, "p", "B"))
  invisible(out)
}



#' @rdname plot-MARS-state
#' @aliases plot_R
#' @details
#' - `plot_R` plots recruitment by stock
#'
#' @export
plot_R <- function(fit, s) {
  var <- "R_ys"

  Dlabel <- get_MARSdata(fit)@Dlabel
  year <- Dlabel@year

  if (missing(s)) {
    x <- fit@report[[var]]
    name <- Dlabel@stock
    ylab <- "Recruitment"
  } else {
    x <- fit@report[[var]][, s, drop = FALSE]
    name <- Dlabel@stock[s]
    ylab <- paste(name, "recruitment")
  }

  color <- make_color(ncol(x), "stock")
  matplot(year, x, xlab = "Year", ylab = ylab, typ = "o", col = color, pch = 16,
          ylim = c(0, 1.1) * range(x, na.rm = TRUE), zero_line = TRUE)
  if (ncol(x) > 1) legend("topleft", legend = name, col = color, lwd = 1, pch = 16, horiz = TRUE)

  x <- structure(x, dimnames = list(year = year, stock = name))
  invisible(array2DF(x, "R"))
}


#' @rdname plot-MARS-state
#' @aliases plot_SRR
#' @param phi Logical, whether to plot unfished replacement line
#' @details
#' - `plot_SRR` plots the stock-recruitment relationship and history (realized recruitment) by stock
#' @importFrom graphics points
#' @export
plot_SRR <- function(fit, s = 1, phi = TRUE) {
  Dlabel <- get_MARSdata(fit)@Dlabel

  S_y <- apply(fit@report$S_yrs[, , s, drop = FALSE], 1, sum)
  R_y <- fit@report$R_ys[, s]
  Rpred_y <- R_y/fit@report$Rdev_ys[, s]

  plot(S_y[order(S_y)], Rpred_y[order(S_y)], typ = "l", lwd = 2, xlab = "Spawning output", ylab = "Recruitment",
       xaxs = "i", yaxs = "i", xlim = c(0, 1.1) * range(S_y), ylim = c(0, 1.1) * range(R_y))
  points(S_y, R_y, pch = 16)

  if (phi) {
    phi_s <- fit@report$phi_s[s]
    abline(a = 0, b = 1/phi_s, lty = 2)
  }
  invisible()
}


#' @rdname plot-MARS-state
#' @aliases plot_Rdev
#' @param log Logical, whether to plot the natural logarithm of the response variable
#' @details
#' - `plot_Rdev` plots recruitment deviations by stock
#' @importFrom graphics arrows
#' @export
plot_Rdev <- function(fit, s = 1, log = TRUE) {

  Dlabel <- get_MARSdata(fit)@Dlabel
  year <- Dlabel@year

  if (log) {
    if (length(fit@SD)) {
      x <- as.list(fit@SD, what = "Estimate")$log_rdev_ys[, s]
      std <- as.list(fit@SD, what = "Std. Error")$log_rdev_ys[, s]
      std[is.na(std)] <- 0
    } else {
      x <- log(fit@report$Rdev_ys[, s])
      std <- numeric(length(x))
    }

    upper <- x + 1.96 * std
    lower <- x - 1.96 * std

    plot(year, x, xlab = "Year", ylab = "log Recruitment deviations", typ = "o", pch = 16,
         ylim = range(lower, upper), lty = 3)
    arrows(x0 = year, y0 = lower, y1 = upper, length = 0)
    abline(h = 0, lty = 2)

  } else {
    x <- fit@report$Rdev_ys[, s]
    plot(year, x, xlab = "Year", ylab = "Recruitment deviations", typ = "o", pch = 16,
         ylim = c(0, 1.1) * range(x, na.rm = TRUE), zero_line = TRUE)

    abline(h = 1, lty = 2)
  }
  invisible()
}

#' @rdname plot-MARS-state
#' @aliases plot_Fstock
#' @details
#' - `plot_Fstock` plots apical instantaneous fishing mortality (per year) by stock
#'
#' @export
plot_Fstock <- function(fit, s) {
  var <- "F_yas"

  Dlabel <- get_MARSdata(fit)@Dlabel
  year <- Dlabel@year

  if (missing(s)) {
    x <- fit@report[[var]]
    name <- Dlabel@stock
    ylab <- "Apical fishing mortality"
  } else {
    x <- fit@report[[var]][, , s, drop = FALSE]
    name <- Dlabel@stock[s]
    ylab <- paste(name, "apical fishing mortality")
  }

  x <- apply(x, c(1, 3), max)
  x[is.infinite(x)] <- NA

  color <- make_color(ncol(x), "stock")
  matplot(year, x, xlab = "Year", ylab = ylab, typ = "o", col = color, pch = 16,
          ylim = c(0, 1.1) * range(x, na.rm = TRUE), zero_line = TRUE)
  if (ncol(x) > 1) legend("topleft", legend = name, col = color, lwd = 1, pch = 16, horiz = TRUE)

  invisible()
}



#' @rdname plot-MARS-state
#' @aliases plot_self
#' @param f Integer for the corresponding fleet
#' @details
#' - `plot_self` plots fishery selectivity
#' @export
plot_self <- function(fit, f = 1) {
  dat <- get_MARSdata(fit)
  Dfishery <- dat@Dfishery

  sel_block <- Dfishery@sel_block_yf[, f]
  sel_b <- Dfishery@sel_f[unique(sel_block)]
  fname <- dat@Dlabel@fleet[f]

  year <- dat@Dlabel@year

  if (all(grepl("length", sel_b))) {
    lmid <- dat@Dmodel@lmid
    x <- fit@report$sel_lf[, unique(sel_block), drop = FALSE]

    color <- make_color(ncol(x), "fleet")
    matplot(lmid, x, xlab = "Length", ylab = paste(fname, "selectivity"),
            typ = "o", col = color, pch = 16,
            ylim = c(0, 1), lty = 1, zero_line = TRUE)
    if (ncol(x) > 1) {
      name <- sapply(unique(sel_block), function(i) {
        y <- year[sel_block == i]
        if (length(y) == 1) {
          return(y)
        } else {
          return(paste(range(y), collapse = "-"))
        }
      })
      legend("topright", legend = name, col = color, lwd = 1, pch = 16)
    }
  } else {

    m <- 1
    s <- 1

    x <- fit@report$sel_ymafs[, m, , f, s]
    xx <- apply(x, 2, diff)

    if (any(xx != 0)) {
      ybreak <- c(1, which(rowSums(xx) > 0) + 1)
      name <- sapply(1:length(ybreak), function(i) {
        if (i == length(ybreak)) {
          y <- c(year[ybreak[i]], year[length(year)])
        } else {
          y <- year[c(ybreak[i], ybreak[i+1] - 1)]
        }
        paste(range(y), collapse = "-")
      })
      x <- x[ybreak, , drop = FALSE]

    } else {
      x <- x[1, , drop = FALSE]
    }
    age <- dat@Dlabel@age

    color <- make_color(nrow(x), "fleet")
    matplot(age, t(x), xlab = "Age", ylab = paste(fname, "selectivity"),
            typ = "o", col = color, pch = 16,
            ylim = c(0, 1), lty = 1, zero_line = TRUE)
    if (nrow(x) > 1) legend("topright", legend = name, col = color, lwd = 1, pch = 16)
  }
  invisible()
}

#' @rdname plot-MARS-state
#' @aliases plot_seli
#' @param i Integer for the corresponding survey
#' @details
#' - `plot_seli` plots index selectivity
#' @export
plot_seli <- function(fit, i = 1) {
  dat <- get_MARSdata(fit)
  sel_i <- dat@Dsurvey@sel_i[i]
  mirror_f <- suppressWarnings(as.numeric(sel_i))

  iname <- dat@Dlabel@index[i]

  if (!is.na(mirror_f)) {
    plot_self(fit, f = mirror_f)
  } else if (grepl("length", sel_i)) {
    lmid <- dat@Dmodel@lmid
    x <- fit@report$sel_li[, i]

    plot(lmid, x, xlab = "Length", ylab = paste(iname, "selectivity"),
         typ = "o", pch = 16,
         ylim = c(0, 1), lty = 1, zero_line = TRUE)
  } else {

    m <- 1
    s <- 1

    x <- fit@report$sel_ymais[, m, , i, s]
    xx <- apply(x, 2, diff)

    if (any(xx != 0)) {
      year <- dat@Dlabel@year
      ybreak <- c(1, which(rowSums(xx) > 0) + 1)
      name <- sapply(1:length(ybreak), function(i) {
        if (i == length(ybreak)) {
          y <- c(year[ybreak[i]], year[length(year)])
        } else {
          y <- year[c(ybreak[i], ybreak[i+1] - 1)]
        }
        paste(range(y), collapse = "-")
      })
      x <- x[ybreak, , drop = FALSE]

    } else {
      x <- x[1, , drop = FALSE]
    }
    age <- dat@Dlabel@age
    color <- make_color(nrow(x), "fleet")

    matplot(age, t(x), xlab = "Age", ylab = paste(iname, "selectivity"),
            typ = "o", col = color, pch = 16,
            ylim = c(0, 1), lty = 1, zero_line = TRUE)
    if (nrow(x) > 1) legend("topright", legend = name, col = color, lwd = 1, pch = 16)
  }
  invisible()
}


#' @rdname plot-MARS-state
#' @aliases plot_selstock
#' @param plot2d Character, plotting function for either a [contour()] or [filled.contour()] plot
#' @param ... Other argument to the base graphics function
#' @details
#' - `plot_selstock` plots the realized selectivity from total catch and total abundance at age
#' @export
#' @importFrom graphics contour filled.contour
plot_selstock <- function(fit, s = 1, plot2d = c("contour", "filled.contour"), ...) {

  plot2d <- match.arg(plot2d)
  plot2d <- match.fun(plot2d)

  dat <- get_MARSdata(fit)
  year <- dat@Dlabel@year
  age <- dat@Dlabel@age

  sel_ya <- fit@report$F_yas[, , s] %>%
    apply(1, function(x) x/max(x)) %>%
    t()

  plot2d(x = year, y = age, xlab = "Year", ylab = "Age", z = sel_ya, levels = seq(0, 1, 0.1), ...)

  invisible()
}

#' @rdname plot-MARS-state
#' @aliases plot_N
#' @param m Integer for the corresponding season
#' @param r Integer for the corresponding region
#' @param ... Other argument to the base graphics function
#' @details
#' - `plot_N` reports total abundance at age
#' @export
#' @importFrom graphics contour filled.contour
plot_N <- function(fit, m = 1, r, s = 1, plot2d = c("contour", "filled.contour"), ...) {
  plot2d <- match.arg(plot2d)
  plot2d <- match.fun(plot2d)

  dat <- get_MARSdata(fit)
  if (missing(r)) r <- 1:dat@Dmodel@nr
  if (length(m) > 1) stop("length(m) should be one")
  ny <- dat@Dmodel@ny
  year <- dat@Dlabel@year
  age <- dat@Dlabel@age

  N_ya <- apply(fit@report$N_ymars[1:ny, m, , r, s, drop = FALSE], c(1, 3), sum)

  dots <- list(...)
  if (!length(dots$nlevels)) dots$nlevels <- 10
  if (!length(dots$levels)) dots$levels <- exp(pretty(log(range(N_ya)), 10))

  dots$x <- year
  dots$y <- age
  dots$xlab = "Year"
  dots$ylab <- "Age"
  dots$z <- N_ya

  do.call(plot2d, dots)

  invisible()
}


#' @rdname plot-MARS-state
#' @aliases plot_V
#' @details
#' - `plot_V` plots vulnerable biomass, availability to the fishery
#' @export
plot_V <- function(fit, f = 1, by = c("stock", "region"), prop = FALSE) {
  by <- match.arg(by)

  var <- "VB_ymfrs"

  Dlabel <- get_MARSdata(fit)@Dlabel
  year <- Dlabel@year
  nm <- max(length(Dlabel@season), 1)

  if (by == "stock") {
    x <- apply(fit@report[[var]][, , f, , , drop = FALSE], c(1, 2, 5), sum)
    name <- Dlabel@stock
  } else {
    x <- apply(fit@report[[var]][, , f, , , drop = FALSE], c(1, 2, 4), sum)
    name <- Dlabel@region
  }

  year <- make_yearseason(year, nm)
  x <- collapse_yearseason(x)

  color <- make_color(ncol(x), type = by)

  fname <- Dlabel@fleet[f]
  ylab <- paste("Biomass available to", fname)

  if (prop) {
    propplot(x, cols = color, leg.names = name, xval = year, ylab = ylab)
  } else {
    matplot(year, x, xlab = "Year", ylab = ylab, typ = "o", col = color, pch = 16,
            ylim = c(0, 1.1) * range(x, na.rm = TRUE), zero_line = TRUE)
    if (ncol(x) > 1) legend("topleft", legend = name, col = color, lwd = 1, pch = 16, horiz = TRUE)
  }

  invisible()
}


#' @rdname plot-MARS-state
#' @aliases plot_Ffleet
#' @details
#' - `plot_Ffleet` plots apical instantaneous fishing mortality (per season) by fleet
#'
#' @export
plot_Ffleet <- function(fit, f = 1) {
  var <- "F_ymfr"

  Dlabel <- get_MARSdata(fit)@Dlabel
  year <- Dlabel@year
  nm <- max(length(Dlabel@season), 1)

  x <- apply(fit@report[[var]][, , f, , drop = FALSE], c(1, 2, 4), identity)
  name <- Dlabel@region

  year <- make_yearseason(year, nm)
  x <- collapse_yearseason(x)

  color <- make_color(ncol(x), type = "region")

  fname <- Dlabel@fleet[f]
  ylab <- paste(fname, "fishing mortality")
  matplot(year, x, xlab = "Year", ylab = ylab, typ = "o", col = color, pch = 16,
          ylim = c(0, 1.1) * range(x, na.rm = TRUE), zero_line = TRUE)
  if (ncol(x) > 1) legend("topleft", legend = name, col = color, lwd = 1, pch = 16, horiz = TRUE)

  invisible()
}

#' @rdname plot-MARS-state
#' @aliases plot_mov
#' @param y Integer, year for plotting the movement matrix
#' @param a Integer, corresponding age for plotting the movement matrix
#' @details
#' - `plot_mov` plots movement matrices and the corresponding equilibrium distribution in multi-area models
#' @export
#' @importFrom graphics title
plot_mov <- function(fit, s = 1, y, a) {

  dat <- get_MARSdata(fit)

  nm <- dat@Dmodel@nm
  nr <- dat@Dmodel@nr
  if (missing(y)) y <- dat@Dmodel@y_phi
  if (missing(a)) a <- 2
  rname <- dat@Dlabel@region
  mname <- dat@Dlabel@season

  mov <- array(fit@report$mov_ymarrs[y, , a, , , s], c(nm, nr, nr))
  dist_eq <- mov_proj(mov, start = fit@report$recdist_rs[, s])

  #dist_start <- ifelse(dat@Dstock@presence_rs[, s], 1/sum(dat@Dstock@presence_rs[, s]), 0)
  #dist_eq <- mov_proj(mov, start = dist_start)

  if (nm > 1) {
    old_mar <- par()$mar
    old_mfrow = par()$mfrow
    par(mar = c(5, 4, 1, 1))
    on.exit(par(mar = old_mar, mfrow = old_mfrow))

    par(mfrow = c(2, ceiling(nm/2)))
  }

  for(m in 1:nm) {
    .plot_mov(m = mov[m, , ], p = dist_eq[m, ], rname = rname, xlab = "", ylab = "")
    if (nm > 1) title(mname[m], font.main = 1)
  }

  par(mfrow = c(1, 1))
  mtext("Destination", side = 1, line = 4)
  mtext("Origin", side = 2, line = 3)

  invisible()
}

#' @rdname plot-MARS-state
#' @aliases plot_recdist
#' @details
#' - `plot_recdist` plots the distribution of recruitment for each stock
#' @export
plot_recdist <- function(fit) {
  dat <- get_MARSdata(fit)

  nr <- dat@Dmodel@nr

  if (nr > 1) {
    ns <- dat@Dmodel@ns
    rname <- dat@Dlabel@region
    sname <- dat@Dlabel@stock

    recdist <- fit@report$recdist

    vcol <- rainbow(100, end = 0.45)

    graphics::plot.default(
      NULL, xlab = "Stock", ylab = "Region", xaxs = "i", yaxs = "i",
      xaxt = "n", yaxt = "n", xlim = c(1, ns+1), ylim = c(1, nr+1)
    )
    for(x in 1:ns) {
      for(y in 1:nr) {
        m_yx <- round(recdist[y, x], 2)
        rect(xleft = x, ybottom = y, xright = x+1, ytop = y+1, col = vcol[100 * m_yx])
        text(x + 0.5, y + 0.5, m_yx)
      }
    }

    axis(1, at = 1:ns + 0.5, labels = as.character(sname), font = 2, cex.axis = 0.75)
    axis(2, at = 1:nr + 0.5, labels = as.character(rname), font = 2, cex.axis = 0.75)
  }

  invisible()
}

#' @importFrom grDevices rainbow
#' @importFrom graphics rect text
.plot_mov <- function(m, p, xlab = "Destination", ylab = "Origin",
                      nr = length(p), rname = paste("Region", 1:nr)) {

  vcol <- rainbow(100, end = 0.45)

  graphics::plot.default(
    NULL, xlab = xlab, ylab = ylab, xaxs = "i", yaxs = "i",
    xaxt = "n", yaxt = "n", xlim = c(1, nr+3), ylim = c(1, nr+1)
  )
  for(x in 1:nr) {
    for(y in 1:nr) {
      m_yx <- round(m[y, x], 2)
      rect(xleft = x, ybottom = y, xright = x+1, ytop = y+1, col = vcol[100 * m_yx])
      text(x + 0.5, y + 0.5, m_yx)
    }
  }
  eq_val <- round(p, 2)
  eq_col <- ifelse(eq_val > 0, vcol[100 * eq_val], NA)
  rect(nr + 2, ybottom = 1:nr, xright = nr + 3, ytop = 1:nr + 1, col = eq_col)
  text(nr + 2.5, 1:nr + 0.5, eq_val)

  axis(2, at = 1:nr + 0.5, labels = as.character(rname), font = 2, cex.axis = 0.75)
  axis(1, at = c(1:nr, nr+2) + 0.5, labels = c(as.character(rname), "Eq."), font = 2, cex.axis = 0.75)

  invisible()
}


mov_proj <- function(x, nm = dim(x)[1], nr = dim(x)[2], start = rep(1/nr, nr), nit = 20) {
  N <- array(NA, c(nit, nm, nr))
  N[1, 1, ] <- start
  for (i in 2:nit - 1) {
    for(m in 1:nm) {
      if (m == nm) {
        N[i+1, 1, ] <- N[i, m, ] %*% x[m, , ]
      } else {
        N[i, m+1, ] <- N[i, m, ] %*% x[m, , ]
      }
    }
  }
  return(array(N[nit-1, , ], c(nm, nr)))
}
