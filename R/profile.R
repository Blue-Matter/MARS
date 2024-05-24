
#' @name profile
#' @aliases profile.MARSassess
#'
#' @title Profile parameters of MARS model
#'
#' @description
#' Evaluate change in objective function and likelihood components for up to 2 parameters.
#'
#' @param fitted [MARSassess-class] object returned by [fit_MARS()]
#' @param p1 Character string that represents the first parameter to be profiled,
#' including the parameter name and index of the vector/array. See "Parameters" section of [make_parameters()].
#' Additionally, this function allows users to specify `R0_s` and `h_s` (in normal units).
#' @param v1 Vector of values corresponding to `p1`
#' @param p2 Character string that represents the optional second parameter to be profiled
#' @param v2 Vector of values corresponding to `p2`
#' @param cores Integer for the number of cores to use for parallel processing (snowfall package)
#' @param ... Not used
#' @return
#' The profile generic returns a data frame of the likelihood values that correspond to
#' fixed values of `p1` and `p2`.
#'
#' - Likelihood `loglike` refers to maximizing the probability of the observed data (higher values for better fit)
#' - Prior `logprior` refers to maximizing the probability of a parameter to their prior distribution (higher values are closer to the prior mode)
#' - Penalty `penalty` are values added to the objective function when parameters exceed model bounds (lower values are better)
#' - `fn` is the objective function returned by RTMB (lower values are better)
#' - `objective` is the objective function returned by the optimizer (lower values are better)
#' @importFrom stats profile
#' @export
profile.MARSassess <- function(fitted, p1, v1, p2, v2, cores = 1, ...) {

  if (cores > 1 && !snowfall::sfIsRunning()) {
    snowfall::sfInit(parallel = TRUE, cpus = cores)
    on.exit(snowfall::sfStop())
  }

  if (snowfall::sfIsRunning()) {
    .lapply <- snowfall::sfLapply
  } else {
    .lapply <- base::lapply
  }

  if (missing(p2)) {
    out <- expand.grid(p1 = v1)
    prof <- .lapply(1:nrow(out), function(i) .prof(fitted, p1, out$p1[i]))
  } else {
    out <- expand.grid(p1 = v1, p2 = v2)
    prof <- .lapply(1:nrow(out), function(i) .prof(fitted, c(p1, p2), c(out$p1[i], out$p2[i])))
  }

  prof_df <- cbind(out, do.call(rbind, prof)) %>%
    structure(class = c("MARSprof", "data.frame"))

  names(prof_df)[1] <- attr(prof_df, "p1") <- p1
  if (!missing(p2)) {
    names(prof_df)[2] <- attr(prof_df, "p2") <- p2
  }

  p <- fitted@obj$env$parList()
  pval1 <- eval(parse(text = paste0("p$", p1)))
  if (is.null(pval1)) pval1 <- eval(parse(text = paste0("fitted@report$", p1)))
  if (is.null(pval1)) pval1 <- NA

  if (missing(p2)) {

    attr(prof_df, "fitted") <- cbind(
      pval1,
      get_likelihood_components(fitted)
    )
    names(attr(prof_df, "fitted"))[1] <- p1

  } else {
    pval2 <- eval(parse(text = paste0("p$", p2)))
    if (is.null(pval2)) pval2 <- eval(parse(text = paste0("fitted@report$", p2)))
    if (is.null(pval2)) pval2 <- NA

    attr(prof_df, "fitted") <- cbind(
      pval1,
      pval2,
      get_likelihood_components(fitted)
    )
    names(attr(prof_df, "fitted"))[1:2] <- c(p1, p2)
  }

  return(prof_df)
}

.prof <- function(fitted, pars, vals) {
  stopifnot(length(pars) == length(vals))

  MARSdata <- get_MARSdata(fitted)
  p <- fitted@obj$env$parList()
  map <- MARSdata@Misc$map
  random <- MARSdata@Misc$random

  for (i in 1:length(pars)) {
    pi <- pars[i]
    vi <- vals[i]

    split_pi <- strsplit(pi, "\\[|\\]")[[1]] # Regex cheat sheet, split by [ or ]
    p_name <- split_pi[1]

    # Transform vi if pi = "R0_s[x]"
    if (grepl("R0_s", pi) && !grepl("t_R0_s", pi)) {
      p_name <- "t_R0_s"
      pi <- sub("R0_s", "t_R0_s", pi)
      vi <- log(vi/MARSdata@Dmodel@scale_s[as.integer(split_pi[2])])
    }

    # Transform vi if pi = "h_s[x]"
    if (grepl("h_s", pi) && !grepl("t_h_s", pi)) {
      p_name <- "t_h_s"
      pi <- sub("h_s", "t_h_s", pi)
      vi <- switch(
        MARSdata@Dstock@SRR_s[as.integer(split_pi[2])],
        "BH" = qlogis((vi - 0.2)/0.8),
        "Ricker" = log(vi - 0.2),
      )
    }

    if (is.null(map[[p_name]])) {
      map[[p_name]] <- if (is.null(dim(p[[p_name]]))) {
        1:length(p[[p_name]])
      } else {
        array(1:length(p[[p_name]]), dim(p[[p_name]]))
      }
    } else if (!is.null(dim(p[[p_name]]))) {
      map[[p_name]] <- array(map[[p_name]], dim(p[[p_name]]))
    }

    map_p <- paste0("map$", pi, " <- NA")
    eval(parse(text = map_p))
    map[[p_name]] <- factor(map[[p_name]])

    assign_p <- paste0("p$", pi, " <- ", vi)
    eval(parse(text = assign_p))
  }

  fit <- fit_MARS(MARSdata, p, map, MARSdata@Misc$random, do_sd = FALSE, silent = TRUE)
  get_likelihood_components(fit)
}

get_likelihood_components <- function(fit) {

  if (length(fit@report)) {
    nm <- names(fit@report)
    nm_like <- grep("loglike", nm, value = TRUE)
    nm_pr <- grep("logprior", nm, value = TRUE)

    nm_out <- c(nm_like, nm_pr, "penalty", "fn")

    out <- lapply(nm_out, function(i) sum(fit@report[[i]])) %>%
      structure(names = nm_out)
  }

  if (length(fit@opt)) {
    out$objective <- fit@opt$objective
  } else {
    out$objective <- NA_real_
  }

  do.call(data.frame, out)
}

#' @rdname profile
#' @aliases plot.MARSprof
#'
#' @param x Output from [profile.MARSassess()]
#' @param component Character for the column in `x` to be plotted
#' @param rel Logical, whether the relative change in `component` is plotted (TRUE) or the raw values (FALSE)
#' @param xlab Optional character for the x-axis label
#' @param ylab Optional character for the y-axis label
#' @param main Optional character for the plot title
#' @param plot2d Character, plotting function for two-dimensional profiling (either a [contour()] or [filled.contour()] plot)
#' @param ... Other argument to the base graphics function, i.e., either plot() or contour()
#' @return
#' The accompanying plot function returns a line plot for a 1-dimensional profile or a contour plot for a two
#' dimensional profile. Will plot the negative log likelihood or negative log prior (better fit with lower values).
#'
#' Relative values are obtained by subtracting from the fitted value. See `attr(x, "fitted")`
#' @importFrom graphics contour filled.contour
#' @importFrom reshape2 acast
#' @export
plot.MARSprof <- function(x, component = "objective", rel = TRUE, xlab, ylab, main,
                          plot2d = c("contour", "filled.contour"), ...) {

  p1 <- attr(x, "p1")
  p2 <- attr(x, "p2")

  if (missing(xlab)) xlab <- p1
  if (is.null(x[[component]])) stop("\"", component, "\" not found in ", substitute(x))

  fitted <- attr(x, "fitted")

  if (is.null(p2)) {
    xplot <- x[[p1]]
    yplot <- x[[component]]
    yfit <- fitted[[component]]
    if (grepl("logprior", component) || grepl("loglike", component)) {
      yplot <- -1 * yplot
      yfit <- -1 * yfit
    }
    if (rel) yplot <- yplot - yfit
    if (missing(ylab)) {
      if (grepl("logprior", component) || grepl("loglike", component)) {
        ylab <- paste("Change in negative", component)
      } else {
        ylab <- paste("Change in", component)
      }
    }
    if (missing(main)) main <- NULL

    plot(xplot, yplot, xlab = xlab, ylab = ylab, typ = "o", main = main, ...)
    abline(v = fitted[[p1]], lty = 2)

  } else {

    plot2d <- match.arg(plot2d)
    plot2d <- match.fun(plot2d)

    names(x)[names(x) == p1] <- "p1"
    names(x)[names(x) == p2] <- "p2"

    zplot <- reshape2::acast(x, list("p1", "p2"), value.var = component)
    zfit <- fitted[[component]]
    if (grepl("logprior", component) || grepl("loglike", component)) {
      zplot <- -1 * zplot
      zfit <- -1 * zfit
    }
    if (rel) zplot <- zplot - zfit
    if (missing(ylab)) ylab <- p2
    if (missing(main)) {
      if (grepl("logprior", component) || grepl("loglike", component)) {
        main <- paste("Change in negative", component)
      } else {
        main <- paste("Change in", component)
      }
    }

    plot2d(
      x = as.numeric(rownames(zplot)), y = as.numeric(colnames(zplot)), z = zplot,
      xlab = xlab, ylab = ylab, main = main,
      ...
    )
    points(fitted[[p1]], fitted[[p2]], col = "red", pch = 16)

  }
  invisible()
}

