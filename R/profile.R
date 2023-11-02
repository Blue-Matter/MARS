
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
#' including the parameter name and index. See "Parameters" section of [make_parameters()]
#' @param v1 Vector of values corresponding to `p1`
#' @param p2 Character string that represents the optional second parameter to be profiled
#' @param v2 Vector of values corresponding to `p2`
#' @param ... Not used
#' @return
#' The profile generic returns a data frame of the likelihood values that correspond to
#' fixed values of `p1` and `p2`.
#'
#' @importFrom stats profile
#' @export
profile.MARSassess <- function(fitted, p1, v1, p2, v2, ...) {

  if (missing(p2)) {
    out <- expand.grid(p1 = v1)
    prof <- lapply(1:nrow(out), function(i) .prof(fitted, p1, out$p1[i]))
  } else {
    out <- expand.grid(p1 = v1, p2 = v2)
    prof <- lapply(1:nrow(out), function(i) .prof(fitted, c(p1, p2), c(out$p1[i], out$p2[i])))
  }

  prof_df <- cbind(out, do.call(rbind, prof)) %>%
    structure(class = c("MARSprof", "data.frame"))

  names(prof_df)[1] <- attr(prof_df, "p1") <- p1
  if (!missing(p2)) {
    names(prof_df)[2] <- attr(prof_df, "p2") <- p2
  }
  attr(prof_df, "fitted") <- get_likelihood_components(fitted)

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

    p_name <- strsplit(pi, "\\[")[[1]][1]

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

  if (length(fit@opt)) out$objective <- fit@opt$objective

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
#' @param ... Other argument to base graphics
#' @return
#' The accompanying plot function returns a line plot for a 1-dimensional profile or a contour plot for a two
#' dimensional profile.
#' @importFrom graphics contour
#' @importFrom reshape2 acast
#' @export
plot.MARSprof <- function(x, component = "objective", rel = TRUE, xlab, ylab, main, ...) {

  p1 <- attr(x, "p1")
  p2 <- attr(x, "p2")

  if (missing(xlab)) xlab <- p1
  if (is.null(x[[component]])) stop("\"", component, "\" not found in ", substitute(x))

  if (is.null(p2)) {
    xplot <- x[[p1]]
    yplot <- x[[component]]
    if (grepl("logprior", component) || grepl("loglike", component)) {
      yplot <- -1 * yplot
    }
    if (rel) yplot <- yplot - min(yplot)
    if (missing(ylab)) ylab <- paste("Change in", component)
    if (missing(main)) main <- NULL

    plot(xplot, yplot, xlab = xlab, ylab = ylab, typ = "o", main = main, ...)

  } else {

    names(x)[names(x) == p1] <- "p1"
    names(x)[names(x) == p2] <- "p2"

    zplot <- reshape2::acast(x, list("p1", "p2"), value.var = component)
    if (grepl("logprior", component) || grepl("loglike", component)) {
      zplot <- -1 * zplot
    }
    if (rel) zplot <- zplot - min(zplot)
    if (missing(ylab)) ylab <- p2
    if (missing(main)) main <- paste("Change in", response)

    contour(x = as.numeric(rownames(zplot)), y = as.numeric(colnames(zplot)),
            z = zplot, xlab = xlab, ylab = ylab, main = main, ...)
  }
  invisible()
}

