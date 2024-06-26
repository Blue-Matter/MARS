% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/report-int-state.R
\name{plot-MARS-state}
\alias{plot-MARS-state}
\alias{plot_S}
\alias{plot_B}
\alias{plot_R}
\alias{plot_SRR}
\alias{plot_Rdev}
\alias{plot_Fstock}
\alias{plot_self}
\alias{plot_seli}
\alias{plot_selstock}
\alias{plot_N}
\alias{plot_V}
\alias{plot_Ffleet}
\alias{plot_mov}
\alias{plot_recdist}
\title{Plotting functions for fitted MARS model}
\usage{
plot_S(fit, by = c("total", "stock", "region"), r, s, prop = FALSE)

plot_B(fit, by = c("total", "stock", "region"), r, s, prop = FALSE)

plot_R(fit, s)

plot_SRR(fit, s = 1, phi = TRUE)

plot_Rdev(fit, s = 1, log = TRUE)

plot_Fstock(fit, s, by = c("annual", "season"))

plot_self(fit, f = 1, type = c("length", "age"))

plot_seli(fit, i = 1)

plot_selstock(
  fit,
  s = 1,
  by = c("annual", "season"),
  plot2d = c("contour", "filled.contour"),
  ...
)

plot_N(fit, m = 1, r, s = 1, plot2d = c("contour", "filled.contour"), ...)

plot_V(fit, f = 1, by = c("stock", "region"), prop = FALSE)

plot_Ffleet(fit, f = 1)

plot_mov(fit, s = 1, y, a)

plot_recdist(fit)
}
\arguments{
\item{fit}{\linkS4class{MARSassess} object returned by \code{\link[=fit_MARS]{fit_MARS()}}}

\item{by}{Character to indicate whether to calculate selectivity from F per year or per season}

\item{r}{Integer for the corresponding region}

\item{s}{Integer for the corresponding stock}

\item{prop}{Logical, whether to plot proportions (TRUE) or absolute numbers}

\item{phi}{Logical, whether to plot unfished replacement line}

\item{log}{Logical, whether to plot the natural logarithm of the response variable}

\item{f}{Integer for the corresponding fleet}

\item{type}{For \code{plot_self}, indicates whether to plot the selectivity by age or length.}

\item{i}{Integer for the corresponding survey}

\item{plot2d}{Character, plotting function for either a \code{\link[=contour]{contour()}} or \code{\link[=filled.contour]{filled.contour()}} plot}

\item{...}{Other argument to the base graphics function}

\item{m}{Integer for the corresponding season}

\item{y}{Integer, year for plotting the movement matrix}

\item{a}{Integer, corresponding age for plotting the movement matrix}
}
\value{
Various base graphics plots
}
\description{
A set of functions to plot state variables (biomass, recruitment time series, etc.)
}
\details{
\itemize{
\item \code{plot_S} plots spawning output by stock or region (either whole numbers or proportions for the latter)
}

\itemize{
\item \code{plot_B} plots total biomass by stock or region (either whole numbers or proportions for the latter)
}

\itemize{
\item \code{plot_R} plots recruitment by stock
}

\itemize{
\item \code{plot_SRR} plots the stock-recruitment relationship and history (realized recruitment) by stock
}

\itemize{
\item \code{plot_Rdev} plots recruitment deviations by stock
}

\itemize{
\item \code{plot_Fstock} plots apical instantaneous fishing mortality (per year or per season) by stock
}

\itemize{
\item \code{plot_self} plots fishery selectivity
}

\itemize{
\item \code{plot_seli} plots index selectivity
}

\itemize{
\item \code{plot_selstock} plots the realized selectivity from total catch and total abundance at age
}

\itemize{
\item \code{plot_N} reports total abundance at age
}

\itemize{
\item \code{plot_V} plots vulnerable biomass, availability to the fishery
}

\itemize{
\item \code{plot_Ffleet} plots apical instantaneous fishing mortality (per season) by fleet
}

\itemize{
\item \code{plot_mov} plots movement matrices and the corresponding equilibrium distribution in multi-area models
}

\itemize{
\item \code{plot_recdist} plots the distribution of recruitment for each stock
}
}
