% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/define-s4-class.R
\docType{class}
\name{MARSassess-class}
\alias{MARSassess-class}
\title{MARSassess S4 object}
\description{
S4 object that returns output from MARS model
}
\section{Slots}{

\describe{
\item{\code{obj}}{RTMB object returned by \link[RTMB]{MakeADFun}}

\item{\code{opt}}{List returned by \link[stats]{nlminb}}

\item{\code{SD}}{List returned by \link[RTMB]{sdreport}}

\item{\code{report}}{List of model output at the parameter estimates, returned by \code{obj$report(obj$env$last.par.best)}}

\item{\code{Misc}}{List, miscellaneous items}
}}

\keyword{MARSassess}