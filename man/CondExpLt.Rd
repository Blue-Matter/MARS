% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/if-statement.R
\name{CondExpLt}
\alias{CondExpLt}
\alias{CondExpLe}
\alias{CondExpGt}
\alias{CondExpGe}
\alias{CondExpEq}
\title{If statements compatible with RTMB}
\usage{
CondExpLt(left, right, if_true, if_false)

CondExpLe(left, right, if_true, if_false)

CondExpGt(left, right, if_true, if_false)

CondExpGe(left, right, if_true, if_false)

CondExpEq(left, right, if_true, if_false)
}
\arguments{
\item{left}{Numeric on left hand side of the evaluation}

\item{right}{Numeric on right hand side of the evaluation}

\item{if_true}{Numeric if expression is true}

\item{if_false}{Numeric if expression is false}
}
\description{
Convenience functions that allow taping of gradients in RTMB with if expressions,
following the corresponding \code{CppAD} functions.
}
\details{
Functions should be vectorized.

\code{CondExpLt} evaluates whether \code{left < right}

\code{CondExpLe} evaluates whether \code{left <= right}

\code{CondExpGt} evaluates whether \code{left > right}

\code{CondExpGe} evaluates whether \code{left >= right}

\code{CondExpEq} evaluates whether \code{left == right}
}
\examples{
library(RTMB)
TapeConfig(comparison = "tape")
f <- function(x) CondExpLt(x, 3, 0, x^2)
g <- MakeTape(f, numeric(1))
x <- seq(0, 5)

# Does not work!
f2 <- function(x) if (x < 3) 0 else x^2
g2 <- MakeTape(f2, numeric(1))

data.frame(x = x, deriv = sapply(x, g$jacobian), deriv2 = sapply(x, g2$jacobian))

}
