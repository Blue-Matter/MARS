
#' @import RTMB
#' @importFrom dplyr %>%
NULL

#' @importFrom methods callNextMethod
setMethod(
  "sapply", signature(X="integer"),
  function (X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE) {
    ans <- callNextMethod()
    if (is.complex(ans))
      class(ans) <- "advector"
    ans
  }
)


#' `sapply2` function
#'
#' An alternate `sapply` function with argument `simplify = "array"` for convenience.
#'
#' @param X,FUN,...,USE.NAMES Same arguments as [sapply()]
#' @export
#' @keywords internal
sapply2 <- function(X, FUN, ..., USE.NAMES = TRUE) {
  sapply(X, FUN, ..., simplify = "array", USE.NAMES = USE.NAMES)
}

#' @name AD
#' @title Additional methods for AD types
#'
#' @description Methods for RTMB AD class
NULL

#' @describeIn AD Matrix product function implemented for mixed AD and non-AD objects with `colSums(x * y)`. See \link[RTMB]{ADmatrix}.
#'
#' @param x AD object
#' @param y Non-AD matrix
#' @aliases %*%,ad,matrix-method matmul
#' @keywords internal
setMethod("%*%",
          signature("ad", "matrix"),
          function(x, y) {
            colSums(x * y)
          })

#' @describeIn AD Finds the maximum value using [CondExpGt()] over a loop
#' @aliases max.advector max
#' @param ... Objects of class advector
#' @param na.rm Not used
#' @export
#' @export max.advector
max.advector <- function(..., na.rm) {
  oldval <- TapeConfig()["comparison"]
  on.exit(TapeConfig(comparison = oldval))
  TapeConfig(comparison = "allow")

  dots <- list(...)
  #c <- ADoverload(x = "c")
  x <- do.call(c, dots)

  xout <- x[1]
  for (i in 2:length(x)) xout <- CondExpGt(xout, x[i-1], xout, x[i-1])
  asS4(xout)
}

#' @describeIn AD Finds the minimum value using [CondExpLt()] over a loop
#' @aliases min.advector min
#' @export
#' @export min.advector
min.advector <- function(..., na.rm) {
  oldval <- TapeConfig()["comparison"]
  on.exit(TapeConfig(comparison = oldval))
  TapeConfig(comparison = "allow")

  dots <- list(...)
  #c <- ADoverload(x = "c")
  x <- do.call(c, dots)

  xout <- x[1]
  for (i in 2:length(x)) xout <- CondExpLt(xout, x[i-1], xout, x[i-1])
  asS4(xout)
}
