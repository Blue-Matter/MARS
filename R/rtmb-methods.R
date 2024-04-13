
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

#' @name matmul
#' @title Matrix multiplication
#'
#' @description Matrix product function implemented for mixed AD and non-AD objects. See \link[RTMB]{ADmatrix}.
#'
#' @param x AD object
#' @param y Matrix
#' @aliases %*%,ad,matrix-method
#' @keywords internal
setMethod("%*%",
          signature("ad", "matrix"),
          function(x, y) {
            colSums(x * y)
          })
