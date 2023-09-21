
#' If statements compatible with RTMB
#'
#' Convenience functions that allow taping of gradients in RTMB with if expressions,
#' following the corresponding `CppAD` functions.
#'
#' @param left Numeric on left hand side of the evaluation
#' @param right Numeric on right hand side of the evaluation
#' @param if_true Numeric if expression is true
#' @param if_false Numeric if expression is false
#' @details Functions should be vectorized.
#'
#' `CondExpLt` evaluates whether `left < right`
#' @aliases CondExpLe CondExpGt CondExpGe CondExpEq
#'
#' @examples
#' library(RTMB)
#' TapeConfig(comparison = "tape")
#' f <- function(x) CondExpLt(x, 3, 0, x^2)
#' g <- MakeTape(f, numeric(1))
#' x <- seq(0, 5)
#'
#' # Does not work!
#' f2 <- function(x) if (x < 3) 0 else x^2
#' g2 <- MakeTape(f2, numeric(1))
#'
#' data.frame(x = x, deriv = sapply(x, g$jacobian), deriv2 = sapply(x, g2$jacobian))
#'
#' @export
CondExpLt <- function(left, right, if_true, if_false) {
  (left < right) * if_true + (left >= right) * if_false
}

#' @rdname CondExpLt
#' @details `CondExpLe` evaluates whether `left <= right`
#' @export
CondExpLe <- function(left, right, if_true, if_false) {
  (left <= right) * if_true + (left > right) * if_false
}

#' @rdname CondExpLt
#' @details `CondExpGt` evaluates whether `left > right`
#' @export
CondExpGt <- function(left, right, if_true, if_false) {
  (left > right) * if_true + (left <= right) * if_false
}

#' @rdname CondExpLt
#' @details `CondExpGe` evaluates whether `left >= right`
#' @export
CondExpGe <- function(left, right, if_true, if_false) {
  (left >= right) * if_true + (left < right) * if_false
}

#' @rdname CondExpLt
#' @details `CondExpEq` evaluates whether `left == right`
#' @export
CondExpEq <- function(left, right, if_true, if_false) {
  (left == right) * if_true + (left != right) * if_false
}
