
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


show <- function(object) methods::show(object)
setMethod("show",
          signature = "MARSassess",
          function(object) {

            cat("Number of parameters:", length(object@obj$par))
            if (length(object@SD) && !is.null(object@SD$gradient.fixed)) {
              gr <- abs(object@SD$gradient.fixed)
            } else {
              gr <- abs(object@obj$gr())
            }
            cat("\nMaximum gradient:", round(max(gr), 4))

            if (length(object@SD) && !is.null(object@SD$gradient.fixed)) {
              gr_large <- gr > 0.1

              if (sum(gr_large)) {
                gr_report <- gr[gr_large]
                gr_names <- make_unique_names(object@SD)[gr_large]

                cat("\nParameters with large gradients (> 0.1):\n")

                gr_order <- order(gr_report, decreasing = TRUE)

                x <- data.frame(Gradient = round(gr_report, 4))
                rownames(x) <- gr_names
                print(x[gr_order, "Gradient", drop = FALSE])
              }
            }

            if (length(object@SD) && !is.null(object@SD$env$hessian)) {
              h <- object@SD$env$hessian
              det_h <- det(h)
              cat(paste("\nDeterminant of Hessian:"), round(det_h, 4))

              zero_rows <- apply(h, 1, function(x) all(x == 0))
              if (any(zero_rows)) {
                cat("\nParameters with all zeros in Hessian:\n")
                par_zero <- names(zero_rows)[zero_rows]
                for(i in par_zero) cat(i, "\n")
              }
            }

            invisible()
          })
