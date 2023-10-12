
#' @importFrom graphics grid matplot abline par
plot.default <- function(..., zero_line = FALSE, mar = c(5, 4, 1, 1)) {
  old_mar <- par()$mar
  on.exit(par(mar = old_mar))
  par(mar = mar)

  if (zero_line) {
    graphics::plot.default(..., panel.first = {graphics::grid(); abline(h = 0, col = "grey60")})
  } else {
    graphics::plot.default(..., panel.first = graphics::grid())
  }
}

matplot <- function(..., zero_line = FALSE, mar = c(5, 4, 1, 1)) {
  old_mar <- par()$mar
  on.exit(par(mar = old_mar))
  par(mar = mar)

  if (zero_line) {
    graphics::matplot(..., panel.first = {graphics::grid(); abline(h = 0, col = "grey60")})
  } else {
    graphics::matplot(..., panel.first = graphics::grid())
  }
}
