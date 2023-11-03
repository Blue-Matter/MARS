
#' @importFrom graphics grid matplot abline par
plot.default <- function(..., zero_line = FALSE) {
  if (zero_line) {
    graphics::plot.default(..., panel.first = {graphics::grid(); abline(h = 0, col = "grey60")})
  } else {
    graphics::plot.default(..., panel.first = graphics::grid())
  }
}

matplot <- function(..., zero_line = FALSE) {
  if (zero_line) {
    graphics::matplot(..., panel.first = {graphics::grid(); abline(h = 0, col = "grey60")})
  } else {
    graphics::matplot(..., panel.first = graphics::grid())
  }
}
