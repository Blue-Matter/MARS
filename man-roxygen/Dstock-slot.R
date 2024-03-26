#' @section Slots inherited from Dstock:
#' \describe{
#' \item{`m_spawn`}{Integer, season of spawning. Defaults to 1.}
#' \item{`m_rec`}{Integer, season of recruitment. Defaults to 1.}
#' \item{`len_ymas`}{Length-at-age. Only needed if `Dmodel@nl > 0`. [calc_growth()] may be a helpful function.}
#' \item{`sdlen_ymas`}{Standard deviation in length-at-age}
#' \item{`LAK_ymals`}{Length-at-age probability array. If empty, values will be calculated by [check_data()] with [calc_LAK()].}
#' \item{`matd_yas`}{Proportion mature by age class. Ignored if maturity ogive is estimated, e.g., when fitting to close-kin genetic data.}
#' \item{`swt_ymas`}{Stock weight-at-age. See [calc_growth()] example.}
#' \item{`fec_yas`}{Fecundity, i.e., spawning output, of mature animals. Default uses stock weight at age.}
#' \item{`Md_yas`}{Natural mortality. Ignored if M is estimated.}
#' \item{`SRR_s`}{Character vector of stock-recruit relationship by stock. See `SRR` argument in [calc_recruitment()] for options.}
#' \item{`delta_s`}{Fraction of season that elapses when spawning occurs, e.g., midseason spawning occurs when `delta_s = 0.5`. Default is zero.}
#' \item{`natal_rs`}{The fraction of the mature stock `s` in region `r` that spawns at
#' time of spawning. See example. Default is 1 for all stocks and regions.}
#' }
#' @examples
#' # Set natal_rs matrix so that the spawning output of stock 1 is
#' # calculated from mature animals present in regions 1, 2.
#' # Similarly for stock 2, spawning output from areas 2 and 3.
#' nr <- 4
#' ns <- 2
#' natal_rs <- matrix(0, nr, ns)
#' natal_rs[1:2, 1] <- natal_rs[2:3, 2] <- 1
#'
