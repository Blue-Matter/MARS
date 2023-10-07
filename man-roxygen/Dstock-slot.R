#' @section Slots inherited from Dstock:
#' \describe{
#' \item{`len_ymas`}{Length-at-age. Only needed if `Dmodel@nl > 0`.}
#' \item{`sdlen_ymas`}{Standard deviation in length-at-age}
#' \item{`LAK_ymals`}{Length-at-age probability array. If empty, values will be calculated in \link{check_data} (see \link{calc_LAK}).}
#' \item{`mat_yas`}{Proportion mature by age class}
#' \item{`fec_yas`}{Fecundity, i.e., spawning output, of mature animals}
#' \item{`swt_ymas`}{Stock weight-at-age}
#' \item{`Md_yas`}{Natural morality. Ignored if M is estimated.}
#' \item{`m_spawn`}{Integer, season of spawning}
#' \item{`m_rec`}{Integer, season of recruitment}
#' \item{`SRR_s`}{Character vector of stock-recruit relationship by stock. See \link{calc_recruitment} for options.}
#' \item{`delta_s`}{Fraction of season that elapses when spawning occurs}
#' \item{`natal_rs`}{Boolean matrix that indicates whether stock `s` spawns when in region `r` at time of spawning.}
#' }
