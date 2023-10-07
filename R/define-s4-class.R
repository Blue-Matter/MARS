

#' @importFrom methods slot setClass setMethod setOldClass
init_fn <- function(.Object, dots = list()) {
  if (length(dots)) {
    for(i in names(dots)) slot(.Object, i) <- dots[[i]]
  }
  attr(.Object, "version") <- paste("MARS", packageVersion("MARS"))
  attr(.Object, "date") <- date()
  attr(.Object, "R.version") <- getRversion()

  return(.Object)
}

#' Dmodel S4 object
#' @template MARSdata-template
#' @template Dmodel-slot
#' @export
setClass(
  "Dmodel",
  slots = c(ny = "numeric", nm = "numeric", na = "numeric", nl = "numeric", nr = "numeric", ns = "numeric",
            lbin = "numeric", lmid = "numeric",
            Fmax = "numeric", nitF = "numeric", dist_type = "character", y_phi = "numeric")
)

#' Dstock S4 object
#' @template MARSdata-template
#' @template Dstock-slot
#' @export
setClass(
  "Dstock",
  slots = c(len_ymas = "array", sdlen_ymas = "array",
            LAK_ymals = "array", fec_yas = "array",
            swt_ymas = "array", Md_yas = "array", m_spawn = "numeric", m_rec = "numeric", SRR_s = "character",
            delta_s = "numeric", natal_rs = "matrix")
)

#' Dfishery S4 object
#' @template MARSdata-template
#' @template Dfishery-slot
#' @export
setClass(
  "Dfishery",
  slots = c(nf = "numeric", Cobs_ymfr = "array", fwt_ymafs = "array", CAAobs_ymafr = "array", CALobs_ymlfr = "array",
            fcomp_like = "character", CAAN_ymfr = "array", CALN_ymfr = "array", CAAtheta_f = "numeric", CALtheta_f = "numeric",
            sel_block_yf = "array", sel_f = "character",
            SC_ymafrs = "array", SC_aa = "matrix", SC_ff = "matrix",
            SC_like = "character", SCN_ymafr = "array", SCtheta_f = "numeric", SCstdev_f = "numeric")
)

#' Dsurvey S4 object
#' @template MARSdata-template
#' @template Dsurvey-slot
#' @export
setClass(
  "Dsurvey",
  slots = c(ni = "numeric", Iobs_ymi = "array", unit_i = "character", IAAobs_ymai = "array", IALobs_ymli = "array",
            icomp_like = "character", IAAN_ymi = "array", IALN_ymi = "array", IAAtheta_i = "numeric", IALtheta_i = "numeric",
            samp_irs = "array", sel_i = "character", delta_i = "numeric")
)

#' DCKMR S4 object
#'
#' Store lists of data frames for parent offspring pairs and half-sibling pairs
#' @template DCKMR-slot
#' @export
setClass(
  "DCKMR",
  slots = c(POP_s = "list", HSP_s = "list", CKMR_like = "character")
)

#' @export
setClass(
  "Dtag",
  slots = c(tag_ymrr = "array", tag_ymr = "array", tag_like = "character")
)

#' MARSdata S4 object
#' @keywords MARSdata
#' @template MARSdata-template
#' @slot Dmodel Class \linkS4class{Dmodel} containing parameters for model structure (number of years, ages, etc.)
#' @slot Dstock Class \linkS4class{Dstock} containing stock parameters (growth, natural mortality, etc.)
#' @slot Dfishery Class \linkS4class{Dfishery} containing fishery data (catch, size and stock composition, etc.)
#' @slot Dsurvey Class \linkS4class{Dsurvey} containing survey data (indices of abundance)
#' @slot DCKMR Class \linkS4class{DCKMR} containing genetic close-kin data
#' @slot Dtag Class \linkS4class{Dtag} containing tagging data
#' @slot Misc List for miscellaneous inputs as needed
#' @template Dmodel-slot
#' @template Dstock-slot
#' @template Dfishery-slot
#' @template Dsurvey-slot
#' @template DCKMR-slot
#' @export
setClass(
  "MARSdata",
  slots = c(Dmodel = "Dmodel", Dstock = "Dstock", Dfishery = "Dfishery", Dsurvey = "Dsurvey",
            DCKMR = "DCKMR", Dtag = "Dtag", Misc = "list")
)
setMethod("initialize", "MARSdata", function(.Object, ...) init_fn(.Object, list(...)))
setMethod("initialize", "Dmodel", function(.Object, ...) init_fn(.Object, list(...)))
setMethod("initialize", "Dstock", function(.Object, ...) init_fn(.Object, list(...)))
setMethod("initialize", "Dfishery", function(.Object, ...) init_fn(.Object, list(...)))
setMethod("initialize", "Dsurvey", function(.Object, ...) init_fn(.Object, list(...)))
setMethod("initialize", "DCKMR", function(.Object, ...) init_fn(.Object, list(...)))
setMethod("initialize", "Dtag", function(.Object, ...) init_fn(.Object, list(...)))



setOldClass("sdreport")

#' MARSassess S4 object
#'
#' S4 object that returns output from MARS model
#'
#' @slot obj RTMB object returned by \link[RTMB]{MakeADFun}
#' @slot opt List returned by \link[stats]{nlminb}
#' @slot SD List returned by \link[RTMB]{sdreport}
#' @slot report List of model output at the parameter estimates, returned by `obj$report(obj$env$last.par.best)`
#' @slot Misc List, miscellaneous items
#'
#' @keywords MARSassess
#'
#' @export
setClass(
  "MARSassess",
  slots = c(obj = "list", opt = "list", SD = "sdreport", report = "list", Misc = "list")
)
setMethod("initialize", "MARSassess", function(.Object, ...) init_fn(.Object, list(...)))
