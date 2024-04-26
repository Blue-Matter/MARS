

#' @importFrom methods slot setClass setMethod setOldClass slot<-
#' @importFrom utils packageVersion
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
            lbin = "numeric", lmid = "numeric", Fmax = "numeric",
            y_phi = "numeric", scale_s = "numeric", nyinit = "numeric",
            condition = "character", nitF = "numeric",
            y_Fmult_f = "numeric", m_Fmult_f = "numeric", r_Fmult_f = "numeric")
)

#' Dstock S4 object
#' @template MARSdata-template
#' @template Dstock-slot
#' @export
setClass(
  "Dstock",
  slots = c(m_spawn = "numeric", m_rec = "numeric",
            len_ymas = "array", sdlen_ymas = "array",
            LAK_ymals = "array", matd_yas = "array", swt_ymas = "array", fec_yas = "array",
            Md_yas = "array", SRR_s = "character",
            delta_s = "numeric", presence_rs = "matrix", natal_rs = "matrix")
)

#' Dfishery S4 object
#' @template MARSdata-template
#' @template Dfishery-slot
#' @export
setClass(
  "Dfishery",
  slots = c(nf = "numeric", Cobs_ymfr = "array", Csd_ymfr = "array", fwt_ymafs = "array", CAAobs_ymafr = "array", CALobs_ymlfr = "array",
            fcomp_like = "character", CAAN_ymfr = "array", CALN_ymfr = "array", CAAtheta_f = "numeric", CALtheta_f = "numeric",
            sel_block_yf = "array", sel_f = "character", Cinit_mfr = "array",
            SC_ymafrs = "array", SC_aa = "matrix", SC_ff = "matrix",
            SC_like = "character", SCN_ymafr = "array", SCtheta_f = "numeric", SCstdev_f = "numeric")
)

#' Dsurvey S4 object
#' @template MARSdata-template
#' @template Dsurvey-slot
#' @export
setClass(
  "Dsurvey",
  slots = c(ni = "numeric", Iobs_ymi = "array", Isd_ymi = "array", unit_i = "character",
            IAAobs_ymai = "array", IALobs_ymli = "array",
            icomp_like = "character", IAAN_ymi = "array", IALN_ymi = "array", IAAtheta_i = "numeric", IALtheta_i = "numeric",
            samp_irs = "array", sel_i = "vector", delta_i = "numeric")
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

#' Dtag S4 object
#' @template MARSdata-template
#' @template Dtag-slot
#' @export
setClass(
  "Dtag",
  slots = c(tag_ymarrs = "array", tag_ymars = "array", tag_yy = "matrix", tag_aa = "matrix", tag_like = "character",
            tagN_ymars = "vector", tagN_ymas = "vector", tagtheta_s = "vector", tagstdev_s = "vector")
)

#' Dlabel S4 object
#'
#' Vectors for labeling plots.
#'
#' @template Dlabel-slot
#' @export
setClass(
  "Dlabel",
  slots = c(year = "vector", season = "vector", age = "vector", region = "vector",
            stock = "vector", fleet = "vector", index = "vector")
)


#' MARSdata S4 object
#' @keywords MARSdata
#' @template MARSdata-template
#' @slot Dmodel Class [Dmodel-class] containing parameters for model structure (number of years, ages, etc.)
#' @slot Dstock Class [Dstock-class] containing stock parameters (growth, natural mortality, etc.)
#' @slot Dfishery Class [Dfishery-class] containing fishery data (catch, size and stock composition, etc.)
#' @slot Dsurvey Class [Dsurvey-class] containing survey data (indices of abundance)
#' @slot DCKMR Class [DCKMR-class] containing genetic close-kin data
#' @slot Dtag Class [Dtag-class] containing tagging data
#' @slot Dlabel Class [Dlabel-class] containing names for various dimensions. Used for plotting.
#' @slot Misc List for miscellaneous inputs as needed
#' @template Dmodel-slot
#' @template Dstock-slot
#' @template Dfishery-slot
#' @template Dsurvey-slot
#' @template DCKMR-slot
#' @template Dtag-slot
#' @template Dlabel-slot
#' @export
setClass(
  "MARSdata",
  slots = c(Dmodel = "Dmodel", Dstock = "Dstock", Dfishery = "Dfishery", Dsurvey = "Dsurvey",
            DCKMR = "DCKMR", Dtag = "Dtag", Dlabel = "Dlabel", Misc = "list")
)
setMethod("initialize", "MARSdata", function(.Object, ...) init_fn(.Object, list(...)))
setMethod("initialize", "Dmodel", function(.Object, ...) init_fn(.Object, list(...)))
setMethod("initialize", "Dstock", function(.Object, ...) init_fn(.Object, list(...)))
setMethod("initialize", "Dfishery", function(.Object, ...) init_fn(.Object, list(...)))
setMethod("initialize", "Dsurvey", function(.Object, ...) init_fn(.Object, list(...)))
setMethod("initialize", "DCKMR", function(.Object, ...) init_fn(.Object, list(...)))
setMethod("initialize", "Dtag", function(.Object, ...) init_fn(.Object, list(...)))
setMethod("initialize", "Dlabel", function(.Object, ...) init_fn(.Object, list(...)))



setOldClass("sdreport")

#' MARSassess S4 object
#'
#' S4 object that returns output from MARS model
#'
#' @slot obj RTMB object returned by [RTMB::MakeADFun()]
#' @slot opt List returned by [stats::nlminb()]
#' @slot SD List returned by [TMB::sdreport()]
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

#' @export
summary.MARSassess <- function(object, ...) {
  if (length(object@SD) > 1) {
    sdreport_int(object@SD)
  } else {
    stop("No SD object found.")
  }
}


#' @importFrom utils globalVariables
if(getRversion() >= "2.15.1") {
  utils::globalVariables(
    c(slotNames("Dmodel"), slotNames("Dstock"), slotNames("Dfishery"),
      slotNames("Dsurvey"), slotNames("DCKMR"), slotNames("Dtag"), slotNames("MARSdata"), slotNames("MARSassess"))
  )
}


#' Generate markdown reports
#'
#' Generate a markdown report of model fits and estimates.
#'
#' @param object An object from MARS.
#' @param ... Additional arguments to render reports.
#'
#' @export
report <- function(object, ...) UseMethod("report")


#' @inheritParams report.MARSretro
#' @param name Optional character string for the model name to include in the report, e.g., model run number. Default
#' uses `substitute(object)`
#' @return
#' `report.MARSassess` returns a HTML markdown report.
#' @rdname report
#' @export
report.MARSassess <- function(object, name, filename = "MARS", dir = tempdir(), open_file = TRUE, render_args = list(), ...) {

  if (missing(name)) name <- substitute(object) %>% as.character()

  dots <- list(...)
  x <- object # Needed for markdown file

  dat <- get_MARSdata(object)

  nm <- dat@Dmodel@nm

  nr <- dat@Dmodel@nr
  rname <- dat@Dlabel@region
  if (!length(rname)) rname <- "Region 1"

  ns <- dat@Dmodel@ns
  sname <- dat@Dlabel@stock
  if (!length(sname)) sname <- "Stock 1"

  rmd <- system.file("include", "MARSreport.Rmd", package = "MARS") %>% readLines()
  rmd_split <- split(rmd, 1:length(rmd))

  name_ind <- grep("NAME", rmd)
  rmd_split[[name_ind]] <- paste("#", name, "{.tabset}")

  fname <- dat@Dlabel@fleet
  nf <- dat@Dfishery@nf
  fishery_ind <- grep("*ADD FISHERY RMD*", rmd)
  rmd_split[[fishery_ind]] <- mapply(make_rmd_fishery, f = 1:nf, fname = fname,
                                     MoreArgs = list(nm = nm, rname = rname)) %>%
    as.character()

  ni <- dat@Dsurvey@ni
  iname <- dat@Dlabel@index
  survey_ind <- grep("*ADD SURVEY RMD*", rmd)
  if (ni > 0) {
    rmd_split[[survey_ind]] <- mapply(make_rmd_survey, i = 1:ni, iname = iname) %>%
      as.character()
  } else {
    rmd_split[[survey_ind]] <- ""
  }

  sc_ind <- grep("*ADD SC RMD*", rmd)
  if (ns > 1 && any(dat@Dfishery@SC_ymafrs > 0, na.rm = TRUE)) {

    r_plot <- apply(dat@Dfishery@SC_ymafrs, 5, sum) > 0

    rmd_split[[sc_ind]] <- local({

      sc_txt <- lapply(1:nr, function(r) {
        if (r_plot[r]) {
          out <- lapply(1:nrow(dat@Dfishery@SC_ff), function(ff) {
            lapply(1:nrow(dat@Dfishery@SC_aa), function(aa) {
              fname <- ifelse(nrow(dat@Dfishery@SC_ff) != ncol(dat@Dfishery@SC_ff), "aggregate fleet", "fleet")
              aname <- ifelse(nrow(dat@Dfishery@SC_aa) != ncol(dat@Dfishery@SC_aa), "age class", "age")
              make_rmd_SC(ff, aa, r, paste(fname, ff), paste(aname, aa), rname[r])
            }) %>% unlist()
          })
          do.call(c, out)
        } else {
          ""
        }
      }) %>%
        unlist() %>%
        as.character()

      c("## Stock composition {.tabset}\n\n", sc_txt)
    })

  } else {
    rmd_split[[sc_ind]] <- ""
  }

  tagmov_ind <- grep("*ADD TAG MOV RMD*", rmd)
  if (nr > 1 && any(dat@Dtag@tag_ymarrs > 0, na.rm = TRUE)) {

    r_plot <- apply(dat@Dfishery@SC_ymafrs, 5, sum) > 0

    rmd_split[[tagmov_ind]] <- local({

      tagmov_txt <- lapply(1:ns, function(s) {
        out <- lapply(1:nrow(dat@Dtag@tag_yy), function(yy) {
          lapply(1:nrow(dat@Dtag@tag_aa), function(aa) {
            yname <- ifelse(nrow(dat@Dtag@tag_yy) != ncol(dat@Dtag@tag_yy), "aggregate year", "year")
            aname <- ifelse(nrow(dat@Dtag@tag_aa) != ncol(dat@Dtag@tag_aa), "age class", "age")
            make_rmd_tagmov(yy, aa, s, paste(yname, yy), paste(aname, aa), sname[s], header = yy == 1 && aa == 1)
          }) %>% unlist()
        })
        do.call(c, out)
      }) %>%
        unlist() %>%
        as.character()

      c("## Tag movement {.tabset}\n\n", tagmov_txt)
    })

  } else {
    rmd_split[[tagmov_ind]] <- ""
  }

  srr_ind <- grep("*ADD SRR RMD*", rmd)
  rmd_split[[srr_ind]] <- mapply(make_rmd_srr, s = 1:ns, sname = sname) %>% as.character()

  stock_ind <- grep("*ADD STOCK REGION RMD*", rmd)
  mov_ind <- grep("*ADD MOVEMENT RMD*", rmd)
  if (nr > 1) {
    year <- dat@Dlabel@year
    rmd_split[[stock_ind]] <- mapply(make_rmd_stock_region, s = 1:ns, sname = sname) %>% as.character()

    if (is.null(dots$ymov)) dots$ymov <- dat@Dmodel@ny
    if (is.null(dots$amov)) dots$amov <- 2

    mov <- expand.grid(a = dots$amov, y = dots$ymov, s = 1:ns)
    rmd_split[[mov_ind]] <- sapply(1:nrow(mov), function(i) {
      make_rmd_mov(s = mov$s[i], y = mov$y[i], a = mov$a[i], yname = dat@Dlabel@year[mov$y[i]],
                   sname = sname[mov$s[i]], header = i == 1)
    }) %>% as.character()

  } else {
    rmd_split[[stock_ind]] <- rmd_split[[mov_ind]] <- ""
  }

  ####### Function arguments for rmarkdown::render
  filename_rmd <- paste0(filename, ".Rmd")

  render_args$input <- file.path(dir, filename_rmd)
  if (is.null(render_args$quiet)) render_args$quiet <- TRUE

  # Generate markdown report
  if (!dir.exists(dir)) {
    message_info("Creating directory: ", dir)
    dir.create(dir)
  }
  write(unlist(rmd_split), file = file.path(dir, filename_rmd))

  # Rendering markdown file
  message_info("Rendering markdown file: ", file.path(dir, filename_rmd))
  output_filename <- do.call(rmarkdown::render, render_args)
  message("Rendered file: ", output_filename)

  if (open_file) browseURL(output_filename)
  invisible(output_filename)
}

