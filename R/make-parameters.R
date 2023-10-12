

check_Dmodel <- function(Dmodel, silent = FALSE) {
  getAllS4(Dmodel)

  ch <- as.character(substitute(Dmodel))
  if (length(ch) > 1) ch <- "Dmodel"

  if (!length(ny)) stop("Need ", ch, "@ny")
  if (!length(nm)) {
    if (!silent) message("Setting ", ch, "@nm to 1")
    Dmodel@nm <- 1L
  }
  if (!length(na)) stop("Need ", ch, "@na")
  if (!length(nl)) {
    if (!silent) message("No length bins in the model.")
    Dmodel@nl <- 0L
  }
  if (!length(nr)) {
    if (!silent) message("Creating one region model.")
    Dmodel@nr <- 1L
  }
  if (!length(ns)) {
    if (!silent) message("Creating one stock model.")
    Dmodel@ns <- 1L
  }
  if (Dmodel@nl > 0) {
    if (length(lbin) != nl+1) stop("Need nl + 1 vector for ", ch, "@lbin")
    if (length(lmid) != nl) stop("Need nl vector for ", ch, "@lmid")
  }
  if (!length(Fmax)) {
    if (!silent) message("Setting ", ch, "@Fmax to 3 (per season)")
    Dmodel@Fmax <- 3
  }
  if (!length(nitF)) {
    if (!silent) message("Setting ", ch, "@nitF to 5 (per season)")
    Dmodel@nitF <- 5
  }
  if (Dmodel@nr > 0 && !length(dist_type)) {
    if (!silent) message("Setting ", ch, "@dist_type to \"dist\"")
    Dmodel@dist_type <- "dist"
  }
  if (!length(y_phi)) {
    if (!silent) message("Setting ", ch, "@y_phi to year 1")
    Dmodel@y_phi <- 1L
  }
  return(Dmodel)
}

check_Dstock <- function(Dstock, Dmodel, silent = FALSE) {
  getAllS4(Dstock, Dmodel)

  ch <- as.character(substitute(Dstock))
  if (length(ch) > 1) ch <- "Dstock"

  if (nl > 0) {
    if (length(LAK_ymals)) {
      dim_LAK <- dim(LAK_ymals) == c(ny, nm, na, nl, ns)
      if (!all(dim_LAK)) stop("dim(LAK_ymals) needs to be: ", c(ny, nm, na, nl, ns) %>% paste(collapse = ", "))
    } else {

      dim_len <- dim(len_ymas) == c(ny, nm, na, ns)
      if (!all(dim_len)) stop("dim(len_ymas) needs to be: ", c(ny, nm, na, ns) %>% paste(collapse = ", "))

      dim_sdlen <- dim(sdlen_ymas) == c(ny, nm, na, ns)
      if (!all(dim_sdlen)) stop("dim(sdlen_ymas) needs to be: ", c(ny, nm, na, ns) %>% paste(collapse = ", "))

      if (!silent) message("Calculating ", ch, "@LAK_ymals array")
      Dstock@LAK_ymals <- sapply2(1:ns, function(s) {
        sapply2(1:nm, function(m) {
          sapply2(1:ny, function(y) calc_LAK(len_ymas[y, m, , s], sdlen_ymas[y, m, , s], lbin))
        })
      }) %>% aperm(c(3, 4, 1, 2, 5))
    }
  }

  dim_mat <- dim(mat_yas) == c(ny, na, ns)
  if (!all(dim_mat)) stop("dim(mat_yas) needs to be: ", c(ny, na, ns) %>% paste(collapse = ", "))

  dim_fec <- dim(fec_yas) == c(ny, na, ns)
  if (!all(dim_fec)) stop("dim(fec_yas) needs to be: ", c(ny, na, ns) %>% paste(collapse = ", "))

  dim_swt <- dim(swt_ymas) == c(ny, nm, na, ns)
  if (!all(dim_swt)) stop("dim(swt_ymas) needs to be: ", c(ny, nm, na, ns) %>% paste(collapse = ", "))

  dim_M <- dim(Md_yas) == c(ny, na, ns)
  if (!all(dim_M)) stop("dim(Md_yas) needs to be: ", c(ny, na, ns) %>% paste(collapse = ", "))

  if (!length(m_spawn) || nm == 1L) {
    if (!silent) message("m_spawn set to 1")
  } else if (m_spawn > nm) {
    stop("m_spawn cannot be greater than nm")
  }
  if (!length(m_rec) || nm == 1L) {
    if (!silent) message("m_rec set to 1")
  } else if (m_rec > nm) {
    stop("m_rec cannot be greater than nm")
  }

  if (length(SRR_s) != ns) stop("SRR_s needs to be length ", ns)
  if (length(delta_s) != ns) stop("delta_s needs to be length ", ns)
  if (any(dim(natal_rs) != c(nr, ns))) stop("dim(natal_rs) needs to be: ", c(nr, ns) %>% paste(collapse = ", "))

  return(Dstock)
}

check_Dfishery <- function(Dfishery, Dstock, Dmodel, silent = FALSE) {
  getAllS4(Dfishery)
  ch <- as.character(substitute(Dfishery))
  if (length(ch) > 1) ch <- "Dfishery"

  if (!length(nf)) stop("Need ", ch, "@nf")

  dim_Cobs <- dim(Cobs_ymfr) == c(ny, nm, nf, nr)
  if (!all(dim_Cobs)) stop("dim(Cobs_ymfr) needs to be: ", c(ny, nm, nf, nr) %>% paste(collapse = ", "))

  if (!length(dim_fwt)) {
    if (!silent) message("Setting fishery weight at age to stock weight at age")
    Dfishery@fwt_ymafs <- sapply2(1:ns, function(s) {
      sapply2(1:nf, function(f) Dstock@swt_ymas[, , , s])
    })
  } else {
    dim_fwt <- dim(fwt_ymafs) == c(ny, nm, na, nf, ns)
    if (!all(dim_fwt)) stop("dim(fwt_ymafs) needs to be: ", c(ny, nm, na, nf, ns) %>% paste(collapse = ", "))
  }

  if (length(CAAobs_ymafr) || length(CALobs_ymlfr)) {
    if (length(fcomp_like)) {
      fcomp_like <- match.arg(fcomp_like, choices = eval(formals(like_comp)$type))
    } else {
      if (!silent) message("Setting ", ch, "@fcomp_like = \"multinomial\"")
      Dfishery@fcomp_like <- "multinomial"
    }
  }

  if (length(CAAobs_ymafr)) {
    dim_CAA <- dim(CAAobs_ymafr) == c(ny, nm, na, nf, nr)
    if (!all(dim_CAA)) stop("dim(CAAobs_ymafr) needs to be: ", c(ny, nm, na, nf, nr) %>% paste(collapse = ", "))

    if (!length(CAAN_ymfr)) {
      if (fcomp_like %in% c("multinomial", "ddirmult1", "ddirmult2") && !silent) {
        message("Setting ", ch, "@CAAN_ymfr from ", ch, "@CAAobs_ymafr")
      }
      Dfishery@CAAN_ymfr <- apply(CAAobs_ymafr, c(1, 2, 4, 5), sum)
    } else {
      dim_CAAN <- dim(CAAN_ymfr) == c(ny, nm, nf, nr)
      if (!all(dim_CAAN)) stop("dim(CAAN_ymfr) needs to be: ", c(ny, nm, nf, nr) %>% paste(collapse = ", "))
    }

    if (!length(CAAtheta_f)) {
      if (grepl("ddirmult", fcomp_like) && !silent) message("Setting ", ch, "@CAAtheta_f to 1 for all fleets")
      Dfishery@CAAtheta_f <- rep(1, nf)
    } else if (length(CAAtheta_f) == 1) {
      Dfishery@CAAtheta_f <- rep(Dfishery@CAAtheta_f, nf)
    } else if (length(CAAtheta_f) != nf) {
      stop("Vector CAAtheta_f needs to be length ", nf)
    }
  }

  if (length(CALobs_ymlfr)) {
    if (!nl) stop("Fishery length composition detected but no length bins found in Dmodel@nl")

    dim_CAL <- dim(CALobs_ymlfr) == c(ny, nm, nl, nf, nr)
    if (!all(dim_CAL)) stop("dim(CALobs_ymlfr) needs to be: ", c(ny, nm, nl, nf, nr) %>% paste(collapse = ", "))

    if (!length(CALN_ymfr)) {
      if (fcomp_like %in% c("multinomial", "ddirmult1", "ddirmult2") && !silent) {
        message("Setting ", ch, "@CALN_ymfr from ", ch, "@CALobs_ymlfr")
      }
      Dfishery@CALN_ymfr <- apply(CALobs_ymlfr, c(1, 2, 4, 5), sum)
    } else {
      dim_CALN <- dim(CALN_ymfr) == c(ny, nm, nf, nr)
      if (!all(dim_CALN)) stop("dim(CALN_ymfr) needs to be: ", c(ny, nm, nf, nr) %>% paste(collapse = ", "))
    }

    if (!length(CALtheta_f)) {
      if (grepl("ddirmult", fcomp_like) && !silent) message("Setting ", ch, "@CALtheta_f to 1 for all fleets")
      Dfishery@CALtheta_f <- rep(1, nf)
    } else if (length(CALtheta_f) == 1) {
      Dfishery@CALtheta_f <- rep(Dfishery@CALtheta_f, nf)
    } else if (length(CAAtheta_f) != nf) {
      stop("Vector CAAtheta_f needs to be length ", nf)
    }
  }

  if (!length(sel_block_yf)) {
    Dfishery@sel_block_yf <- matrix(1:nf, ny, nf, byrow = TRUE)
  } else {
    dim_sb <- dim(sel_block_yf) == c(ny, nf)
    if (!all(dim_sb)) stop("dim(sel_block_yf) needs to be: ", c(ny, nf) %>% paste(collapse = ", "))
  }
  nb <- max(Dfishery@sel_block_yf)
  if (length(sel_f) != nb) stop("Vector sel_f should be length ", nf)

  if (!length(SC_ymafrs)) {
    dim_SC <- dim(SC_ymafrs)
    if (length(dim_SC) != 6) stop("SC_ymafrs should be a six dimensional array")

    if (any(dim_SC[c(1, 2, 5, 6)] != c(ny, nm, nr, ns))) {
      stop("dim(SC_ymafrs) should be ", c(ny, nm, dim_SC[3:4], nr, ns) %>% paste(collapse = ", "))
    }

    if (dim_SC[3] == na) {
      if (!length(SC_aa)) Dfishery@SC_aa <- diag(1, na)
    } else if (any(dim(SC_aa) != c(dim_SC[3], na))) {
      stop("dim(SC_aa) should be: ", c(dim_SC[3], na) %>% paste(collapse = ", "))
    }

    if (dim_SC[4] == nf) {
      if (!length(SC_ff)) Dfishery@SC_ff <- diag(1, nf)
    } else if (any(dim(SC_ff) != c(dim_SC[4], nf))) {
      stop("dim(SC_ff) should be: ", c(dim_SC[4], nf) %>% paste(collapse = ", "))
    }

    if (!length(SC_like)) {
      if (!silent) message("Setting ", ch, "@SC_like = \"multinomial\"")
      Dfishery@SC_like <- "multinomial"
    } else {
      SC_like <- match.arg(fcomp_like, choices = eval(formals(like_comp)$type))
    }

    if (!length(SCN_ymafr)) {
      if (SC_like %in% c("multinomial", "ddirmult1", "ddirmult2") && !silent) {
        message("Setting ", ch, "@SCN_ymafr from ", ch, "@SC_ymafrs")
      }
      Dfishery@SCN_ymafr <- apply(SC_ymafrs, 1:5, sum)
    } else {
      dim_SCN <- dim(SCN_ymafr) == c(ny, nm, dim_SC[3], dim_SC[4], nr)
      if (!all(dim_SCN)) stop("dim(SCN_ymafr) needs to be: ", c(ny, nm, dim_SC[3], dim_SC[4], nr) %>% paste(collapse = ", "))
    }

    if (!length(SCtheta_f)) {
      if (grepl("ddirmult", SC_like) && !silent) message("Setting ", ch, "@SCtheta_f to 1 for all fleets")
      Dfishery@SCtheta_f <- rep(1, dim_SC[4])
    } else if (length(SCtheta_f) == 1) {
      Dfishery@SCtheta_f <- rep(Dfishery@SCtheta_f, dim_SC[4])
    } else if (length(SCtheta_f) != dim_SC[4]) {
      stop("Vector SCtheta_f needs to be length ", dim_SC[4])
    }

    if (!length(SCstdev_f)) {
      if (grepl("log", SC_like) && !silent) message("Setting ", ch, "@SCstdev_f to 0.1 for all fleets")
      Dfishery@SCstdev_f <- rep(0.1, dim_SC[4])
    } else if (length(SCstdev_f) == 1) {
      Dfishery@SCstdev_f <- rep(Dfishery@SCstdev_f, dim_SC[4])
    } else if (length(SCstdev_f) != dim_SC[4]) {
      stop("Vector SCstdev_f needs to be length ", dim_SC[4])
    }
  }
  return(Dfishery)
}

check_Dsurvey <- function(Dsurvey, Dmodel, silent = FALSE) {
  getAllS4(Dsurvey, Dmodel)

  ch <- as.character(substitute(Dsurvey))
  if (length(ch) > 1) ch <- "Dsurvey"

  if (!length(ni)) {
    if (length(Iobs_ymi) || length(IAAobs_ymai) || length(IALobs_ymli)) {
      stop("Need ", ch, "@ni")
    } else {
      Dsurvey@ni <- ni <- 0
    }
  }

  if (ni > 0) {
    dim_Iobs <- dim(Iobs_ymi) == c(ny, nm, ni)
    if (!all(dim_Iobs)) stop("dim(Iobs_ymi) needs to be: ", c(ny, nm, ni) %>% paste(collapse = ", "))

    if (!length(unit_i)) {
      if (!silent) message("unit_i set to biomass for all indices")
      Dsurvey@unit_i <- rep("B", ni)
    }

    if (length(IAAobs_ymai) || length(IALobs_ymli)) {
      if (length(icomp_like)) {
        icomp_like <- match.arg(icomp_like, choices = eval(formals(like_comp)$type))
      } else {
        if (!silent) message("Setting ", ch, "@icomp_like = \"multinomial\"")
        Dsurvey@icomp_like <- "multinomial"
      }
    }

    if (length(IAAobs_ymai)) {
      dim_IAA <- dim(IAAobs_ymai) == c(ny, nm, na, ni)
      if (!all(dim_IAA)) stop("dim(IAAobs_ymai) needs to be: ", c(ny, nm, na, ni) %>% paste(collapse = ", "))

      if (!length(IAAN_ymi)) {
        if (icomp_like %in% c("multinomial", "ddirmult1", "ddirmult2") && !silent) {
          message("Setting ", ch, "@IAAN_ymi from ", ch, "@IAAobs_ymai")
        }
        Dfishery@IAAN_ymi <- apply(IAAobs_ymai, c(1, 2, 4), sum)
      } else {
        dim_IAAN <- dim(IAAN_ymi) == c(ny, nm, ni)
        if (!all(dim_IAAN)) stop("dim(IAAN_ymi) needs to be: ", c(ny, nm, ni) %>% paste(collapse = ", "))
      }

      if (!length(IAAtheta_i)) {
        if (grepl("ddirmult", icomp_like) && !silent) message("Setting ", ch, "@IAAtheta_i to 1 for all indices")
        Dsurvey@IAAtheta_i <- rep(1, ni)
      } else if (length(IAAtheta_i) == 1) {
        Dsurvey@IAAtheta_i <- rep(Dsurvey@IAAtheta_i, ni)
      } else if (length(IAAtheta_i) != ni) {
        stop("Vector IAAtheta_i needs to be length ", ni)
      }
    }

    if (length(IALobs_ymli)) {
      if (!nl) stop("Index length composition detected but no length bins found in Dmodel@nl")

      dim_IAL <- dim(IALobs_ymli) == c(ny, nm, nl, ni)
      if (!all(dim_IAL)) stop("dim(IALobs_ymli) needs to be: ", c(ny, nm, nl, ni) %>% paste(collapse = ", "))

      if (!length(IALN_ymi)) {
        if (icomp_like %in% c("multinomial", "ddirmult1", "ddirmult2") && !silent) {
          message("Setting ", ch, "@IALN_ymi from ", ch, "@IALobs_ymli")
        }
        Dsurvey@IALN_ymi <- apply(IALobs_ymli, c(1, 2, 4), sum)
      } else {
        dim_IALN <- dim(IALN_ymi) == c(ny, nm, ni)
        if (!all(dim_IALN)) stop("dim(IALN_ymi) needs to be: ", c(ny, nm, ni) %>% paste(collapse = ", "))
      }

      if (!length(IALtheta_i)) {
        if (grepl("ddirmult", icomp_like) && !silent) message("Setting ", ch, "@IALtheta_i to 1 for all indices")
        Dsurvey@IALtheta_i <- rep(1, ni)
      } else if (length(IALtheta_i) == 1) {
        Dsurvey@IALtheta_i <- rep(Dsurvey@IALtheta_i, ni)
      } else if (length(IALtheta_i) != ni) {
        stop("Vector IALtheta_i needs to be length ", ni)
      }
    }

    if (!length(samp_irs)) {
      if (!silent) message("Setting samp_irs such that all indices operate in all regions and sample all stocks")
      Dsurvey@samp_irs <- array(1, c(ni, nr, ns))
    } else {
      dim_samp <- dim(samp_irs) == c(ni, nr, ns)
      if (!all(dim_samp)) stop("dim(samp_irs) needs to be: ", c(ni, nr, ns) %>% paste(collapse = ", "))
    }

    if (length(sel_i) == 1) {
      Dsurvey@sel_i <- rep(sel_i, ni)
    } else if(length(sel_i) != ni) {
      stop("Vector sel_i needs to be length ", ni)
    }

    if (!length(delta_i)) {
      if (!silent) message("Setting delta_i = 0 for all indices")
      Dsurvey@delta_i <- rep(0, ni)
    } else if (length(delta_i) == 1) {
      Dsurvey@delta_i <- rep(delta_i, ni)
    } else if(length(delta_i) != ni) {
      stop("Vector delta_i needs to be length ", ni)
    }
  }

  return(Dsurvey)
}

check_DCKMR <- function(DCKMR, Dmodel, silent = FALSE) {
  getAllS4(DCKMR, Dmodel)

  if (!length(POP_s)) {
    vars_POP <- c("a", "t", "y", "n", "m")
    check_POP <- sapply(1:ns, function(s) {
      if (nrow(POP_s[[s]])) {
        all(vars_POP %in% names(POP_s[[s]]))
      } else {
        TRUE
      }
    })
    if (any(!check_POP)) {
      stop("Missing columns in close-kin POP data frames. See: help(\"MARSdata-class\")")
    }
  }

  if (!length(HSP_s)) {
    vars_HSP <- c("yi", "yj", "n", "m")
    check_HSP <- sapply(1:ns, function(s) {
      if (nrow(HSP_s[[s]])) {
        all(vars_HSP %in% names(HSP_s[[s]]))
      } else {
        TRUE
      }
    })
    if (any(!check_HSP)) {
      stop("Missing columns in close-kin HSP data frames. See: help(\"MARSdata-class\")")
    }
  }

  if (!length(POP_s) || !length(HSP_s)) {
    if (!length(CKMR_like)) {
      if (!silent) message("Setting close-kin likelihood to \"binomial\"")
      DCKMR@CKMR_like <- "binomial"
    }
  }
  return(DCKMR)
}

check_Dtag <- function(Dtag, silent = FALSE) {

  return(Dtag)
}

#' Check dimensions and inputs in MARSdata object
#'
#' Ensures that data inputs are of proper dimension. Whenever possible, default values are added to missing items.
#'
#' @param MARSdata S4 object containing data inputs. See \linkS4class{MARSdata}
#' @param silent Logical, whether or not to report default values to the console
#' @returns An updated \linkS4class{MARSdata} object.
#' @seealso \link{MARSdata-class}
#' @export
check_data <- function(MARSdata, silent = FALSE) {
  MARSdata@Dmodel <- check_Dmodel(MARSdata@Dmodel, silent)
  MARSdata@Dstock <- check_Dstock(MARSdata@Dstock, MARSdata@Dmodel, silent)
  MARSdata@Dfishery <- check_Dfishery(MARSdata@Dfishery, MARSdata@Dstock, MARSdata@Dmodel, silent)
  MARSdata@Dsurvey <- check_Dsurvey(MARSdata@Dsurvey, MARSdata@Dstock, silent)

  MARSdata@DCKMR <- check_DCKMR(MARSdata@DCKMR, MARSdata@Dmodel, silent)
  #MARSdata@Dtag <- check_Dtag(MARSdata@Dtag, MARSdata@Dmodel, silent)
  return(MARSdata)
}

make_parameters <- function(MARSdata, start = list()) {

}

make_map <- function() {

}

check_parameters <- function(p = list()) {

}
