
make_rmd_CAA <- function(f, r, fname, rname) {
  cap <- paste0("Observed (black) and predicted (red) catch at age for ", fname, " in ", rname)
  cap2 <- paste0("Pearson residuals for the catch at age for ", fname, " in ", rname)
  cap3 <- paste0("Residual histogram for the catch at age for ", fname, " in ", rname)

  rmd <- c(
    paste0("```{r caa-f", f, "-r", r, ", fig.cap=\"", cap, ".\"}"),
    paste0("if (length(dat@Dfishery@CAAobs_ymafr) && any(dat@Dfishery@CAAobs_ymafr[, , , ", f, ", ", r, "] > 0, na.rm = TRUE)) {"),
    paste0("  plot_CAA(x, f = ", f, ", r = ", r, ")"),
    "}",
    "```",
    "",
    paste0("```{r caa-resid-f", f, "-r", r, ", fig.cap=\"", cap2, ".\"}"),
    paste0("if (length(dat@Dfishery@CAAobs_ymafr) && any(dat@Dfishery@CAAobs_ymafr[, , , ", f, ", ", r, "] > 0, na.rm = TRUE)) {"),
    paste0("  plot_resid_CAA(x, f = ", f, ", r = ", r, ")"),
    "}",
    "```",
    "",
    paste0("```{r caa-hist-resid-f", f, "-r", r, ", fig.cap=\"", cap3, ".\"}"),
    paste0("if (length(dat@Dfishery@CAAobs_ymafr) && any(dat@Dfishery@CAAobs_ymafr[, , , ", f, ", ", r, "] > 0, na.rm = TRUE)) {"),
    paste0("  plot_resid_CAA(x, f = ", f, ", r = ", r, ", do_hist = TRUE)"),
    "}",
    "```",
    ""
  )

  return(rmd)
}

make_rmd_CAL <- function(f, r, fname, rname) {
  cap <- paste0("Catch at length for ", fname, " in ", rname)
  cap2 <- paste0("Pearson residuals for the catch at length for ", fname, " in ", rname)
  cap3 <- paste0("Residual histogram for the catch at length for ", fname, " in ", rname)

  rmd <- c(
    paste0("```{r cal-f", f, "-r", r, ", fig.cap=\"", cap, ".\"}"),
    paste0("if (length(dat@Dfishery@CALobs_ymlfr) && any(dat@Dfishery@CALobs_ymlfr[, , , ", f, ", ", r, "] > 0, na.rm = TRUE)) {"),
    paste0("  plot_CAL(x, f = ", f, ", r = ", r, ")"),
    "}",
    "```",
    "",
    paste0("```{r cal-resid-f", f, "-r", r, ", fig.cap=\"", cap2, ".\"}"),
    paste0("if (length(dat@Dfishery@CALobs_ymlfr) && any(dat@Dfishery@CALobs_ymlfr[, , , ", f, ", ", r, "] > 0, na.rm = TRUE)) {"),
    paste0("  plot_resid_CAL(x, f = ", f, ", r = ", r, ")"),
    "}",
    "```",
    "",
    paste0("```{r cal-hist-resid-f", f, "-r", r, ", fig.cap=\"", cap3, ".\"}"),
    paste0("if (length(dat@Dfishery@CALobs_ymlfr) && any(dat@Dfishery@CALobs_ymlfr[, , , ", f, ", ", r, "] > 0, na.rm = TRUE)) {"),
    paste0("  plot_resid_CAL(x, f = ", f, ", r = ", r, ", do_hist = TRUE)"),
    "}",
    "```",
    ""
  )

  return(rmd)
}



make_rmd_fishery <- function(f, fname, nm = 1, rname, nr = length(rname)) {
  Frate <- ifelse(nm > 1, "season", "year")

  rmd_CAA <- sapply(1:nr, function(i) {
    make_rmd_CAA(f = f, r = i, fname = fname, rname = rname[i])
  }) %>% as.character()

  rmd_CAL <- sapply(1:nr, function(i) {
    make_rmd_CAL(f = f, r = i, fname = fname, rname = rname[i])
  }) %>% as.character()

  rmd <- c(
    paste("###", fname, "{.tabset}"),
    "",
    "#### Fits and residuals",
    "",
    paste0("```{r catch-", f, ", fig.cap=\"Catch of ", fname, ".\"}"),
    paste0("plot_catch(x, f = ", f, ")"),
    "```",
    "",
    paste0("```{r catch-resid-", f, ", fig.cap=\"Catch residuals of ", fname, ".\"}"),
    paste0("plot_resid_Cobs(x, f = ", f, ")"),
    "```",
    "",
    rmd_CAA,
    rmd_CAL,
    "#### Estimates",
    "",
    paste0("```{r sel-", f, ", fig.cap=\"Selectivity.\"}"),
    paste0("plot_self(x, f = ", f, ")"),
    "```",
    "",
    paste0("```{r Ffleet-", f, ", fig.cap=\"Apical fishing mortality of ", fname, " (instantaneous rate, per ", Frate, ").\"}"),
    paste0("plot_Ffleet(x, f = ", f, ")"),
    "```",
    "",
    paste0("```{r catch-region-prop-", f, ", fig.cap=\"Catch by region of ", fname, ".\"}"),
    "if (nr > 1) {",
    paste0("  plot_catch(x, f = ", f, ", by = \"region\", prop = TRUE)"),
    "}",
    "```",
    "",
    paste0("```{r catch-stock-prop-", f, ", fig.cap=\"Catch by stock of ", fname, ".\"}"),
    "if (ns > 1) {",
    paste0("  plot_catch(x, f = ", f, ", by = \"stock\", prop = TRUE)"),
    "}",
    "```",
    "",
    paste0("```{r VB-", f, ", fig.cap=\"Vulnerable biomass.\"}"),
    paste0("plot_V(x, f = ", f, ", by = \"stock\")"),
    "```",
    "",
    paste0("```{r VB-stock-prop-", f, ", fig.cap=\"Proportion vulnerable biomass by stock.\"}"),
    "if (ns > 1) {",
    paste0("  plot_V(x, f = ", f, ", by = \"stock\", prop = TRUE)"),
    "}",
    "```",
    "",
    #paste0("```{r VB-region-", f, ", fig.cap=\"Vulnerable biomass by region.\"}"),
    #"if (nr > 1) {",
    #paste0("  plot_V(x, f = ", f, ", by = \"region\")"),
    #"}",
    #"```",
    #"",
    paste0("```{r VB-region-prop-", f, ", fig.cap=\"Proportion vulnerable biomass by region.\"}"),
    "if (nr > 1) {",
    paste0("  plot_V(x, f = ", f, ", by = \"region\", prop = TRUE)"),
    "}",
    "```",
    ""
  )

  return(rmd)
}



make_rmd_survey <- function(i, iname) {

  rmd <- c(
    ifelse(i == 1, "## Surveys {.tabset}\n", ""),
    paste("###", iname, "{.tabset}"),
    "",
    paste0("```{r isel-", i, ", fig.cap=\"Selectivity for ", iname, ".\"}"),
    paste0("plot_seli(x, i = ", i, ")"),
    "```",
    "",
    paste0("```{r index-", i, ", fig.cap=\"Predicted (red) and observed (black) values for ", iname, ".\"}"),
    paste0("plot_index(x, i = ", i, ")"),
    "```",
    "",
    paste0("```{r IAA-", i, ", fig.cap=\"Age composition from ", iname, ".\"}"),
    paste0("plot_IAA(x, i = ", i, ")"),
    "```",
    "",
    paste0("```{r IAA-resid-", i, ", fig.cap=\"Pearson residuals for the age composition from ", iname, ".\"}"),
    paste0("plot_resid_IAA(x, i = ", i, ")"),
    "```",
    "",
    paste0("```{r IAA-hist-resid-", i, ", fig.cap=\"Residual histogram for the age composition from ", iname, ".\"}"),
    paste0("plot_resid_IAA(x, i = ", i, ", do_hist = TRUE)"),
    "```",
    "",
    paste0("```{r IAL-", i, ", fig.cap=\"Length composition from ", iname, ".\"}"),
    paste0("plot_IAL(x, i = ", i, ")"),
    "```",
    "",
    paste0("```{r IAL-resid-", i, ", fig.cap=\"Pearson residuals for the length composition from ", iname, ".\"}"),
    paste0("plot_resid_IAL(x, i = ", i, ")"),
    "```",
    "",
    paste0("```{r IAL-hist-resid-", i, ", fig.cap=\"Residual histogram for the length composition from ", iname, ".\"}"),
    paste0("plot_resid_IAL(x, i = ", i, ", do_hist = TRUE)"),
    "```",
    ""
  )

  return(rmd)
}



make_rmd_srr <- function(s, sname) {
  rmd <- c(
    paste0("```{r Rdev-", s, ", fig.cap=\"Recruitment deviations of ", sname, ".\"}"),
    paste0("plot_Rdev(x, s = ", s, ")"),
    "```",
    "",
    paste0("```{r SRR-", s, ", fig.cap=\"Stock recruit relationship, and historical stock-recruit pairs, of ", sname, ". The dotted line is the unfished replacement line.\"}"),
    paste0("plot_SRR(x, s = ", s, ")"),
    "```",
    ""
  )

  return(rmd)
}

make_rmd_stock_region <- function(s, sname) {
  rmd <- c(
    ifelse(s == 1, "### By region\n\n", ""),
    paste0("```{r SB-r", s, ", fig.cap=\"Spawning output of ", sname, " by region at the spawning season.\"}"),
    paste0("plot_S(x, by = \"region\", s = ", s, ")"),
    "```",
    "",
    paste0("```{r SB-rp", s, ", fig.cap=\"Proportion spawning output of ", sname, " by region.\"}"),
    "plot_S(x, by = \"region\", s = ", s, ", prop = TRUE)",
    "```",
    ""
  )

  return(rmd)
}


make_rmd_mov <- function(s, y, a, yname, sname, header = TRUE) {
  rmd <- c(
    ifelse(header, "### Movement\n\n", ""),
    paste0("```{r fig.cap=\"Movement of ", sname, " for year ", yname, ", age ", a, " and the corresponding equilibrium distribution.\"}"),
    paste0("plot_mov(x, s = ", s, ", y = ", y, ", a = ", a, ")"),
    "```",
    ""
  )
  return(rmd)
}


