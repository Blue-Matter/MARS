
make_rmd_fishery <- function(f, fname, nm = 1) {
  Frate <- ifelse(nm > 1, "season", "year")

  rmd <- c(
    paste("###", fname),
    "",
    paste0("```{r sel-", f, ", fig.cap=\"Selectivity.\"}"),
    paste0("plot_self(x, f = ", f, ")"),
    "```",
    "",
    paste0("```{r Ffleet-", f, ", fig.cap=\"Apical fishing mortality of ", fname, " (instantaneous rate, per ", Frate, ").\"}"),
    paste0("plot_Ffleet(x, f = ", f, ")"),
    "```",
    "",
    paste0("```{r VB-", f, ", fig.cap=\"Vulnerable biomass.\"}"),
    paste0("plot_V(x, f = ", f, ", by = \"stock\")"),
    "```",
    "",
    paste0("```{r VB-sprop-", f, ", fig.cap=\"Proportion vulnerable biomass by stock.\"}"),
    "if (ns > 1) {",
    paste0("  plot_V(x, f = ", f, ", by = \"stock\", prop = TRUE)"),
    "}",
    "```",
    "",
    paste0("```{r VB-r-", f, ", fig.cap=\"Vulnerable biomass by region.\"}"),
    "if (nr > 1) {",
    paste0("  plot_V(x, f = ", f, ", by = \"region\")"),
    "}",
    "```",
    "",
    paste0("```{r VB-rprop-", f, ", fig.cap=\"Proportion vulnerable biomass by region.\"}"),
    "if (nr > 1) {",
    paste0("  plot_V(x, f = ", f, ", by = \"region\", prop = TRUE)"),
    "}",
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
    paste0("```{r fig.cap=\"Movement of ", sname, " for year ", yname, ", age ", a, ".\"}"),
    paste0("plot_mov(x, s = ", s, ", y = ", y, ", a = ", a, ")"),
    "```",
    ""
  )
  return(rmd)
}


