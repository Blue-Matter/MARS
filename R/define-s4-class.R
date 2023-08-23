

setOldClass("sdreport")

MSAdata <- setClass(
  "MSAdata",
  slots = c(
    Name = "character",
    Catch = "array",
    Csd = "array",
    Index = "array",
    Isd = "array"
  )
)

MSAassess <- setClass(
  "MSAassess",
  slots = c(
    Name = "character",
    obj = "list",
    opt = "list",
    SD = "sdreport",
    report = "list",
    Misc = "list"
  )
)

setMethod("initialize", "MSAassess", function(.Object, ...) {
  dots <- list(...)
  if (length(dots)) {
    for(i in names(dots)) slot(.Object, i) <- dots[[i]]
  }
  attr(.Object, "version") <- paste("MSA", packageVersion("MSA"))
  attr(.Object, "date") <- date()
  attr(.Object, "R.version") <- getRversion()

  return(.Object)
})
