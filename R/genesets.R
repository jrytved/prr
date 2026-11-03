#' Read a .gmt file containing genesets
#' @param path A path to the .gmt file
#' @return A named list containing the gene sets
#' @export
read_gmt <- function(path) {
  lines  <- readLines(path, warn = FALSE)
  fields <- strsplit(lines, "\t", fixed = TRUE)
  fields <- fields[lengths(fields) >= 3]
  sets   <- lapply(fields, function(x) { g <- x[-(1:2)]; g[g != ""] })
  setNames(sets, vapply(fields, `[`, "", 1))
}