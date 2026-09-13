#'@export
internal_trypsin_site_regex <- "(?<=[KR])(?!P).(?<!^)"

#' Count number of missed clevages in string of amino acids/peptide sequence
#' @importFrom stringr str_extract
#' @description Counts the number of internal R/K not suffixed by proline in an AA-string
#' @param seq A string representation of amino acid sequence.
#' @returns An integer
#' @examples count_missed_cleavage_string("AAAALLLRE")
#' @export
count_missed_cleavage_string <- function(seq){
  internal_trypsin_site_regex <- "(?<=[KR])(?!P).(?<!^)"
  count <- stringr::str_count(string = seq, pattern = internal_trypsin_site_regex)
  return(count)
}