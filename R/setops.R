
#' Compute Jaccard Index
#' @description Given a named list of sets, returns a matrix with the Jaccard indices of each pairwise intersection / union.
#' @param named.setlist A named list of sets (proteins, peptides, african frog species...)
#' @returns A matrix of pairwise jaccard indices
#' @export
compute_jaccard_index <- function(named.setlist) {
  
  n <- length(named.setlist)
  nms <- names(named.setlist)
  
  mat <- matrix(0, nrow = n, ncol = n, dimnames = list(nms, nms))
  
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      a <- named.setlist[[i]]
      b <- named.setlist[[j]]
      intersection <- length(intersect(a, b))
      union <- length(union(a, b))
      mat[i, j] <- if (union == 0) 0 else intersection / union
    }
  }
  
  return(mat)
}
