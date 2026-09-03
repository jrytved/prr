
#' Compute Jaccard Index
#' Given a named list of sets, returns a matrix with the Jaccard indices of each pairwise intersection / union.
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

#' Get unique sets of items in a target column, grouped by a grouping column.
#' 
#' @export
get_unique_column_set <- function(df, groupcol = "uid", targetcol ="Stripped.Sequence"){
  
  
  ids <- df %>%
  pull(groupcol)%>%
  unique()

  n.ids <- length((ids))
  sets <- vector(mode="list", length=length(ids))
  names(sets) <- ids
  
  for(i in seq_along(ids)){
    
    id <- ids[[i]]
    
    unique.value <- df %>%
      filter(.data[[groupcol]] == id)%>%
      pull(targetcol)%>%
      unique()
    
    sets[[id]] <- unique.value
    
  }
  
  return(sets)
  
}
