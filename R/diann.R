
#' Read Parquet Files from Multiple Directories
#'
#' Given a list of directories, finds the \code{report.parquet} file in each,
#' reads it using \code{arrow::read_parquet}, and row-binds them into a single
#' data frame.
#'
#' @param dirlist A character vector of directory paths to search for
#'   \code{report.parquet} files.
#'
#' @return A data frame containing the row-bound contents of all discovered
#'   \code{report.parquet} files.
#'
#' @examples
#' \dontrun{
#' dirs <- c("path/to/run1", "path/to/run2")
#' report <- read_parquet_pathlist(dirs)
#' }
#'
#' @importFrom arrow read_parquet
#' @importFrom dplyr bind_rows
#' @export
read_parquet_pathlist <- function(dirlist){
  
  # Given a list of directories, finds the report.parquet in each, reads it
  # and joins them all together.
  
  paths <- lapply(dirlist, list.files, pattern="^report\\.parquet$", full.names=T)
  
  n_matches <- paths |> length()
  
  s <- sprintf("Found %i match(es)", n_matches)
  print(s)
  
  out <- lapply(paths, arrow::read_parquet)
  out.bound <- bind_rows(out)
  return(out.bound)
  
}

read_parquet_sitereport_pathlist <- function(dirlist){
  
  # Given a list of directories, finds the report.sitre_report.parquet in each, reads it
  # and joins them all together.
  
  paths <- lapply(dirlist, list.files, pattern="^report\\.site_report\\.parquet$", full.names=T)
  
  n_matches <- paths |> length()
  
  s <- sprintf("Found %i match(es)", n_matches)
  print(s)
  
  out <- lapply(paths, arrow::read_parquet)
  out.bound <- bind_rows(out)
  return(out.bound)
  
}



#' Extract sample ID and technical ID from a parquet report
#' @description 
#' Apply a regex pattern containing two groups to the 'Run' column of a .parquet report
#' to extract a sample id and a replicate id.
#'
#' @param t A tibble with a DIA-NN parquet report
#' @param patter A regex pattern with exactly two groups, the first being the sample id and the second being the replicate id.
#' @returns An altered tibble which contains two new columns, Sample.ID and Replicate.ID
#' @examples
#' tibble %>% extract_sample_ids(., pattern="P167_R2011_X3961_[A-Z][0-9]+_S([A-Z])([0-9])")
#' 
#' @export
extract_sample_ids <- function(t, pattern){
  
  r <- t %>%
    mutate(
      Sample.ID = str_extract(Run, pattern, group = 1),
      Replicate.ID = str_extract(Run, pattern, group = 2)
    )

  sample.id.missing <- any(is.na(report$Sample.ID))
  replicate.id.missing <- any(is.na(report$Replicate.ID))

  if(sample.id.missing || replicate.id.missing){
    print("At least one sample id or replicate id could not be extracted. Double check your data if this is not expected.")
  }
  
  return(r)
  
}

#'Summarize parquet report to protein groups and their intensities
#'@param t A tibble with a DIA-NN output parquet report
#'@param groups A character vector of columns to group the proteins by
#'@param log2.transform Return protein intensities (i) as log2(i+1)
#'@examples
#'tibble %>% summarize_to_pogs(groups=c("Sample.ID", "Replicate.ID"))
#'@returns A summarized report with proteins and their intensities grouped by the 'groups' parameter. 
#'@export

summarise_to_pgs <- function(t, groups, log2.transform=TRUE){
  
  all.groups <- c(groups, "Protein.Names", "Protein.Group")
  
  r <- t %>%
    group_by(!!!syms(all.groups))%>%
    summarise(PG.MaxLFQ = first(PG.MaxLFQ), .groups = "drop")
  
  if(log2.transform){
    r <- r %>%
      mutate(PG.MaxLFQ = log2(PG.MaxLFQ+1))
  }
  
  return(r)
  
}

#'Summarize parquet report to precursors and their intensities
#'@param t A tibble with a DIA-NN output parquet report
#'@param groups A character vector of columns to group the proteins by
#'@param log2.transform Return protein intensities (i) as log2(i+1)
#'@examples
#'tibble %>% summarize_to_pogs(groups=c("Sample.ID", "Replicate.ID"))
#'@returns A summarized report with proteins and their intensities grouped by the 'groups' parameter. 
#'@export

summarise_to_precursors <- function(t, groups, log2.transform=TRUE){
  
  all.groups <- c(groups, "Protein.Group", "Protein.Names", "Precursor.Id", "Stripped.Sequence", "Modified.Sequence")
  
  r <- t %>%
    group_by(!!!syms(all.groups))%>%
    summarise(Precursor.Normalised = first(Precursor.Normalised), .groups = "drop")
  
  if(log2.transform){
    r <- r %>%
      mutate(Precursor.Normalised = log2(Precursor.Normalised+1))
  }
  
  return(r)
  
}

#' Summarize a DIA-NN parquet report to peptide level using MaxLFQ algorithm
#' 
#' @param t A tibble with a DIA-NN parquet report
#' @param sample.id Name of column containing a sample id
#' @param replicate.id Name of column containing a replicate id
#' @export
summarise_to_peptides <- function(t, sample.id, replicate.id){
  
  temp_tibble <- t %>%
    mutate(temp_id = paste0(!!sym(sample.id), !!sym(replicate.id)))
  
  r <- temp_tibble %>%
    maxlfq(sample.header="temp_id", group.header="Stripped.Sequence")%>%
    data.frame()%>%
    rownames_to_column("Stripped.Sequence")%>%
    pivot_longer(-Stripped.Sequence,names_to="Temporary.Id", values_to = "Peptide.MaxLFQ")
  
  return(r)
  
  
}



#' Add a column 'Missed.Cleavages' with integer number of missed cleavages
#'
#' The calculation is based on internal R/K not suffixed by P. The calculation defaults to use the 'Stripped.Sequence' column of the .parquet report passed
#' @param parquet A DIA-NN outputted .parquet report
#' @return A parquet report with an added column 'Missed.Cleavages'
#' @examples
#' \dontrun{p <- add_missed_cleavage_column(parquet)}
#' @export
add_missed_cleavage_column <- function(parquet){

  p <- parquet %>%
    mutate(Missed.Cleavages = count_missed_cleavage_string(Stripped.Sequence))

  return(p)

}


