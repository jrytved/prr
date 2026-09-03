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
  
  # Given a list of directories, finds the report.parquet in each, reads it
  # and joins them all together.
  
  paths <- lapply(dirlist, list.files, pattern="^report\\.site_report\\.parquet$", full.names=T)
  
  n_matches <- paths |> length()
  
  s <- sprintf("Found %i match(es)", n_matches)
  print(s)
  
  out <- lapply(paths, arrow::read_parquet)
  out.bound <- bind_rows(out)
  return(out.bound)
  
}



#' Extract Run ID from a DIA-NN Report
#'
#' Applies a regular expression to the \code{Run} column of a DIA-NN parquet
#' report, extracts the specified capture group, and stores the result in a new
#' \code{uid} column.
#'
#' @param report A data frame containing at minimum a \code{Run} column, as
#'   produced by \code{read_parquet_pathlist}.
#' @param re A regular expression string containing at least one capture group.
#' @param g An integer indicating which capture group to extract.
#'
#' @return The input data frame with an additional \code{uid} column containing
#'   the extracted values.
#'
#' @examples
#' \dontrun{
#' report <- extract_run_id(report, re = "sample_(\\d+)_", g = 1)
#' }
#'
#' @importFrom dplyr mutate
#' @importFrom stringr str_extract
#' @export
extract_run_id <- function(report, re, g){
  
  # Given a DIA-NN parquet report and a re-pattern, maps the pattern over the Run
  # column, extracts the given group and adds it to a new 'uid'-column.
  report%>%mutate(uid = str_extract(Run, pattern=re, group=g))
  
}



#' Plot Identification Counts Across Conditions
#'
#' Joins a DIA-NN report with sample metadata, computes per-replicate
#' identification counts at the protein group, precursor, peptide, and
#' peptidoform levels, and returns a faceted bar chart with per-condition means,
#' individual replicate jitter points, and CV annotations.
#'
#' @param df A data frame containing DIA-NN output columns including
#'   \code{Protein.Group}, \code{Precursor.Id}, \code{Stripped.Sequence},
#'   \code{Modified.Sequence}, \code{Lib.PG.Q.Value}, \code{Lib.Q.Value}, and
#'   \code{Lib.Peptidoform.Q.Value}.
#' @param meta A data frame of sample metadata to join onto \code{df}. Must
#'   contain columns matching \code{condition_col} and \code{replicate_col}.
#' @param condition_col A string specifying the column name in \code{meta} that
#'   identifies the experimental condition. Defaults to \code{"Condition"}.
#' @param replicate_col A string specifying the column name in \code{meta} that
#'   identifies the replicate. Defaults to \code{"Replicate"}.
#'
#' @return A \code{ggplot} object with one facet per identification level
#'   (\code{N.Protein.Groups}, \code{N.Precursors}, \code{N.Peptide},
#'   \code{N.Peptidoform}), bars representing condition means, jittered points
#'   for individual replicates, and text labels showing the coefficient of
#'   variation.
#'
#' @examples
#' \dontrun{
#' p <- plot_identifications(report, meta, condition_col = "Condition", replicate_col = "Replicate")
#' print(p)
#' }
#'
#' @importFrom dplyr left_join group_by summarise n_distinct across all_of
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_bar facet_wrap scale_fill_brewer
#'   theme_minimal theme element_text geom_jitter geom_text
#' @importFrom RColorBrewer brewer.pal
#' @export
plot_identifications <- function(df, meta, condition_col = "Condition", replicate_col = "Replicate") {
  
  df.joined <- df %>% left_join(meta)
  
  summarise_measures <- function(df){
    df %>% summarise(
      N.Protein.Groups = n_distinct(Protein.Group[Lib.PG.Q.Value <= 0.01]),
      N.Precursors = n_distinct(Precursor.Id[Lib.Q.Value <= 0.01]),
      N.Peptide = n_distinct(Stripped.Sequence[Lib.Peptidoform.Q.Value <= 0.01]),
      N.Peptidoform = n_distinct(Modified.Sequence[Lib.Peptidoform.Q.Value <= 0.01]),
      .groups = "drop"
    )
  }
  
  condition_replicate_summarised <- df.joined %>% 
    group_by(across(all_of(c(condition_col, replicate_col)))) %>%
    summarise_measures() %>%
    pivot_longer(-all_of(c(condition_col, replicate_col)), 
                 names_to = "Level", values_to = "Value")
  
  condition_summarised <- condition_replicate_summarised %>%
    group_by(across(all_of(c(condition_col, "Level")))) %>%
    summarise(Mean = mean(Value), .groups = "drop")
  
  condition_var_summarised <- condition_replicate_summarised %>%
    group_by(across(all_of(c(condition_col, "Level")))) %>%
    summarise(
      CV = sd(Value) / mean(Value),
      Mean = mean(Value),
      .groups = "drop"
    )
  
  p <- ggplot(condition_summarised, aes(x = .data[[condition_col]], y = Mean)) +
    geom_bar(stat = "identity", color = "black", aes(fill = Level)) +
    facet_wrap(~Level, scales = "free_y") +
    scale_fill_brewer(palette = "GnBu") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = -90)) +
    geom_jitter(condition_replicate_summarised, 
                mapping = aes(x = .data[[condition_col]], y = Value), 
                width = 0.1) +
    geom_text(condition_var_summarised,
              mapping = aes(x = .data[[condition_col]], y = Mean + 0.1 * Mean, 
                           label = sprintf("CV=%.1f%%", CV * 100)))
  
  return(p)
}


#' Summarise Identification Counts Across Conditions
#'
#' Joins a DIA-NN report with sample metadata and computes identification counts
#' at the protein group, precursor, peptide, and peptidoform levels, returning
#' either replicate-level or condition-level summaries.
#'
#' @param df A data frame containing DIA-NN output columns including
#'   \code{Protein.Group}, \code{Precursor.Id}, \code{Stripped.Sequence},
#'   \code{Modified.Sequence}, \code{Lib.PG.Q.Value}, \code{Lib.Q.Value}, and
#'   \code{Lib.Peptidoform.Q.Value}.
#' @param meta A data frame of sample metadata to join onto \code{df}. Must
#'   contain columns matching \code{condition_col} and \code{replicate_col}.
#' @param condition_col A string specifying the column name in \code{meta} that
#'   identifies the experimental condition. Defaults to \code{"Condition"}.
#' @param replicate_col A string specifying the column name in \code{meta} that
#'   identifies the replicate. Defaults to \code{"Replicate"}.
#' @param get_replicate_level Logical. If \code{TRUE} (default), returns a
#'   long-format data frame at the replicate level. If \code{FALSE}, returns
#'   condition-level means.
#'
#' @return A long-format data frame with columns for condition, replicate (if
#'   \code{get_replicate_level = TRUE}), \code{Level} (the identification type),
#'   and either \code{Value} (replicate-level) or \code{Mean}
#'   (condition-level).
#'
#' @examples
#' \dontrun{
#' rep_summary <- summarise_identifications(report, meta)
#' cond_summary <- summarise_identifications(report, meta, get_replicate_level = FALSE)
#' }
#'
#' @importFrom dplyr left_join group_by summarise n_distinct across all_of
#' @importFrom tidyr pivot_longer
#' @export
summarise_identifications <- function(df, meta, condition_col = "Condition", replicate_col = "Replicate", get_replicate_level = TRUE) {
  
  df.joined <- df %>% left_join(meta)
  
  summarise_measures <- function(df){
    df %>% summarise(
      N.Protein.Groups = n_distinct(Protein.Group[Lib.PG.Q.Value <= 0.01]),
      N.Precursors = n_distinct(Precursor.Id[Lib.Q.Value <= 0.01]),
      N.Peptide = n_distinct(Stripped.Sequence[Lib.Peptidoform.Q.Value <= 0.01]),
      N.Peptidoform = n_distinct(Modified.Sequence[Lib.Peptidoform.Q.Value <= 0.01]),
      .groups = "drop"
    )
  }
  
  condition_replicate_summarised <- df.joined %>% 
    group_by(across(all_of(c(condition_col, replicate_col)))) %>%
    summarise_measures() %>%
    pivot_longer(-all_of(c(condition_col, replicate_col)), 
                 names_to = "Level", values_to = "Value")
  
  condition_summarised <- condition_replicate_summarised %>%
    group_by(across(all_of(c(condition_col, "Level")))) %>%
    summarise(Mean = mean(Value), .groups = "drop")
  
  if(get_replicate_level){
    return(condition_replicate_summarised)
  } else {
    return(condition_summarised)
  }
  
}


#' Compute Precursor-Level Coefficient of Variation
#'
#' Joins a DIA-NN report with sample metadata, computes the mean, standard
#' deviation, and coefficient of variation (CV) of \code{Precursor.Normalised}
#' intensities per condition and precursor, and flags precursors whose CV falls
#' at or below a given threshold.
#'
#' @param data A data frame containing at minimum \code{Precursor.Id} and
#'   \code{Precursor.Normalised} columns from a DIA-NN report.
#' @param meta A data frame of sample metadata to join onto \code{data}. Must
#'   contain a column matching \code{condition_col}.
#' @param condition_col A string specifying the column name in \code{meta} that
#'   identifies the experimental condition.
#' @param replicate_col A string specifying the column name in \code{meta} that
#'   identifies the replicate. Passed for consistency but not used directly in
#'   grouping.
#' @param n_min Integer. Minimum number of observations required for a
#'   precursor to be retained. Defaults to \code{3}.
#' @param threshold Numeric. CV threshold below which a precursor is considered
#'   low-variability. Defaults to \code{0.2} (20\%).
#'
#' @return A data frame with one row per condition–precursor combination
#'   (filtered to \code{N >= n_min}), containing columns \code{N},
#'   \code{mean}, \code{sd}, \code{CV}, and \code{LT} (logical; \code{TRUE}
#'   when \code{CV <= threshold}).
#'
#' @examples
#' \dontrun{
#' cv <- compute_precursor_cv(report, meta, condition_col = "Condition",
#'                            replicate_col = "Replicate", n_min = 3, threshold = 0.2)
#' }
#'
#' @importFrom dplyr left_join group_by summarise mutate filter across all_of n
#' @export
compute_precursor_cv <- function(data, meta, condition_col, replicate_col, n_min=3, threshold = 0.2){
  
  cv.data <- data %>%
  left_join(meta)%>%
  group_by(across(all_of(c(condition_col, "Precursor.Id"))))%>%
  summarise(
    N = n(),
    mean = mean(Precursor.Normalised+0.1),
    sd = sd(Precursor.Normalised+0.1),
    .groups="drop"
  ) %>%
  mutate(CV = sd/mean)%>%
  filter(N >= n_min)%>%
  mutate(LT = CV <= threshold)
    
  return(cv.data)
  
}


#' Compute Protein Group-Level Coefficient of Variation
#'
#' Joins a DIA-NN report with sample metadata, computes the mean, standard
#' deviation, and coefficient of variation (CV) of \code{PG.MaxLFQ} intensities
#' per condition and protein group, and flags protein groups whose CV falls at
#' or below a given threshold.
#'
#' @param data A data frame containing at minimum \code{Protein.Group} and
#'   \code{PG.MaxLFQ} columns from a DIA-NN report.
#' @param meta A data frame of sample metadata to join onto \code{data}. Must
#'   contain a column matching \code{condition_col}.
#' @param condition_col A string specifying the column name in \code{meta} that
#'   identifies the experimental condition.
#' @param replicate_col A string specifying the column name in \code{meta} that
#'   identifies the replicate. Passed for consistency but not used directly in
#'   grouping.
#' @param n_min Integer. Minimum number of observations required for a protein
#'   group to be retained. Defaults to \code{3}.
#' @param threshold Numeric. CV threshold below which a protein group is
#'   considered low-variability. Defaults to \code{0.2} (20\%).
#'
#' @return A data frame with one row per condition–protein group combination
#'   (filtered to \code{N >= n_min}), containing columns \code{N},
#'   \code{mean}, \code{sd}, \code{CV}, and \code{LT} (logical; \code{TRUE}
#'   when \code{CV <= threshold}).
#'
#' @examples
#' \dontrun{
#' cv <- compute_pg_cv(report, meta, condition_col = "Condition",
#'                     replicate_col = "Replicate", n_min = 3, threshold = 0.2)
#' }
#'
#' @importFrom dplyr left_join group_by summarise mutate filter across all_of n
#' @export
compute_pg_cv <- function(data, meta, condition_col, replicate_col, n_min=3, threshold = 0.2){
  
  cv.data <- data %>%
  left_join(meta)%>%
  group_by(across(all_of(c(condition_col, "Protein.Group"))))%>%
  summarise(
    N = n(),
    mean = mean(PG.MaxLFQ+0.1),
    sd = sd(PG.MaxLFQ+0.1),
    .groups="drop"
  ) %>%
  mutate(CV = sd/mean)%>%
  filter(N >= n_min)%>%
  mutate(LT = CV <= threshold)
    
  return(cv.data)
  
}

#' @export
remove_crap <- function(df){
  df %>%
    filter(!grepl("cRAP", Protein.Ids))
}

#' @export
force_proteotypic <- function(df){
  df %>%
    filter(Proteotypic == T)
}

#' @export
is_phosphorylated <- function(v){
  return(grepl("UniMod:21", v))
}