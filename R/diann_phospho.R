#' Combine site-level and global parquet reports appropriately
#' @param site_parquet An object containing the site-level parquet report output by dia-nn.
#' @param full_parquet An object containing the full (report.parquet) output by dia-nn.
#' @return A data frame containing the joined site-level report (by precursor lib index and auxillary columns)
#' @export
#'
join_site_report <-function(site_parquet, full_parquet, join_by = c("Sample.ID", "Replicate.ID")){

    join_by = c(join_by, "Precursor.Lib.Index")

    dplyr::left_join(
      site.report, full_parquet,
      by=c("Precursor.Lib.Index", "Sample.ID", "Replicate.ID")
    )

}
