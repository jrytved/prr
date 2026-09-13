#' Make a gt table presentable
#' @description Style the gt-table in a latex booktabs-ish style.
#' @param gt_tbl a gt table
#' @param font Typeface
#' @param title Plot title
#' @param subtitle Plot subtitle
#' @param footnote Plot footnote
#' @param source_note Plot sourcenote
#' @param  add_horizontal_rules Separate rows by horizontal rules
#' @export 
make_presentable_table <- function(gt_tbl, 
                                   font = "Computer Modern",
                                   title = NULL,
                                   subtitle = NULL,
                                   footnote = NULL,
                                   source_note = NULL,
                                   add_horizontal_rules = FALSE) {
  
  tbl <- gt_tbl %>%
    
    gt::tab_options(
      table.font.names = font,
      table.font.size = gt::px(13),
      
      column_labels.font.size = gt::px(13),
      heading.title.font.size = gt::px(15),
      heading.title.font.weight = "bold",
      heading.subtitle.font.size = gt::px(13),
      
      table.border.top.style = "solid",
      table.border.top.width = gt::px(2),
      table.border.top.color = "black",
      
      table.border.bottom.style = "solid",
      table.border.bottom.width = gt::px(2),
      table.border.bottom.color = "black",
      
      column_labels.border.top.style = "solid",
      column_labels.border.top.width = gt::px(1),
      column_labels.border.top.color = "black",
      
      column_labels.border.bottom.style = "solid",
      column_labels.border.bottom.width = gt::px(1),
      column_labels.border.bottom.color = "black",
      
      table_body.hlines.style = if (add_horizontal_rules) "solid" else "none",
      table_body.hlines.width = if (add_horizontal_rules) gt::px(0.5) else gt::px(0),
      table_body.hlines.color = if (add_horizontal_rules) "grey80" else "black",
      
      table_body.vlines.style = "none",
      column_labels.vlines.style = "none",
      
      row_group.border.top.style = "solid",
      row_group.border.top.width = gt::px(1),
      row_group.border.top.color = "black",
      row_group.border.bottom.style = "solid",
      row_group.border.bottom.width = gt::px(1),
      row_group.border.bottom.color = "black",
      
      table.background.color = "white",
      
      data_row.padding = gt::px(4),
      column_labels.padding = gt::px(4),
      
      footnotes.font.size = gt::px(11),
      source_notes.font.size = gt::px(11),
      footnotes.marks = "standard" 
    ) %>%
    
    gt::cols_align(align = "center") %>%
    
    gt::cols_align(align = "left", columns = 1)

  if (!is.null(title)) {
    tbl <- tbl %>% gt::tab_header(title = title, subtitle = subtitle)
  }
  
  if (!is.null(footnote)) {
    tbl <- tbl %>% gt::tab_footnote(footnote = footnote)
  }
  
  if (!is.null(source_note)) {
    tbl <- tbl %>% gt::tab_source_note(source_note = source_note)
  }
  
  tbl
}