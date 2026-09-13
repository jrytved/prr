#' Flip x labels of any plot with 90 degrees
#' @export
flip_xlab <- function(p){
  theme(axis.text.x=element_text(angle=90))
}

#' Make a figure presentable
#' @export
make_presentable <- function(base_size = 12,margin=15,font = "Arial", legend.position="top") {
  
  theme_classic(base_size = base_size, base_family = font) +
    theme(
      
      axis.text.x = element_text(size = base_size, color = "black"),
      axis.text.y = element_text(size = base_size, color = "black"),
      axis.title.x = element_text(size = base_size+2, color = "black", face = "plain"),
      axis.title.y = element_text(size = base_size+2, color = "black", face = "plain"),
      
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      axis.ticks.length = unit(0.15, "cm"),
      
      legend.title = element_blank(),
      legend.text = element_text(size = base_size, color = "black"),
      legend.background = element_blank(),
      legend.key = element_blank(),
      
      plot.title = element_text(size = base_size+4, face = "bold", hjust = 0.5),
      
      panel.grid = element_blank(),
      panel.background = element_blank(),
      plot.background = element_blank(),
      plot.margin = margin(margin,margin,margin,margin),
      legend.position = legend.position
    )
}