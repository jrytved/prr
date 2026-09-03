#' Flip x labels of any plot with 90 degrees
#' @export
flip_xlab <- function(p){
  theme(axis.text.x=element_text(angle=90))
}

#' Get the sweet Johan style on your current ggplot
#' @export
style_gg_johan <- function(){
  theme_classic()+
  theme(
    axis.text.x=element_text(size=10),
    axis.text.y=element_text(size=10),
    axis.title.x=element_text(size=14)
  )
}