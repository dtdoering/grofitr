#'
#' @import ggplot2

plot_plates <- function(data){
  ggplot(data, aes(x = time, y = OD)) +
    geom_point(size = 0.2) +
    geom_line() +
    facet_grid(`Well Row` ~ `Well Col`, switch = "y") +
    scale_y_continuous(position = "right") +
    theme_bw() +
    theme(panel.grid = element_blank())
}
