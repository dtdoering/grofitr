#' Plot tidy plate reader data with ggplot2
#'
#' Uses ggplot2 to plot an entire plate of data in a 12 x 8 grid.
#'
#' @param data Plate reader data in tidy (long) format.
#'
#' @import ggplot2
#'
#' @export
#'

plot_plates <- function(data){
  ggplot(data, aes(x = time, y = OD)) +
    geom_point(size = 0.2) +
    geom_line() +
    facet_grid(`Well Row` ~ `Well Col`, switch = "y") +
    scale_y_continuous(position = "right") +
    theme_bw() +
    theme(panel.grid = element_blank())
}
