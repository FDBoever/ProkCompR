#' Adds custom theme to ggplot
#'
#' @param aspect.ratio sets aspect.ratio in ggplot2::theme()
#'
#' @keywords fdb_style
#' @export
#' @examples
#' ggplot(data=mtcars, aes(cyl,mpg)) +
#' geom_point() +
#' geom_smooth(method='lm')+
#' ggtitle('mtcars') +
#' fdb_style()

fdb_style <- function(aspect.ratio=1) {
  #font <- "Helvetica"
  theme_classic() + ggplot2::theme(
    #Set aspect ratio, this is changable
    aspect.ratio = aspect.ratio,

    #Format title
    plot.title = ggplot2::element_text(hjust = 0.5),

    #Axis format
    axis.title = ggplot2::element_text(size = 11, colour = '#000000'),
    axis.text = ggplot2::element_text(size = 10, colour = '#000000'),

    #Legend format
    legend.justification = c(1, 1),
    legend.key.width = unit(0.25, 'cm'),
    legend.key.height = unit(0.55, 'cm'),
    legend.text = ggplot2::element_text(size = 10),
    legend.title = ggplot2::element_text(size = 11),
    #legend.position='none',

    #Blank background
    #This sets the panel background as blank, removing the standard grey ggplot background colour from the plot
    panel.background = ggplot2::element_blank(),

    #Strip background (#This sets the panel background for facet-wrapped plots to white, removing the standard grey ggplot background colour and sets the title size of the facet-wrap title to font size 22)
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(size  = 11,  hjust = 0)
   # theme(strip.background = element_blank())

    #Grid lines
    # panel.grid.minor = ggplot2::element_blank(),
    # panel.grid.major.y = ggplot2::element_line(color="#cbcbcb"),
    # panel.grid.major.x = ggplot2::element_blank(),
    )
}
