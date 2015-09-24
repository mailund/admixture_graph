#' Make a data frame an f4 statistics object.
#' 
#' This is mostly just a convinience function to set the class of a data frame such that 
#' we plot data as error bars in a meaningful way for statistics for admixture graphs.
#' 
#' @param x       Data frame with observed D (f4) statistics
#' 
#' @export
f4stats <- function(x) {
  class(x) <- c("f4stats", class(x))
  x
}

#' Plot the fit of a graph to data.
#' 
#' @param x       Data frame with observed D (f4) statistics
#' @param sigma   How many sigmas the error bars should be wide.
#' @param ...     Additional parameters.
#' 
#' @import ggplot2
#' @import dplyr
#' @export
plot.f4stats <- function(x, sigma = 6, ...) {
  
  # I know this is not the 'dplyr' way of doing it, but package check doesn't like
  # non standard evaluation, so this is what is needed.
  x$stderr <- with(x, D / Z.value)
  x$error_bar_start <- with(x, D - sigma/2*stderr)
  x$error_bar_end   <- with(x, D + sigma/2*stderr)
  x$test <- with(x, paste("D(",W,",",X,";",Y,",",Z,")"))
  x$test <- factor(x$test, levels=arrange(x, D)$test)
  
  ggplot(x) +
    geom_hline(yintersect = 0, linetype = 'dashed', color = 'gray') +
    geom_errorbar(aes_string(x = 'test', ymin = 'error_bar_start', ymax = 'error_bar_end'), color='black') +
    geom_point(aes_string(x = 'test', y = 'D'), color='black') +
    xlab('') + ylab('') + 
    coord_flip() +
    theme_classic() +
    theme(legend.position="none") +
    theme(axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks = element_blank()) +
    theme(panel.grid.major = element_line(color = 'white'),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(color = '#eeeeee'))
}
