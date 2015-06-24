#' Plot the fit of a graph to data.
#' 
#' @param x       Fitted graph object.
#' @param sigma   How many sigmas the error bars should be wide.
#' @param ...     Additional parameters.
#' 
#' @import dplyr
#' @import ggplot2
#' @export
plot.agraph_fit <- function(x, sigma = 5, ...) {
  
  # I know this is not the 'dplyr' way of doing it, but package check doesn't like
  # non standard evaluation, so this is what is needed.
  fit <- x$fit_data
  fit$stderr <- with(fit, D / Z.value)
  fit$error_bar_start <- with(fit, D - sigma/2*stderr)
  fit$error_bar_end   <- with(fit, D + sigma/2*stderr)
  fit$test <- with(fit, paste("D(",W,",",X,";",Y,",",Z,")"))
  fit$test <- factor(fit$test, levels=arrange(fit, D)$test)
  fit$hit <- with(fit, Vectorize(between)(graph_f4, error_bar_start, error_bar_end))

  ggplot(fit) +
    geom_hline(yintersect = 0, linetype = 'dashed', color = 'gray') +
    geom_segment(aes_string(x = 'test', xend = 'test', y = 'D', yend = 'graph_f4', color = 'hit'), 
                 linetype = 'dashed') +
    geom_errorbar(aes_string(x = 'test', ymin = 'error_bar_start', ymax = 'error_bar_end'), color='black') +
    geom_point(aes_string(x = 'test', y = 'D'), color='black') +
    geom_point(aes_string(x = 'test', y = 'graph_f4', color = 'hit')) +
    scale_color_manual(values = c("red", "green")) +
    xlab('') + ylab('') + 
    coord_flip() +
    theme_classic() +
    theme(legend.position="none") +
    theme(axis.text.y = element_text(size = 4))
    theme(axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks = element_blank()) +
    theme(panel.grid.major = element_line(color = 'white'),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(color = '#eeeeee'))
}
