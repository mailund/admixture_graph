#' Make a data frame an f_4 statistics object.
#' 
#' This is mostly just a convinience function to set the class of a data frame such that 
#' we plot data as error bars in a meaningful way for statistics for admixture graphs.
#' 
#' @param x  Data frame with observed \eqn{D} (\eqn{f_4}) statistics.
#' 
#' @return Something about classes.
#' 
#' @export
f4stats <- function(x) {
  class(x) <- c("f4stats", class(x))
  x
}

#' Plot the fit of a graph to data.
#' 
#' @param x      Data frame with observed \eqn{D} (\eqn{f_4}) statistics
#' @param sigma  How many sigmas the error bars should be wide.
#' @param ...    Additional parameters.
#' 
#' @export
plot.f4stats <- function(x, sigma = 6, ...) {
  # I know this is not the 'dplyr' way of doing it, but package check doesn't like
  # non standard evaluation, so this is what is needed.
  D <- x$D
  x$stderr <- with(x, D / Z.value)
  x$error_bar_start <- with(x, D - sigma/2*stderr)
  x$error_bar_end   <- with(x, D + sigma/2*stderr)
  x$test <- with(x, paste("D(",W,",",X,";",Y,",",Z,")"))
  x$test <- factor(x$test, levels=dplyr::arrange(x, D)$test)
  
  ggplot2::ggplot(x) +
    ggplot2::geom_hline(yintercept = 0, linetype = 'dashed', color = 'gray') +
    ggplot2::geom_errorbar(ggplot2::aes_string(x = 'test', ymin = 'error_bar_start', ymax = 'error_bar_end'), color='black') +
    ggplot2::geom_point(ggplot2::aes_string(x = 'test', y = 'D'), color='black') +
    ggplot2::xlab('') + ggplot2::ylab('') + 
    ggplot2::coord_flip() +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position="none") +
    ggplot2::theme(axis.line = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank()) +
    ggplot2::theme(panel.grid.major = ggplot2::element_line(color = 'white'),
                   panel.grid.major.x = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_line(color = '#eeeeee'))
}