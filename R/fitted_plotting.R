#' Plot the fit of a graph to data.
#' 
#' @param x       Fitted graph object.
#' @param sigma   How many sigmas the error bars should be wide.
#' @param ...     Additional parameters.
#' 
#' @import ggplot2
#' @import dplyr
#' @export
plot.agraph_fit <- function(x, sigma = 6, ...) {
  
  fit <- fitted(x)
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
    (if (all(fit$hit))        scale_color_manual(values = c("green"))
     else if (all (!fit$hit)) scale_color_manual(values = c("red"))
     else                     scale_color_manual(values = c("red", "green"))) +
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

#' A contour plot of the cost function around the best fit with respect to 2
#' admix variables specified by the user.
#' 
#' Returning the cost, an example edge solution of an optimal fit, and linear 
#' relations describing the set of all edge solutions. Operating with the full
#' edge optimisation matrix.
#' 
#' @param object      The fitted object.
#' @param X           An admix variable name (remember quotation marks) or number.
#' @param Y           An admix variable name (remember quotation marks) or number.
#' @param resolution  How densely is the function evaluated.
#' @param ...         Additional parameters passed to the plotting function filled.contour
#'   
#' @return  Just a contour plot with FIERY colours.
#'   
#' @export
contour_plot <- function(object, X, Y, resolution = 10, ...) {
  if (!requireNamespace("grDevices", quietly = TRUE)) {
    stop("This function requires grDevices to be installed.")
  }
  
  fitted_parameters <- coef(object)
  if (! X %in% names(fitted_parameters)) {
    stop(paste(X, "is not a parameter of the fitted graph."))
  }
  
  x <- seq(0, resolution)
  y <- seq(0, resolution)
  z <- matrix(0, resolution, resolution)
  point <- object$best_fit
  for (i in seq (1, resolution)) {
    for (j in seq(1, resolution)) {
      point[X] <- i/resolution
      point[Y] <- j/resolution
      z[i, j] <- cost_function(object$data, object$matrix$column_reduced,
                               object$graph, object$parameters)(point)
    }  
  }
  x <- 1:nrow(z)/resolution
  y <- 1:ncol(z)/resolution
  
  filled.contour(x, y, z, xlab = X, ylab = Y, color.palette = grDevices::heat.colors, ...)
}


