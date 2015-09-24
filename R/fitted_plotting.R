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



unpack_environment <- function(parameters, x) {
  n_edges <- length(parameters$edges)
  n_admix_prop <- length(parameters$admix_prop)
  edges <- x[1:n_edges]
  admix_prop <- x[ (n_edges + 1) : (n_edges + n_admix_prop)]
  graph_environment(parameters, edges, admix_prop)
}

make_predict_function <- function(data, graph, parameters = extract_graph_parameters(graph)) {
  force(data)
  force(graph)
  force(parameters)
  
  goal <- data$D
  expressions <- Map(function(W,X,Y,Z) sf4(graph, W, X, Y, Z), data$W, data$X, data$Y, data$Z)
  
  function(x) {
    env <- unpack_environment(parameters, x)
    predictions <- unlist(Map(function(expression) eval(expression, env), expressions), use.names = FALSE)
    predictions
  }
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
#' @param show_fit    Should the function plot the number of tests where the graph fits
#'                    the data instead of the sum of squared errors?
#' @param sigma       If the function plots the fit, it is considered a hit if the difference
#'                    between a prediction and the observed statistics is no more than sigma/2.
#' @param ...         Additional parameters passed to the plotting function filled.contour
#'   
#' @return  The matrix of values computed and plotted.
#'   
#' @export
contour_plot <- function(object, X, Y, resolution = 10, show_fit = FALSE, sigma = 6, ...) {
  
  fitted_parameters <- coef(object)
  if (! X %in% names(fitted_parameters)) {
    stop(paste("'", X,"' is not a parameter of the fitted graph.", sep = ""))
  }
  if (! Y %in% names(fitted_parameters)) {
    stop(paste("'", Y,"' is not a parameter of the fitted graph.", sep = ""))
  }
  
  fitted_parameters <- coef(object)
  best_x <- fitted_parameters[X]
  best_y <- fitted_parameters[Y]
  
  x <- seq(0, resolution)
  y <- seq(0, resolution)
  z <- matrix(0, resolution, resolution)
  
  predict_f4 <- make_predict_function(object$data, object$graph)
  evaluate_point_SSE <- function(point) {
    target <- object$data$D
    prediction <- predict_f4(point)
    sum((target-prediction)**2)
  }
  evaluate_point_fits <- function(point) {
    target <- object$data$D
    prediction <- predict_f4(point)
    sum(abs(target-prediction)/object$data$Z.value > sigma/2)
  }
  
  if (show_fit) {
    evaluate_point <- evaluate_point_fits
  } else {
    evaluate_point <- evaluate_point_SSE
  }
  
  point <- coef(object)# object$best_fit
  for (i in seq (1, resolution)) {
    for (j in seq(1, resolution)) {
      point[X] <- i/resolution
      point[Y] <- j/resolution
      z[i, j] <- evaluate_point(point)
    }  
  }
  x <- 1:nrow(z)/resolution
  y <- 1:ncol(z)/resolution
  
  image(x, y, z, xlab = X, ylab = Y, col = rev(heat.colors(12)), ...)
  contour(x, y, z, add = TRUE, ...)
  points(best_x, best_y, pch=3)
  invisible(z)
}


