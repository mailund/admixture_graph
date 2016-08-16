#' Plot the fit of a graph to data.
#' 
#' Plot the fit of a graph to data.
#' 
#' @param x          Fitted graph object.
#' @param sigma      How many standard deviations the error bars should be wide.
#' @param grayscale  Should the plot be in black and white?
#' @param ...        Additional parameters.
#' 
#' @export
plot.agraph_fit <- function(x, sigma = 6, grayscale = FALSE, ...) {
  
  fit <- stats::fitted(x)
  D <- fit$D
  fit$stderr <- with(fit, D / Z.value)
  fit$error_bar_start <- with(fit, D - sigma/2*stderr)
  fit$error_bar_end   <- with(fit, D + sigma/2*stderr)
  fit$test <- with(fit, paste("D(",W,",",X,";",Y,",",Z,")"))
  fit$test <- factor(fit$test, levels=dplyr::arrange(fit, D)$test)
  fit$hit <- with(fit, Vectorize(dplyr::between)(graph_f4, error_bar_start, error_bar_end))

  if (grayscale) {
    hitcolor <- "black"
    misscolor <- "black"
  } else {
    hitcolor <- "green"
    misscolor <- "red"
  }
  
  ggplot2::ggplot(fit) +
    ggplot2::geom_hline(yintercept = 0, linetype = 'dashed', color = 'gray') +
    ggplot2::geom_segment(ggplot2::aes_string(x = 'test', xend = 'test', y = 'D', 
                                              yend = 'graph_f4', color = 'hit'), 
                          linetype = 'dashed') +
    ggplot2::geom_errorbar(ggplot2::aes_string(x = 'test', 
                                               ymin = 'error_bar_start', 
                                               ymax = 'error_bar_end'), color='black') +
    ggplot2::geom_point(ggplot2::aes_string(x = 'test', y = 'D'), color='black', shape = 3) +
    ggplot2::geom_point(ggplot2::aes_string(x = 'test', y = 'graph_f4', color = 'hit'), shape = 16) +
    (if (all(fit$hit))        ggplot2::scale_color_manual(values = c(hitcolor))
     else if (all (!fit$hit)) ggplot2::scale_color_manual(values = c(misscolor))
     else                     ggplot2::scale_color_manual(values = c(misscolor, hitcolor))) +
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
    predictions <- unlist(Map(function(expression) eval(expression, env), expressions), 
                          use.names = FALSE)
    predictions
  }
}

#' A contour plot of the cost function.
#'
#' A contour plot of the cost function with respect to two admix variables specified by the user.
#' 
#' @param object      The fitted object.
#' @param X           An admix variable name (remember quotation marks).
#' @param Y           An admix variable name (remember quotation marks).
#' @param resolution  How densely is the function evaluated.
#' @param show_fit    Should the function plot the number of statistics where the graph fits
#'                    the data instead of the \code{\link{cost_function}}?
#' @param sigma       If \code{show_fit} is \code{TRUE} then each statistic is considered fitted
#'                    if the difference between a prediction and the observation statistics is no
#'                    more than \eqn{D*\sigma/(2*Z)}. Notice that even when plotting the number of
#'                    fitted statistics, we have no guarantee that the chosen variables maximize
#'                    this number as the fitting function still optimizes \code{\link{cost_function}}.
#' @param grayscale   Should the figure be plotted in grayscale or in colour?
#' @param ...         Additional parameters passed to the plotting function
#'                    \code{\link{contour}}.
#'   
#' @return The matrix of values computed and plotted.
#'
#' @seealso \code{\link{contour}}
#' @seealso \code{\link{plot_fit_1}}
#'
#' @export
plot_fit_2 <- function(object, X, Y, 
                       resolution = 10, 
                       show_fit = FALSE, 
                       sigma = 6, 
                       grayscale = FALSE, 
                       ...) {
  
  if (show_fit == TRUE) {
    if("Z.value" %in% colnames(data) == TRUE) {
      zvalues <- data$Z.value
    } else {
      cat("No column called \"Z.value\", using default value one for each. \n")
      zvalues <- rep(1, length(data$D))
    }
    tol <- data$D/zvalues
    tol <- sigma*tol/2
  }
  
  data <- object$data
  reduced_matrix <- object$matrix$column_reduced
  graph <- object$graph
  
  epsilon <- 1e-5
  parameters <- object$parameters
  min <- rep(epsilon, length(object$best_fit))
  names(min) <- names(object$best_fit)
  max <- rep(1 - epsilon, length(object$best_fit))
  names(max) <- names(object$best_fit)
  point <- list(min, max)
  
  
  evaluate_point_cost <- function(point) { # Plotting the cost function:
    fast_fit(data, graph, point)$best_error
  }
  evaluate_point_fits <- function(point) { # Plotting the number of fits:
    residuals <- residuals(fit_graph(data, graph, point))
    number <- 0
    for (i in seq(1, length(tol))) {
      if (abs(residuals[i]) < tol[i]) {
        number <- number + 1
      }
    }
    number
  }
  if (show_fit == TRUE) {
    evaluate_point <- evaluate_point_fits
  } else {
    evaluate_point <- evaluate_point_cost
  }
  
  x <- seq(epsilon, 1-epsilon, length.out = resolution)
  y <- seq(epsilon, 1-epsilon, length.out = resolution)
  z <- matrix(NA, resolution, resolution)
  
  for (i in seq_along(x)) {
    for (j in seq_along(y)) {
      point[[1]][X] <- x[i]
      point[[2]][X] <- x[i]
      point[[1]][Y] <- y[j]
      point[[2]][Y] <- y[j]
      z[i, j] <- evaluate_point(point)
    }
  }
  
  if (grayscale) {
    palette <- rev(grDevices::gray(1:10 / 12))
  } else {
    palette <- rev(grDevices::heat.colors(12))
  }
  
  graphics::image(x, y, z, xlab = X, ylab = Y, col = palette, ...)
  graphics::contour(x, y, z, add = TRUE, ...)
  
  fitted_parameters <- object$best_fit # We draw a little plus sign on the best point.
  best_x <- fitted_parameters[X]
  best_y <- fitted_parameters[Y]
  
  graphics::points(best_x, best_y, pch = 3)
  
  invisible(z)
}

#' A plot of the cost function or number of fitted statistics.
#'
#' A plot of the cost function with respect to one admix variable specified by the user.
#' Sorry about the name, all the good ones were taken and the fact that the word "graph" means
#' two different things doesn't help any.
#' 
#' @param object      The fitted object.
#' @param X           An admix variable name (remember quotation marks).
#' @param resolution  How densely is the function evaluated.
#' @param show_fit    Should the function plot the number of statistics where the graph fits
#'                    the data instead of the \code{\link{cost_function}}?
#' @param sigma       If \code{show_fit} is \code{TRUE} then each statistic is considered fitted
#'                    if the difference between a prediction and the observation statistics is no
#'                    more than \eqn{D*\sigma/(2*Z)}. Notice that even when plotting the number of
#'                    fitted statistics, we have no guarantee that the chosen variables maximize
#'                    this number as the fitting function still optimizes 
#'                    \code{\link{cost_function}}.
#' @param ...         Additional parameters.
#'
#' @return Values for optimal cost function for values of \code{X} between zero and one, plotted.
#'
#' @seealso \code{\link{plot_fit_2}}
#'
#' @export
plot_fit_1 <- function(object, X, resolution = 100, show_fit = FALSE, sigma = 6, ...) {

  fitted_parameters <- object$best_fit # We draw a little plus sign on the best point.
  
  data <- object$data
  reduced_matrix <- object$matrix$column_reduced
  graph <- object$graph
  
  epsilon <- 1e-5
  parameters <- object$parameters
  min <- rep(epsilon, length(object$best_fit))
  names(min) <- names(object$best_fit)
  max <- rep(1 - epsilon, length(object$best_fit))
  names(max) <- names(object$best_fit)
  point <- list(min, max)
  
  evaluate_point_cost <- function(point) { # Plotting the cost function:
    fast_fit(data, graph, point)$best_error
  }
  
  if (show_fit == TRUE) {
    if("Z.value" %in% colnames(data) == TRUE) {
      zvalues <- data$Z.value
    } else {
      cat("No column called \"Z.value\", using default value one for each. \n")
      zvalues <- rep(1, length(data$D))
    }
    tol <- data$D/zvalues
    tol <- sigma*tol/2
  }
  evaluate_point_fits <- function(point) { # Plotting the number of fits:
    residuals <- residuals(fit_graph(data, graph, point))
    number <- 0
    for (i in seq(1, length(tol))) {
      if (abs(residuals[i]) < tol[i]) {
        number <- number + 1
      }
    }
    number
  }
  
  if (show_fit == TRUE) {
    ep <- evaluate_point_fits
    ylabel <- "Number of fitted statistics"
  } else {
    ep <- evaluate_point_cost
    ylabel <- "Cost function"
  }
  evaluate_point <- Vectorize(function(x) {
    p <- point
    p[[1]][X] <- x
    p[[2]][X] <- x
    ep(p)
  })
  
  x <- seq(epsilon, 1-epsilon, length.out = resolution)
  y <- evaluate_point(x)
  
  best_x <- fitted_parameters[X]
  best_y <- evaluate_point(best_x)
  
  graphics::plot(x, y, xlab = X, ylab = ylabel, type = "h", ...)
  graphics::points(best_x, best_y, pch = 3)
  invisible(y)
}
