
pack_environment <- function(parameters, env) {
  n_edges <- length(parameters$edges)
  n_admix_prop <- length(parameters$admix_prop)
  x <- rep(NA, n_edges + n_admix_prop)
  for (i in seq_along(parameters$edges)) {
    x[i] <- get(parameters$edges[i], env)
  }
  for (i in seq_along(parameters$admix_prop)) {
    x[i + n_edges] <- get(parameters$admix_prop[i], env)
  }
  x
}

unpack_environment <- function(parameters, x) {
  n_edges <- length(parameters$edges)
  n_admix_prop <- length(parameters$admix_prop)
  edges <- x[1:n_edges]
  admix_prop <- x[ (n_edges + 1) : (n_edges + n_admix_prop)]
  graph_environment(parameters, edges, admix_prop)
}

make_cost_function <- function(data, graph,
                               parameters = extract_graph_parameters(graph)) {
  force(data)
  force(graph)
  force(parameters)
  
  goal <- data$D
  expressions <- Map(function(W,X,Y,Z) sf4(graph, W, X, Y, Z), data$W, data$X, data$Y, data$Z)
  
  function(x) {
    env <- unpack_environment(parameters, x)
    predictions <- unlist(Map(function(expression) eval(expression, env), expressions), use.names = FALSE)
    sum( (goal - predictions) ** 2)
  }
}

#' Fit the graph parameters to a data set.
#' 
#' Tries to minimize the squared distance between statistics in \code{data} and 
#' statistics given by the graph.
#' 
#' The data frame, \code{data}, must contain columns \code{W}, \code{X}, 
#' \code{Y}, and \code{Z}. The function then computes the \eqn{f_4(W,X;Y,Z)}
#' statistics for all rows from these to obtain the prediction made by the
#' graph.
#' 
#' The data frame must also contain a column, \code{D}, containing the 
#' statistics observed in the data. The fitting algorithm attempts to minimize 
#' the distance from this column and the predictions made by the graph.
#' 
#' @param data  The data set.
#' @param graph The admixture graph.
#' @param optimisation_options  Options to the optimisation algorithm.
#'   
#' @return A list containing the best fitted values (in an environment) and the 
#'   data extended with a column containing the graph predictions.
#'   
#' @seealso \code{\link[neldermead]{optimset}}
#'   
#' @export
fit_graph <- function(data, graph, optimisation_options = NULL) {
  if (!requireNamespace("neldermead", quietly = TRUE)) {
    stop("This function requires neldermead to be installed.")
  }

  params <- extract_graph_parameters(graph)
  init_env <- graph_environment(params)
  x0 <- pack_environment(params, init_env)
  cfunc <- make_cost_function(data, graph, params)

  opti <- neldermead::fminbnd(cfunc, x0 = x0, xmin = rep(0, length(x0)),
                              xmax = rep(1, length(x0)),
                              options = optimisation_options)

  best_fit <- neldermead::neldermead.get(opti, "xopt")
  best_fit_env <- unpack_environment(params, best_fit)
  best_fit_data <- add_graph_f4(data, graph, best_fit_env)
  
  structure(list(
      call = sys.call(),
      graph = graph,
      params = params,
      error = with(best_fit_data, sum((D-graph_f4)**2)),
      fit_env = best_fit_env, 
      fit_data = best_fit_data),
    class = "agraph_fit")
}

#' Print function for a fitted graph.
#' 
#' Print summary of the result of a fit.
#' 
#' @param x  The fitted object.
#' @param ...     Additional arguments.
#'  
#' @export
print.agraph_fit <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat("Sum of squared error:", x$error, "\n")
  invisible(x)
}

fitted_parameters <- function(object) {
  edges <- unlist(Map(function(param) get(param, object$fit_env), object$param$edges))
  admixture_proportions <- unlist(Map(function(param) get(param, object$fit_env), object$param$admix_prop))
  list(edges = edges, admixture_proportions = admixture_proportions)
}

#' Get fitted parameters for a fitted graph.
#' 
#' Extract the graph parameters for a graph fitted to data.
#' 
#' @param object  The fitted object.
#' @param ...     Additional arguments.
#' 
#' @export
coef.agraph_fit <- function(object, ...) {
  parameters <- fitted_parameters(object)
  c(parameters$edges, parameters$admixture_proportions)
}

#' Print function for a fitted graph.
#' 
#' Print summary of the result of a fit.
#' 
#' @param object  The fitted object.
#' @param ...     Additional arguments.
#' 
#' @export
summary.agraph_fit <- function(object, ...) {
  result <- fitted_parameters(object)
  print(result)
  invisible(result)
}


#' Extract fitted data for a fitted graph.
#' 
#' Get the predicted f4 statistics for a fitted graph.
#' 
#' @param object  The fitted object.
#' @param full    Should the fitted values include the full data used for fitting?
#' @param ...     Additional arguments.
#' 
#' @export
fitted.agraph_fit <- function(object, full = TRUE, ...) {
  if (full) object$fit_data
  else      object$fit_data$graph_f4
}

#' Extract the individual errors in a fitted graph.
#' 
#' Get D - graph_f4 for each data point used in the fit.
#' 
#' @param object  The fitted object.
#' @param ...     Additional arguments.
#' 
#' @export
residuals.agraph_fit <- function(object, ...) {
  with(object$fit_data, D - graph_f4)
}

#' Predict statistics on new data.
#' 
#' Predict expected f4 statistics. If \code{newdata} is not specified this function just
#' returns the predicted values on the original data.
#' 
#' @param object  The fitted object.
#' @param newdata New data frame to predict values for.
#' @param ...     Additional arguments.
#' 
#' @export
predict.agraph_fit <- function(object, newdata = NULL, ...) {
  if (is.null(newdata)) {
    fitted(object, full = TRUE)
  } else {
    add_graph_f4(newdata, object$graph, object$fit_env)
  }
}
