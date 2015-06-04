## NB! I HAVE NOT MADE ANY ATTEMPTS AT OPTIMISING THIS CODE. THE PACKING AND UNPACKING IS COSTLY
## AS IS COMPUTING THE SYMBOLIC REPRESENTATION OF STATISTICS FOR EACH STATISTICS. THERE IS PLENTY
## OF  WAYS THIS COULD BE OPTIMISED -- FOR NOW I JUST AIM AT GETTING IT TO WORK.


pack_environment <- function(parameters, env) {
  n_edges <- length(parameters$edges)
  n_admix_prop <- length(parameters$admix_prop)
  x <- rep(NA, n_edges + n_admix_prop)
  for (i in seq_along(parameters$edges)) {
    x[i] <- get(parameters$edges[i], env)
  }
  for (i in seq_along(parameters$admix_prop)) {
    x[i+n_edges] <- get(parameters$admix_prop[i], env)
  }
  x
}

unpack_environment <- function(parameters, x) {
  n_edges <- length(parameters$edges)
  n_admix_prop <- length(parameters$admix_prop)
  edges <- x[1:n_edges]
  admix_prop <- x[(n_edges+1):(n_edges+n_admix_prop)]
  graph_environment(parameters, edges, admix_prop)
}

make_cost_function <- function(data, graph, parameters = extract_graph_parameters(graph)) {
  force(data)
  force(graph)
  force(parameters)
  function(x) {
    env <- unpack_environment(parameters, x)
    predictions <- with(data, Vectorize(function(W,X,Y,Z) evaluate_f4(graph, env, W, X, Y, Z))(W, X, Y, Z))
    sum((data$D - predictions)**2)
  }
}

#' Fit the graph parameters to a data set.
#' 
#' Tries to minimize the squared distance between statistics in \code{data} and statistics given
#' by the graph.
#' 
#' The data frame, \code{data}, must contain columns \code{W}, \code{X}, \code{Y}, and \code{Z}. The function then computes
#' the f4(W,X;Y,Z) statistics for all rows from these to obtain the prediction made by the graph.
#' 
#' The data frame must also contain a column, \code{D}, containing the statistics observed in the data. The fitting
#' algorithm attempts to minimize the distance from this column and the predictions made by the graph.
#' 
#' @param data  The data set.
#' @param graph The admixture graph.
#' @param optimisation_options  Options to the optimisation algorithm. 
#' 
#' @return FIXME
#' 
#' @seealso neldermead::optimset
fit_graph <- function(data, graph, optimisation_options = NULL) {
  if (!requireNamespace("neldermead", quietly = TRUE)) {
    stop("This function requires neldermead to be installed.")
  }
  
  params <- extract_graph_parameters(graph)
  init_env <- graph_environment(params)
  x0 <- pack_environment(params, init_env)
  cfunc <- make_cost_function(data, graph, params)
  
  opti <- fminbnd(cfunc, x0=x0, xmin=rep(0,length(x0)), xmax=rep(1,length(x0)), options=optimisation_options)
  
  best_fit <- neldermead.get(opti, 'xopt')
  best_fit_env <- unpack_environment(params, best_fit)
  
  add_graph_f4(data, graph, best_fit_env)
}

