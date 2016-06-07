#' Extract all the parameters a graph contains.
#' 
#' The graph is parameterized by edge lengths and admixture proportions. This 
#' function extracts these parameters.
#' 
#' @param graph  The admixture graph.
#'   
#' @return A list containing two values: \code{edges}, a vector of edges and 
#'         \code{admix_prop}, a vector containing admixture proportions.
#'         
#' @export
extract_graph_parameters <- function(graph) {
  edges <- c()
  for (node in graph$nodes) {
    children <- graph$nodes[graph$children[node,]]
    if (length(children) > 0) {
      edges <- c(edges, paste("edge_", node, "_", children, sep=""))
    }
  }

  admix_prop_variable <- function(x) nchar(x) > 0 && substr(x, 1, 1) != "("
  admix_prop <- unique(Filter(admix_prop_variable, graph$probs))

  list(edges = edges, admix_prop = admix_prop)
}

#' Build an environment in which f statistics can be evaluated.
#' 
#' Constructs an environment in which the \eqn{f} statistics for a graph can be
#' evaluted, based on the parameters in a graph and values for edge lengths and
#' admixture proportions (with defaults if not specified).
#' 
#' @param parameters    The parameters of a graph as returned by
#'                      \code{\link{extract_graph_parameters}}.
#' @param edge_lengths  If specified, a vector of edge lengths. Otherwise
#'                      defaults are used.
#' @param admix_prop    If specified, a vector of admixture proportions.
#'                      Otherwise defaults are used.
#'   
#' @return A list containing two values: \code{edges}, a vector of edges and
#'         \code{admix_prop}, a vector containing admixture proportions.
#'   
#' @export
graph_environment <- function(parameters,
                              edge_lengths = NULL,
                              admix_prop = NULL) {

  if (is.null(edge_lengths)) {
    edge_lengths <- rep(0.01, length(parameters$edges))
  }
  if (is.null(admix_prop)) {
    admix_prop <- rep(0.5, length(parameters$admix_prop))
  }

  env <- new.env()
  for (i in seq_along(parameters$edges)) {
    assign(parameters$edges[i], edge_lengths[i], env)
  }
  for (i in seq_along(parameters$admix_prop)) {
    assign(parameters$admix_prop[i], admix_prop[i], env)
  }
  env
}

#' Evaluates an f_4 statistics in a given environment.
#' 
#' Evaluates an \eqn{f_4} statistics in a given environment.
#' 
#' @param graph  The admixture graph.
#' @param env    The environment containing the graph parameters.
#' @param W      First population/sample.
#' @param X      Second population/sample.
#' @param Y      Third population/sample.
#' @param Z      Fourth population/sample.
#'   
#' @return The \eqn{f_4} value specified by the graph and the environment.
#'   
#' @export
evaluate_f4 <- function(graph, env, W, X, Y, Z) {
  eval(sf4(graph, W, X, Y, Z), env)
}

#' Evalutes the f_4 statistics for all rows in a data frame and extends 
#' the data frame with the graph f_4.
#' 
#' The data frame, \code{data}, must contain columns \code{W}, \code{X}, 
#' \code{Y}, and \code{Z}. The function then computes the \eqn{f_4(W, X; Y, Z)}
#' statistics (also known as the \eqn{D} statistics) for all rows and adds these
#' as a column, \code{graph_f4}, to the data frame.
#' 
#' @param data   The data frame to get the labels to compute the \eqn{f_4} statistics from.
#' @param graph  The admixture graph.
#' @param env    The environment to evaluate the \eqn{f_4} statistics in.
#'   
#' @return A data frame identical to \code{data} except with an additional 
#'         column, \code{graph_f4}, containing the \eqn{f_4} values as determined by
#'         the graph and the environment.
#'   
#' @export
add_graph_f4 <- function(data, graph, env) {
  data$graph_f4 <- unlist(Map(function(W,X,Y,Z) evaluate_f4(graph, env, W, X, Y, Z),
                          data$W, data$X, data$Y, data$Z),
                          use.names = FALSE)
  data
}