
#' Extract all the parameters a graph contains.
#' 
#' The graph is parameterized by edge lengths and admixture proportions. This function extracts these parameters.
#'  
#' @param graph The admixture graph
#' 
#' @return A list containing two values: "edges", a vector of edges and "admix_prop", a vector containing admixture proportions.
#' @export
extract_graph_parameters <- function(graph) {
  edges <- c()
  for (node in graph$nodes) {
    children <- graph$nodes[graph$children[node,]]
    if (length(children) > 0) {
      edges <- c(edges, paste("edge_", node, "_", children, sep=''))
    }
  }
  
  admix_prop <- unique(Filter(function(x) nchar(x) > 0 && substr(x, 1, 1) != '(', graph$probs))
  
  list(edges = edges, admix_prop = admix_prop)
}

#' Build an environment in which f statistics can be evaluated.
#' 
#' Constructs an environment in which the f statistics for a graph can be evaluted, based on the parameters in a graph
#' and values for edge lengths and admixture proportions (with defaults if not specified).
#'  
#' @param parameters     The parameters of a graph as returned by \code{extract_graph_parameters}
#' @param edge_lengths   If specified, a vector of edge lengths. Otherwise defaults are used.
#' @param admix_prop     If specified, a vector of admixture proportions. Otherwise defaults are used.
#' 
#' @return A list containing two values: "edges", a vector of edges and "admix_prop", a vector containing admixture proportions.
#' @export
graph_environment <- function(parameters, edge_lengths = NULL, admix_prop = NULL) {
  
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
