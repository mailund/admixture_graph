
#' Extract all the parameters a graph contains.
#' 
#' The graph is parameterized by edge lengths and admixture proportions. This function extracts these parameters.
#'  
#' @param graph The admixture graph
#' 
#' @return A list containing two values: "edges", a vector of edges and "admix_prob", a vector containing admixture proportions.
#' @export
extract_graph_parameters <- function(graph) {
  edges <- c()
  for (node in graph$nodes) {
    children <- graph$nodes[graph$children[node,]]
    if (length(children) > 0) {
      edges <- c(edges, paste("edge_", node, "_", children, sep=''))
    }
  }
  
  admix_prob <- unique(Filter(function(x) nchar(x) > 0 && substr(x, 1, 1) != '(', graph$probs))
  
  list(edges = edges, admix_prob = admix_prob)
}
