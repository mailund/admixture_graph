

#' Build a parent incidene matrix from an edge list
#' 
#' @param nodes The name of the nodes in the admxture graph
#' @param parent_edges An n x 2 matrix where the first column is the child node and the second the parent.
#'                          
#' @return An incidence matrix for the parent structure of an admixture graph.
agraph_parents <- function(nodes, parent_edges) {
  n <- length(nodes)
  parents <- matrix(FALSE, n, n)
  rownames(parents) <- colnames(parents) <- nodes
  for (row in 1:nrow(parent_edges)) {
    parents[edges[row,1], edges[row,2]] = TRUE
  }
  parents
}

#' Build matrix of admixture proportions from an edge list
#' 
#' @param nodes The name of the nodes in the admxture graph
#' @param admixture_weights An n x 3 matrix where the first column is the child node, the second the parent
#'                          and the third the admixture weight on that edge.
#'                          
#' @return A matrix containing the admixture weights.
agraph_weights <- function(nodes, admixture_weights) {
  n <- length(nodes)
  weights <- matrix("", n, n)
  rownames(weights) <- colnames(weights) <- nodes
  for (row in 1:nrow(admixture_weights)) {
    weights[admixture_weights[row,1], admixture_weights[row,2]] = admixture_weights[row,3]
  }
  weights
}

#' Build the children incidence matrix from the parents incidence matrix.
#' 
#' @param parents The incidence matrix for the parent structure in an admixture graph.
#' 
#' @return An incidence matrix for the children structure.
agraph_children <- function(parents) {
  nodes <- rownames(parents)
  n <- length(nodes)
  children <- matrix(FALSE, n, n)
  rownames(children) <- colnames(children) <- nodes
  
  for (child in nodes) {
    node_parents <- parents[child,]
    node_parents <- names(node_parents[node_parents > 0])
    for (p in node_parents) {
      children[p, child] = parents[child, p]
    }
  }
  
  children
}

#' Create an admixture graph object.
#' 
#' @param nodes The name of the nodes in the admxture graph
#' @param parent_edges An n x 2 matrix where the first column is the child node and the second the parent.
#' @param admixture_weights An n x 3 matrix where the first column is the child node, the second the parent
#'                          and the third the admixture weight on that edge.
#'
#' @return An admixture graph object.
#' 
#' @examples
#' nodes <- c("A", "B", "C", "AB", "BC", "ABC", "R", "O")
#' edges <- matrix(ncol = 2, byrow=TRUE,
#'                 data = c("A", "AB",
#'                          "C", "BC", 
#'                          "B", "AB", 
#'                          "B", "BC", 
#'                          "AB", "ABC",
#'                          "BC", "ABC",
#'                          "ABC", "R",
#'                          "O", "R"))
#' admixture_proportions <- matrix(ncol = 3, byrow=TRUE,
#'                                 data = c("B", "AB", "a", 
#'                                          "B", "BC", "(1-a)"))
#' 
#' 
#' graph <- agraph(nodes, edges, admixture_proportions)
#' 
#' @export
agraph <- function(nodes, parent_edges, admixture_proportions) {
  parents <- agraph_parents(nodes, parent_edges)
  weigths <- agraph_weights(nodes, admixture_proportions)
  
  graph <- structure(list(parents = parents, probs = weights,
                          children = agraph_children(parents)), 
                     class = "agraph")
  graph
}

