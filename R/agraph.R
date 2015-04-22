

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
    parents[parent_edges[row,1], parent_edges[row,2]] = TRUE
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
  cat(weights, "\n")
  for (row in 1:nrow(admixture_weights)) {
    weights[admixture_weights[row,1], admixture_weights[row,2]] = admixture_weights[row,3]
    weights[admixture_weights[row,2], admixture_weights[row,1]] = admixture_weights[row,3]
  }
  weights
}

#' Build a children incidene matrix from an parent edge list
#' 
#' @param nodes The name of the nodes in the admxture graph
#' @param parent_edges An n x 2 matrix where the first column is the child node and the second the parent.
#'                          
#' @return An incidence matrix for the children structure of an admixture graph.
agraph_children <- function(nodes, parent_edges) {
  n <- length(nodes)
  children <- matrix(FALSE, n, n)
  rownames(children) <- colnames(children) <- nodes
  for (row in 1:nrow(parent_edges)) {
    children[parent_edges[row,2], parent_edges[row,1]] = TRUE
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
  children <- agraph_children(nodes, parent_edges)
  admixture_probs <- agraph_weights(nodes, admixture_proportions)
  structure(list(parents = parents, probs = admixture_probs, children = children), class = "agraph")
}

