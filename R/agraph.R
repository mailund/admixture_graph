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
#' @export
agraph <- function(parents, probs) {
  graph <- structure(list(parents = parents, probs = probs,
                          children = agraph_children(parents)), 
                     class = "agraph")
  graph
}

