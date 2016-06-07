#' Build the parent incidence matrix from an edge list.
#' 
#' Build the parent incidence matrix from an edge list.
#' 
#' @param nodes         Nodes in the admxture graph.
#' @param parent_edges  An \eqn{n \times 2}{n x 2} matrix where the first column is the
#'                      child node and the second is the parent.
#'   
#' @return An incidence matrix for the parent structure of an admixture graph.
#' 
#' @seealso \code{\link{agraph_children}}
agraph_parents <- function(nodes, parent_edges) {
  n <- length(nodes)
  parents <- matrix(FALSE, n, n)
  rownames(parents) <- colnames(parents) <- nodes
  for (row in 1:nrow(parent_edges)) {
    if (!(parent_edges[row,1] %in% nodes)) {
      stop(paste("The node", parent_edges[row,1], "is used in the edges but is not specified as a node."))
    }
    if (!(parent_edges[row,2] %in% nodes)) {
      stop(paste("The node", parent_edges[row,2], "is used in the edges but is not specified as a node."))
    }
    parents[parent_edges[row,1], parent_edges[row,2]] <- TRUE
  }
  parents
}

#' Build the child incidene matrix from an parent edge list.
#' 
#' Build the child incidene matrix from an parent edge list.
#' 
#' @param nodes         Nodes in the admxture graph.
#' @param parent_edges  An \eqn{n \times 2}{n x 2} matrix where the first column is the 
#'                      child node and the second is the parent.
#'   
#' @return An incidence matrix for the child structure of an admixture graph.
#' 
#' @seealso \code{\link{agraph_parents}}
agraph_children <- function(nodes, parent_edges) {
  n <- length(nodes)
  children <- matrix(FALSE, n, n)
  rownames(children) <- colnames(children) <- nodes
  for (row in 1:nrow(parent_edges)) {
    if (!(parent_edges[row,1] %in% nodes)) {
      stop(paste("The node", parent_edges[row,1], "is used in the edges but is not specified as a node."))
    }
    if (!(parent_edges[row,2] %in% nodes)) {
      stop(paste("The node", parent_edges[row,2], "is used in the edges but is not specified as a node."))
    }
    children[parent_edges[row,2], parent_edges[row,1]] <- TRUE
  }
  children
}

#' Build the matrix of admixture proportions from an edge list.
#' 
#' Build the matrix of admixture proportions from an edge list.
#' 
#' @param nodes              The name of the nodes in the admxture graph.
#' @param admixture_weights  An \eqn{n \times 3}{n x 3} matrix where the first column is
#'                           the child node, the second isthe parent and the third
#'                           is the admixture weight on that edge.
#' @param parents            The parent edge list. Used for checking graph consistency.
#'   
#' @return A matrix containing the admixture weights.
agraph_weights <- function(nodes, admixture_weights, parents) {
  n <- length(nodes)
  weights <- matrix("", n, n)
  rownames(weights) <- colnames(weights) <- nodes
  if (!is.null(admixture_weights)) {
    for (row in 1:nrow(admixture_weights)) {
      child <- admixture_weights[row,1]
      parent <- admixture_weights[row,2]
      admix_weight <- admixture_weights[row,3]
      
      if (! parents[child,parent]) {
        stop(paste("There is no edge in the graph from", child, "to", parent))
      }
      
      weights[child, parent] <- admix_weight
      weights[parent, child] <- admix_weight
    }
  }
  weights
}

#' Create an edge from a child to a parent.
#' 
#' Syntactic suggar for constructing edges in an admixture graph.
#' 
#' @param child   The name of the child node.
#' @param parent  The name of the parent node.
#' 
#' @export
edge <- function(child, parent) c(child, parent)

#' Create an admixture edge from a child to two parents.
#' 
#' Syntactic suggar for constructing edges in an admixture graph.
#' 
#' @param child    The name of the child node.
#' @param parent1  The name of the parent node.
#' @param parent2  The name of the parent node.
#' 
#' @export
admixture_edge <- function(child, parent1, parent2) c(child, parent1, child, parent2)

#' Specify the proportions in an admixture event.
#' 
#' Syntactic suggar for constructing edges in an admixture graph.
#' 
#' @param child    The child node.
#' @param parent1  The first parent.
#' @param parent2  The second parent.
#' @param prop     The admixture proportions coming from the first parent.
#' 
#' @export
admix_props <- function(child, parent1, parent2, prop)
  c(child, parent1, prop, child, parent2, paste("(1 - ", prop, ")", sep=""))

#' Create the list of edges for an admixture graph.
#' 
#' Syntactic suggar for constructing edges in an admixture graph.
#' 
#' @param edges  List of edges.
#' 
#' @export
parent_edges <- function(edges) matrix(ncol = 2, byrow = TRUE, data = edges)

#' Create the list of admixture proportions for an admixture graph.
#' 
#' Syntactic suggar for constructing edges in an admixture graph.
#' 
#' @param admix_props  The admixture proportions.
#' 
#' @export
admixture_proportions <- function(admix_props) matrix(ncol = 3, byrow = TRUE, data = admix_props)

#' Create an admixture graph object.
#' 
#' Create an admixture graph object, an acyclic directed graph.
#' 
#' @param leaves                 The names of the leaves in the admixture graph.
#' @param inner_nodes            The name of the inner nodes in the admxture graph.
#' @param parent_edges           The list of edges in the graph, created by 
#'                               \code{\link{parent_edges}}.
#' @param admixture_proportions  The list of admixture proportions; created by 
#'                               \code{\link{admixture_proportions}}.
#'
#' @return An admixture graph object.
#' 
#' @seealso \code{\link{edge}}
#' @seealso \code{\link{admixture_edge}}
#' @seealso \code{\link{admix_props}}
#' @seealso \code{\link{parent_edges}}
#' @seealso \code{\link{admixture_proportions}}
#' @seealso \code{\link{plot.agraph}}
#'   
#' @examples
#' leaves <- c("A", "B", "C", "D")
#' inner_nodes <- c("ab", "b", "bc", "abc", "abcd")
#' edges <- parent_edges(c(edge("A", "ab"),
#'                         edge("B", "b"),
#'                         edge("C", "bc"),
#'                         edge("D", "abcd"),
#'                         edge("ab", "abc"),
#'                         edge("bc", "abc"),
#'                         edge("abc", "abcd"),
#'                         admixture_edge("b", "ab", "bc")))
#' admixtures <- admixture_proportions(c(admix_props("b", "ab", "bc", "x")))
#' 
#' graph <- agraph(leaves, inner_nodes, edges, admixtures)
#' 
#' @export
agraph <- function(leaves, inner_nodes, parent_edges, admixture_proportions) {
  nodes <- c(leaves, inner_nodes)
  parents <- agraph_parents(nodes, parent_edges)
  children <- agraph_children(nodes, parent_edges)
  admixture_probs <- agraph_weights(nodes, admixture_proportions, parents)
  structure(list(leaves = leaves, inner_nodes = inner_nodes,
                 nodes = nodes,
                 parents = parents,
                 probs = admixture_probs,
                 children = children),
            class = "agraph")
}