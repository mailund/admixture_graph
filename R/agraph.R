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
      stop(paste("The node", parent_edges[row,1], 
                 "is used in the edges but is not specified as a node."))
    }
    if (!(parent_edges[row,2] %in% nodes)) {
      stop(paste("The node", parent_edges[row,2], 
                 "is used in the edges but is not specified as a node."))
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
      stop(paste("The node", parent_edges[row,1], 
                 "is used in the edges but is not specified as a node."))
    }
    if (!(parent_edges[row,2] %in% nodes)) {
      stop(paste("The node", parent_edges[row,2], 
                 "is used in the edges but is not specified as a node."))
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
edge <- function(child, parent) c(child, parent, NA)

#' Create an admixture edge from a child to two parents.
#' 
#' Syntactic suggar for constructing edges in an admixture graph.
#' 
#' @param child    The name of the child node.
#' @param parent1  The name of the parent node.
#' @param parent2  The name of the parent node.
#' @param prop     Admixture proportions from \code{parent1} to \code{child}.
#'                 If this parameter is not provided, you must explicitly
#'                 specify the admixture proportion parameters in the
#'                 \code{agraph} function call.
#' 
#' @export
admixture_edge <- function(child, parent1, parent2, prop = NA) {
  if (is.na(prop)) {
    other_prop <- NA
  } else {
    other_prop <- paste0("(1 - ", prop, ")")
  }
  c(child, parent1, prop, child, parent2, other_prop)
}

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
  c(child, parent1, prop, child, parent2, paste0("(1 - ", prop, ")"))

#' Extract the admixture proportion parameter from edge specifications.
#' 
#' This function is simply selecting the edges with admixture proportion specifications
#' so these can be handled when building a graph using \code{agraph}. It is not a function
#' you would need to call explicitly, rather it is there to allow people \emph{not} to use
#' it to provide admixture proportions explicitly (which we normally wouldn't recommend).
#' 
#' @param parent_edges  Matrix created with the \code{agraph_parents} function.
#' @return The parents edges reduced to the rows with admixture proportions.
extract_admixture_proportion_parameters <- function(parent_edges) {
  result <- parent_edges[!is.na(parent_edges[,3]),]
  if(nrow(result) > 0) result else NULL
}

#' Create the list of edges for an admixture graph.
#' 
#' Syntactic suggar for constructing edges in an admixture graph.
#' 
#' @param edges  List of edges.
#' 
#' @export
parent_edges <- function(edges) matrix(ncol = 3, byrow = TRUE, data = edges)

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
agraph <- function(leaves, inner_nodes, parent_edges, 
                   admixture_proportions = extract_admixture_proportion_parameters(parent_edges)) {
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

#' Extract trees
#'
#' Extracts all the trees embedded in an agraph object
#' 
#' @param graph  An agraph object
#' 
#' @return  A list of trees
#' 
#' @export
extract_trees <- function(graph) {
  to_do_list <- list(graph)
  tree_list <- list()
  while (length(to_do_list) > 0) {
    G <- to_do_list[[1]]
    to_do_list[[1]] <- NULL
    if (length(extract_graph_parameters(G)$admix_prop) > 0) {
      new_graphs <- split_first_admixture(G)
      to_do_list[[length(to_do_list) + 1]] <- new_graphs[[1]]
      to_do_list[[length(to_do_list) + 1]] <- new_graphs[[2]]
    } else {
      tree_list[[length(tree_list) + 1]] <- remove_joints_from_a_tree(G)
    }
  }
  return(tree_list)
}

split_first_admixture <- function(graph) {
  split <- FALSE
  nodes <- graph$nodes
  leaves <- graph$leaves
  inner_nodes <- graph$inner_nodes
  edge_vector1 <- character(0)
  edge_vector2 <- character(0)
  for (i in seq(1, NROW(graph$parents))) {
    row <- graph$parents[i, ]
    if (length(row[row == TRUE]) == 0) {
      # This is the root row.
    } else if (length(row[row == TRUE]) == 1) {
      # Keep the normal edge.
      edge_vector1 <- c(edge_vector1, edge(nodes[i], nodes[which(row == TRUE)[1]]))
      edge_vector2 <- c(edge_vector2, edge(nodes[i], nodes[which(row == TRUE)[1]]))
    } else if (length(row[row == TRUE]) == 2) {
      # Either split the admixture into two or keep it.
      if (split == FALSE) {
        # Splitting.
        edge_vector1 <- c(edge_vector1, edge(nodes[i], nodes[which(row == TRUE)[1]]))
        edge_vector2 <- c(edge_vector2, edge(nodes[i], nodes[which(row == TRUE)[2]]))
        split <- TRUE
      } else {
        # The names of the admixture proportions can be forgotten!
        edge_vector1 <- c(edge_vector1, admixture_edge(nodes[i], nodes[which(row == TRUE)[1]],
                                                       nodes[which(row == TRUE)[2]], i))
        edge_vector2 <- c(edge_vector2, admixture_edge(nodes[i], nodes[which(row == TRUE)[1]],
                                                       nodes[which(row == TRUE)[2]], i))
      }
    }
  }
  graph1 <- agraph(leaves, inner_nodes, parent_edges(edge_vector1))
  graph2 <- agraph(leaves, inner_nodes, parent_edges(edge_vector2))
  # The graphs might now contain joints, inner nodes with one parent and one child.
  # They are easier to remove once we have no admixtures left, so for now we leave them be.
  return(list(graph1, graph2))
}

remove_joints_from_a_tree <- function(tree) {
  # The joints must be removed one at a time.
  work_left <- TRUE
  while (work_left == TRUE) {
    nodes <- tree$nodes
    count <- 0
    for (i in seq(1, NROW(tree$parents))) {
      row <- tree$parents[i, ]
      column <- tree$parents[, i]
      if (length(row[row == TRUE]) == 1 && length(column[column == TRUE]) == 1) {
        # There is at least one joint left; we call it remove and remove it during this iteration.
        remove <- nodes[i]
        count <- count + 1
      }
    }
    if (count == 0) {
      work_left <- FALSE # All clear.
    } else {
      leaves <- tree$leaves
      inner_nodes <- tree$inner_nodes
      inner_nodes <- inner_nodes[-which(inner_nodes == remove)]
      edge_vector <- character(0)
      for (i in seq(1, NROW(tree$parents))) {
        row <- tree$parents[i, ]
        column <- tree$parents[, i]
        if (nodes[i] == remove) {
          edge_vector <- c(edge_vector, edge(nodes[which(column == TRUE)[1]],
                                             nodes[which(row == TRUE)[1]]))
        } else {
          if (length(row[row == TRUE]) == 1 && row[remove] == FALSE) {
            edge_vector <- c(edge_vector, edge(nodes[i], nodes[which(row == TRUE)[1]]))
          }
        }
      }  
    }
    tree <- agraph(leaves, inner_nodes, parent_edges(edge_vector))
  }
  return(tree)
}