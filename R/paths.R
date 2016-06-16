#' Create a path data frame from a list of nodes.
#' 
#' Creates a path data frame from a list of nodes.
#'
#' @param graph  The admixture graph the path is in.
#' @param nodes  A list of nodes on a path.
#'
#' @return A data frame capturing the path and the probabilities/weights on the edges.
format_path <- function(graph, nodes) {
  path <- data.frame(from = graph$nodes[nodes[1 : length(nodes) - 1]],
                     to = graph$nodes[nodes[2 : length(nodes)]],
                     stringsAsFactors = FALSE)
  path$prob <- mapply(function(f,t) graph$probs[f, t], path$from, path$to)
  path
}

#' Compute all paths from one leaf to another.
#'
#' Computes all paths from one leaf to another.
#'
#' @param graph  The admixture graph.
#' @param src    The starting leaf.
#' @param dst    The destination leaf.
#'
#' @return A list containing all the paths from \code{src} to \code{dst}.
#' 
#' @export
all_paths <- function(graph, src, dst) {
  if (!(src %in% graph$leaves)) stop(paste(src, "is not a leaf in the graph."))
  if (!(dst %in% graph$leaves)) stop(paste(dst, "is not a leaf in the graph."))
  
  src_idx <- which(src == rownames(graph$parents))
  dst_idx <- which(dst == rownames(graph$parents))
  PATHS <- list()
  X <- 1

  recurse_up <- function(node, path = c()) {
    up_edges <- which(graph$parents[node,] != 0)
    down_edges <- which(graph$children[node,] != 0)

    lapply(up_edges, function(e) recurse_up(e, c(path, node)))
    lapply(down_edges, function(e) recurse_down(e, c(path, node)))
  }

  recurse_down <- function(node, path) {
    if (node %in% path)  return(NULL)
    if (node == dst_idx) {
      PATHS[[X]] <<- format_path(graph, c(path, node))
      X <<- X + 1
    } else {
      down_edges <- which(graph$children[node,] != 0)
      lapply(down_edges, function(e) recurse_down(e, c(path, node)))
    }
  }

  recurse_up(src_idx)
  PATHS
}