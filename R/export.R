

no_parents <- function(graph) function(node) {
  length(which(graph$parents[node, ]))
} 

get_root <- function(graph) {
  roots <- which(Map(no_parents(graph), graph$nodes) == 0)
  if (length(roots) > 1) stop("Don't know how to handle more than one root")
  sprintf("root\t%s", names(roots[1]))
}

get_labels <- function(graph) {
  paste("label", graph$leaves, graph$leaves, sep = "\t", collapse = "\n")
}

mk_edge_label <- function(from, to) {
  paste0("edge_", from, "_", to, sep = "")
}

get_edges <- function(graph) {
  edge_statements <- c()
  nodes <- rownames(graph$parents)
  for (n in nodes) {
    parents <- names(which(graph$parents[n, ]))
    if (length(parents) == 0) next  # root
    if (length(parents) == 1) {  # normal node
      edge_label <- mk_edge_label(parents[1], n)
      edge_statements <- c(edge_statements, 
                           sprintf("edge\t%s\t%s\t%s",
                                   edge_label, parents[1], n))
      next
    }
    if (length(parents) == 2) {  # admixture node
      left <- paste0(parents[1], "_left", sep = "")
      right <- paste0(parents[2], "_right", sep = "")  
      left_label <- mk_edge_label(left, n)
      right_label <- mk_edge_label(right, n)
      edge_statements <- c(edge_statements,
                           sprintf("edge\t%s\t%s\t%s",
                                   left_label, left, n),
                           sprintf("edge\t%s\t%s\t%s",
                                   right_label, right, n)
                           )
      next
    }
    stop(paste0("Unexpected number of parents for node ", n))
  }
  edge_statements
}

get_admixtures <- function(graph) {
  admix_nodes <- names(which(rowSums(graph$parents) > 1))
  admix_statements <- vector("character", length = length(admix_nodes))
  for (i in seq_along(admix_nodes)) {
    node <- admix_nodes[i]
    parents <- names(which(graph$parents[node, ]))
    if (length(parents) != 2) stop(paste0("Unexpected parents for ", node))
    left <- paste0(parents[1], "_left", sep = "")
    right <- paste0(parents[2], "_right", sep = "")
    admix_statements[i] <- sprintf("admix\t%s\t%s\t%s\t50\t50", node, left, right)
  }
  admix_statements
}

#' Export to Patterson's qpGraph format.
#' 
#' This function writes a graph to a file-object, f, in the format used by the
#' qpGraph tool. This format takes admixture proportions as part of the specification,
#' but since we do not hold these proportions in our graphs, but only in fitted data,
#' the export functions puts all admixture proportions as 50%/50%. Edit the output
#' file by hand if you want to change this.
#' 
#' @param f      File object, e.g. stdout()
#' @param graph  A graph to export.
#' @export
export_to_qpGraph <- function(f, graph) {
  root <- get_root(graph)
  labels <- get_labels(graph)
  edges <- get_edges(graph)
  admixtures <- get_admixtures(graph)

  lines <- c(root, "", labels, "", edges, "", admixtures)
  writeLines(lines, con = f)
}
