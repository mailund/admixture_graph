

agraph_children <- function(parents) {
  nodes <- rownames(parents)
  n <- length(nodes)
  children <- matrix(0, n, n)
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

agraph <- function(parents) {
  structure(list(parents = parents, children = agraph_children(parents)), 
            class = "agraph")
}


all_paths <- function(graph, src, dst) {
  src_idx <- which(src == rownames(graph$parents))
  dst_idx <- which(dst == rownames(graph$parents))
  PATHS <- list()
  X <- 1
  
  recurse_up <- function(node, path = c()) {
    up_edges <- which(parents[node,] > 0)
    down_edges <- which(children[node,] > 0)
    
    sapply(up_edges, function(e) recurse_up(e, c(path, node)))
    sapply(down_edges, function(e) recurse_down(e, c(path, node)))
  }
  
  recurse_down <- function(node, path) {
    if (node %in% path)  return(NULL)
    if (node == dst_idx) {
      PATHS[[X]] <<- path
      X <<- X + 1
    } else {
      down_edges <- which(children[node,] > 0)
      sapply(down_edges, function(e) recurse_down(e, c(path, node)))
    }
  }
  
  recurse_up(src_idx)
  PATHS
}




alpha <- 0.1
parents <- matrix(0.0, 6,6)
rownames(parents) <- colnames(parents) <- c("A", "B", "C", "AB", "BC", "R")
parents["A","AB"] <- 1.0
parents["C","BC"] <- 1.0
parents["B","AB"] <- alpha
parents["B","BC"] <- 1.0 - alpha
parents["AB","R"] <- 1.0
parents["BC","R"] <- 1.0

G <- agraph(parents)


AtoC <- all_paths(G, "A", "C")
AtoB <- all_paths(G, "A", "B")

