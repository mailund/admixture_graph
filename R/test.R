

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

agraph <- function(parents, probs) {
  graph <- structure(list(parents = parents, probs = probs,
                     children = agraph_children(parents)), 
                     class = "agraph")
  graph
}

format_path <- function(graph, nodes) {
  path <- data.frame(from = nodes[1:length(nodes)-1], to = nodes[2:length(nodes)])
  path$prob <- mapply(function(f,t) graph$probs[f,t], path$from, path$to)
  path
}

all_paths <- function(graph, src, dst) {
  src_idx <- which(src == rownames(graph$parents))
  dst_idx <- which(dst == rownames(graph$parents))
  PATHS <- list()
  X <- 1
  
  recurse_up <- function(node, path = c()) {
    up_edges <- which(graph$parents[node,] != 0)
    down_edges <- which(graph$children[node,] != 0)
    
    sapply(up_edges, function(e) recurse_up(e, c(path, node)))
    sapply(down_edges, function(e) recurse_down(e, c(path, node)))
  }
  
  recurse_down <- function(node, path) {
    if (node %in% path)  return(NULL)
    if (node == dst_idx) {
      PATHS[[X]] <<- format_path(graph, c(path, node))
      X <<- X + 1
    } else {
      down_edges <- which(graph$children[node,] != 0)
      sapply(down_edges, function(e) recurse_down(e, c(path, node)))
    }
  }
  
  recurse_up(src_idx)
  PATHS
}

path_probability <- function(path) {
  paste(path$prob, sep=" ", collapse="")
}

path_overlap <- function(path1, path2) {
  probability = paste(path_probability(path1), path_probability(path2), sep='', collapse='')
  rev_path2 <- data.frame(from = path2$to, to = path2$from, prob = path2$prob)
  
  list(prob = probability,
       positive = merge(path1, path2)[,1:2],
       negative = merge(path1, rev_path2)[,1:2])
}

path_overlap_weight <- function(overlap) {
  if (nrow(overlap$positive) > 0) {
    positive = paste(paste('[', overlap$positive$from, ':', overlap$positive$to, ']', sep=''),
                     collapse=' + ')  
  } else {
    positive = ''
  }
  if (nrow(overlap$negative) > 0) {
    negative = paste(paste('[', overlap$negative$from, ':', overlap$negative$to, ']', sep=''),
                     collapse=' - ')  
  } else {
    negative = ''
  }
  
  if (positive != '' && negative != '') {
    paste(overlap$prob, ' * (', positive, ' - ', negative, ')', sep='')  
  } else if (positive == '' && negative != '') {
    paste('-', overlap$prob, ' * (', negative, ')', sep='')  
  } else if (positive != '' && negative == '') {
    paste(overlap$prob, ' * (', positive, ')', sep='')  
  } else {
    ''
  }  
}

all_overlaps_weights <- function(paths1, paths2) {
  overlap_weights <- c()
  for (i in 1:length(paths1)) {
    for (j in 1:length(paths2)) {
      overlap <- path_overlap(paths1[[i]], paths2[[j]])
      overlap_weights <- c(path_overlap_weight(overlap), overlap_weights)
    }
  }
  paste(Filter(function(x) x != "", overlap_weights), sep=" + ")
}

f4 <- function(graph, W, X, Y, Z) {
  WX <- all_paths(graph, W, X)
  YZ <- all_paths(graph, Y, Z)
  all_overlaps_weights(WX, YZ)
}

f3 <- function(graph, A, B, C) f4(graph, A, B, A, C)
f2 <- function(graph, A, B)    f4(graph, A, B, A, B)

parents <- matrix(FALSE, 8, 8)
rownames(parents) <- colnames(parents) <- c("A", "B", "C", "AB", "BC", "ABC", "R", "O")
parents["A","AB"] <- TRUE
parents["C","BC"] <- TRUE
parents["B","AB"] <- TRUE
parents["B","BC"] <- TRUE
parents["AB","ABC"] <- TRUE
parents["BC","ABC"] <- TRUE
parents["ABC","R"] <- TRUE
parents["O","R"] <- TRUE

weights <- matrix("", 8, 8)
rownames(weights) <- colnames(weights) <- c("A", "B", "C", "AB", "BC", "ABC", "R", "O")
weights["B","AB"] <- weights["AB","B"] <- 'a'
weights["B","BC"] <- weights["BC","B"] <- '(1-a)'


graph <- agraph(parents, weights)
f4(graph, "O", "A", "B", "C")
f3(graph, "O", "A", "C")
f3(graph, "O", "A", "B")