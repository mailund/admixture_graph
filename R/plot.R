
#' Plot an admixture graph.
#' 
#' This is a basic drawing routine for visualising the graph. For publication quality graphs a lot
#' more tweaking is probably needed.
#' 
#' @param ordered_leaves The leaf-nodes in the left to right order they should be drawn.
#'                       I don't have a good algorithm for figuring out that order so for now
#'                       it is required as a function argument.
#'                       
#' @export
plot.agraph <- function(graph, ordered_leaves) {
  
  dfs <- function(node, basis, step) {
    result <- rep(NA, length(graph$nodes))
    names(result) <- graph$nodes
    
    dfs <- function(node) {
      children <- which(graph$children[node,])
      if (length(children) == 0) {
        result[node] <<- basis(node)
      } else {
        result[node] <<- step(vapply(children, dfs, numeric(1)))
      } 
    }
    dfs(node)
    result
  }
  
  root <- "R"
  ypos <- dfs(root, basis = function(x) 0.0, step = function(x) max(x) + 1.0)
  
  leaf_index <- function(n) which(graph$nodes[n] == ordered_leaves)
  left_x  <- dfs(root, basis = leaf_index, step = min)
  right_x <- dfs(root, basis = leaf_index, step = max)
  xpos <- left_x + (right_x - left_x)/2.0
  
  # Start the actual drawing of the graph...
  plot(xpos, ypos, type='n', axes=FALSE, frame.plot=FALSE, xlab='', ylab='', ylim=c(-1,max(ypos)+0.5))
  
  for (node in graph$nodes) {
    parents <- graph$nodes[graph$parents[node,]]
    if (length(parents) == 1) {
      lines(c(xpos[node],xpos[parents]), c(ypos[node], ypos[parents]))
      
    } else if (length(parents) == 2) {
      break_y <- ypos[node]
      break_x_left <- xpos[node] - 0.3
      break_x_right <- xpos[node] + 0.3
      
      lines(c(xpos[parents[1]], break_x_left),  c(ypos[parents[1]], break_y))
      lines(c(xpos[parents[2]], break_x_right), c(ypos[parents[2]], break_y))
      segments(break_x_left, break_y, xpos[node],  ypos[node], col='red')
      segments(break_x_right, break_y, xpos[node], ypos[node], col='red')
      text(break_x_left, break_y, graph$probs[parents[[1]],node], cex=0.5, pos=1, col='red', offset=0.1)
      text(break_x_right, break_y, graph$probs[parents[[2]],node], cex=0.5, pos=1, col='red', offset=0.1)
    }
  }
  
  is_inner <- Vectorize(function(n) sum(graph$children[n,]) > 0)
  inner_nodes <- which(is_inner(graph$nodes))
  leaves <- which(!is_inner(graph$nodes))
  
  text(xpos[inner_nodes], ypos[inner_nodes], labels = graph$nodes[inner_nodes], cex=0.6, col='blue', pos=3)
  text(xpos[leaves], ypos[leaves], labels = graph$nodes[leaves], cex=0.7, col='black', pos=1)
  #points(xpos, ypos, pch=20)
  
  inner_nodes
}

plot(graph, c("BLK", "APB", "PB", "BC", "A", "YB", "BB", "EBB"))
