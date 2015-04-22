
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

#' Calculate the f4(W,X;Y,Z) statistics.
#'  
#' @param graph The admixture graph
#' @param W     A leaf node
#' @param X     A leaf node
#' @param Y     A leaf node
#' @param Z     A leaf node
#' 
#' @return A symbolic representation of the equation for the f4 statistics given by the admixture graph.
#' @export
f4 <- function(graph, W, X, Y, Z) {
  WX <- all_paths(graph, W, X)
  YZ <- all_paths(graph, Y, Z)
  all_overlaps_weights(WX, YZ)
}

#' Calculate the f4(A;B,C) statistics.
#'  
#' @param graph The admixture graph
#' @param A     A leaf node
#' @param B     A leaf node
#' @param C     A leaf node
#' 
#' @return A symbolic representation of the equation for the f3 statistics given by the admixture graph.
#' @export
f3 <- function(graph, A, B, C) f4(graph, A, B, A, C)

#' Calculate the f2(A,B) statistics.
#'  
#' @param graph The admixture graph
#' @param A     A leaf node
#' @param B     A leaf node
#' 
#' @return A symbolic representation of the equation for the f2 statistics given by the admixture graph.
#' @export
f2 <- function(graph, A, B) f4(graph, A, B, A, B)

