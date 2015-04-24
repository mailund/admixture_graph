#' Collects the postive and negative overlap between two paths.
#' 
#' @param path1 The first path
#' @param path2 The second path
#' @return The (admixture) probability of seeing the two paths together with the positive
#'         and negative edges in the overlap.
#' @export
path_overlap <- function(path1, path2) {
  rev_path2 <- data.frame(from = path2$to, to = path2$from, prob = path2$prob)
  list(prob = c(path_probability(path1), path_probability(path2)),
       positive = merge(path1, path2)[,1:2],
       negative = merge(path1, rev_path2)[,1:2])
}

#' Get the list of overlaps of all paths.
#' 
#' @param paths1 Paths between one pair of leaves
#' @param paths2 Paths between another pair of leaves
#' @return A list of the overlaps of all combinations of paths from \code{paths1} and \code{paths2}.
#' 
#' @export
all_path_overlaps <- function(paths1, paths2) {
  n1 <- length(paths1)
  n2 <- length(paths2)
  overlaps <- list(rep(NA), n1*n2)
  
  idx <- 1
  for (i in 1:n1) {
    for (j in 1:n2) {
      overlaps[[idx]] <- path_overlap(paths1[[i]], paths2[[j]])
      idx <- idx + 1
    }
  }
  overlaps
}


path_probability <- function(path) Filter(function(x) x != "", path$prob)
format_edge <- function(from, to) paste('[',from,':',to,']',sep='')

format_path_overlap <- function(overlap) {
  weight <- NULL
  if (length(overlap$prob) > 0) {
    weight <- paste(overlap$prob, collapse = " * ")
  }
  
  positive <- unlist(Map(format_edge, overlap$positive$from, overlap$positive$to), use.names = FALSE)
  negative <- unlist(Map(format_edge, overlap$negative$from, overlap$negative$to), use.names = FALSE)
  
  format_list <- c()
  
  if (length(positive) > 0) {
    format_list <- c(paste(positive, collapse = ' + '))
  }
  if (length(negative) > 0) {
    format_list <- c(format_list, ' - ', paste(negative, collapse = ' - '))
  }
  
  if (length(format_list) > 0 && !is.null(weight)) {
    format_list <- c(weight, ' * ', '(', format_list, ')')
  } else {
    format_list <- c('0')
  }
  paste(format_list, collapse = "")
}

format_all_overlaps <- function(overlaps) {
  overlaps <- vapply(overlaps, format_path_overlap, character(1))
  result <- paste(Filter(function(x) x != "0", overlaps), collapse = " + ")
  if (length(result) > 0) result else '0'
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
  overlaps <- all_path_overlaps(WX, YZ)
  format_all_overlaps(overlaps)
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

