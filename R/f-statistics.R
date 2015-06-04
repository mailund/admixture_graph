#' Calculate the f4(W,X;Y,Z) statistics.
#'
#' @param graph The admixture graph
#' @param W     A leaf node
#' @param X     A leaf node
#' @param Y     A leaf node
#' @param Z     A leaf node
#'
#' @return The overlaps between paths from W to X and paths from Y to Z.
#' @export
f4 <- function(graph, W, X, Y, Z) {
  WX <- all_paths(graph, W, X)
  YZ <- all_paths(graph, Y, Z)
  all_path_overlaps(WX, YZ)
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
