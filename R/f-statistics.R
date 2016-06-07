#' Calculate the f_4(W, X; Y, Z) statistics.
#' 
#' Calculate the \eqn{f_4(W, X; Y, Z)} statistics.
#' 
#' @param graph  The admixture graph.
#' @param W      A leaf node.
#' @param X      A leaf node.
#' @param Y      A leaf node.
#' @param Z      A leaf node.
#'   
#' @return The overlaps between paths from \code{W} to \code{X} and paths from
#'         \code{Y} to \code{Z}.
#'         
#' @export
f4 <- function(graph, W, X, Y, Z) {
  if (!(W %in% graph$leaves)) stop(paste(W, "is not a leaf in the graph."))
  if (!(X %in% graph$leaves)) stop(paste(X, "is not a leaf in the graph."))
  if (!(Y %in% graph$leaves)) stop(paste(Y, "is not a leaf in the graph."))
  if (!(Z %in% graph$leaves)) stop(paste(Z, "is not a leaf in the graph."))
  
  WX <- all_paths(graph, W, X)
  YZ <- all_paths(graph, Y, Z)
  all_path_overlaps(WX, YZ)
}

#' Calculate the f_3(A; B, C) statistics.
#' 
#' Calculate the \eqn{f_3(A; B, C)} statistics.
#' 
#' @param graph  The admixture graph.
#' @param A      A leaf node.
#' @param B      A leaf node.
#' @param C      A leaf node.
#'   
#' @return A symbolic representation of the equation for the \eqn{f_3}
#'         statistics given by the admixture graph.
#'   
#' @export
f3 <- function(graph, A, B, C) f4(graph, A, B, A, C)

#' Calculate the f_2(A, B) statistics.
#' 
#' Calculate the \eqn{f_2(A, B)} statistics.
#' 
#' @param graph  The admixture graph.
#' @param A      A leaf node.
#' @param B      A leaf node.
#'   
#' @return A symbolic representation of the equation for the \eqn{f_2}
#'         statistics given by the admixture graph.
#'         
#' @export
f2 <- function(graph, A, B) f4(graph, A, B, A, B)