path_zero <- function(overlap)
  nrow(overlap$positive) == 0 && nrow(overlap$negative) == 0
path_positive <- function(overlap)
  nrow(overlap$positive) > 0 && nrow(overlap$negative) == 0
path_negative <- function(overlap)
  nrow(overlap$positive) == 0 && nrow(overlap$negative) > 0

path_non_negative <- function(overlap)
  path_zero(overlap) || path_positive(overlap)
path_non_positive <- function(overlap)
  path_zero(overlap) || path_negative(overlap)

#' All overlaps are empty.
#' 
#' All overlaps are empty.
#' 
#' @param overlaps  Data frame representing path overlaps, typically generated
#'                  by \code{\link{all_path_overlaps}}.
#'   
#' @export
is_zero <- function(overlaps) {
  all(unlist(Map(path_zero, overlaps)))
}

#' All overlaps are either empty or have a positive weight.
#' 
#' All overlaps are either empty or have a positive weight.
#' 
#' @param overlaps  Data frame representing path overlaps, typically generated
#'                  by \code{\link{all_path_overlaps}}.
#'
#' @export
is_positive <- function(overlaps) {
  !is_zero(overlaps) && all(unlist(Map(path_non_negative, overlaps)))
}

#' All overlaps are either empty or have a negative weight.
#' 
#' All overlaps are either empty or have a negative weight.
#' 
#' @param overlaps  Data frame representing path overlaps, typically generated
#'                  by \code{\link{all_path_overlaps}}.
#'
#' @export
is_negative <- function(overlaps) {
  !is_zero(overlaps) && all(unlist(Map(path_non_positive, overlaps)))
}

#' Overlapping edges have both positive and negative contributions.
#' 
#' Overlapping edges have both positive and negative contributions.
#' 
#' @param overlaps  Data frame representing path overlaps, typically generated
#'                  by \code{\link{all_path_overlaps}}.
#'
#' @export
is_unknown <- function(overlaps) {
  !is_zero(overlaps) && !is_positive(overlaps) && !is_negative(overlaps)
}

#' Get the sign of overlapping paths.
#' 
#' Get the sign of overlapping paths.
#' 
#' @param overlaps  Data frame representing path overlaps, typically generated
#'                  by \code{\link{all_path_overlaps}}.
#'
#' @export
overlaps_sign <- function(overlaps) {
  if (is_zero(overlaps)) return(0)
  if (is_negative(overlaps)) return(-1)
  if (is_positive(overlaps)) return(+1)
  else return(NA)
}

#' Extracts the sign for the f_4 statistics predicted by the graph.
#' 
#' Extracts the sign for the \eqn{f_4} statistics predicted by the graph.
#' 
#' @param graph  The admixture graph.
#' @param W      First population/sample.
#' @param X      Second population/sample.
#' @param Y      Third population/sample.
#' @param Z      Fourth population/sample.
#'   
#' @return The sign of the \eqn{f_4} specified by the graph (or \code{NA} when it
#'         cannot be determined without knowing the graph parameters).
#'   
#' @export
get_graph_f4_sign <- function(graph, W, X, Y, Z) {
  overlaps_sign(f4(graph, W, X, Y, Z))
}

#' Extend a data frame with f_4 statistics predicted by a graph.
#' 
#' Extracts the sign for the \eqn{f_4} statistics predicted by the graph for all
#' rows in a data frame and extends the data frame with the graph \eqn{f_4}.
#' 
#' The data frame, \code{data}, must contain columns \code{W}, \code{X}, 
#' \code{Y}, and \code{Z}. The function then computes the sign of the 
#' \eqn{f4(W, X; Y, Z)} statistics for all rows and adds these as a column, 
#' \code{graph_f4_sign}, to the data frame.
#' 
#' @param data   The data frame to get the labels to compute the \eqn{f_4} 
#'               statistics from.
#' @param graph  The admixture graph.
#'   
#' @return A data frame identical to \code{data} except with an additional 
#'         column, \code{graph_f4_sign}, containing the sign of the \eqn{f_4}
#'         statistics as determined by the graph.
#'   
#' @export
add_graph_f4_sign <- function(data, graph) {
  signs <- unlist(Map(function(W,X,Y,Z) get_graph_f4_sign(graph, W, X, Y, Z),
                      data$W, data$X, data$Y, data$Z),
                  use.names = FALSE)
  data$graph_f4_sign <- signs
  data
}