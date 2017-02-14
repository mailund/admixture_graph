path_probability <- function(path) Filter(function(x) x != "", path$prob)

format_edge <- function(graph) {
  force(graph)
  function(from, to) {
    if (graph$children[from,to]) {
      paste("edge_", from, "_", to, sep="")
    } else {
      paste("edge_", to, "_", from, sep="")
    }
  }
}

format_path_overlap <- function(graph) function(overlap) {
  weight <- NULL
  if (length(overlap$prob) > 0) {
    weight <- paste(overlap$prob, collapse = " * ")
  }

  positive <- unlist(Map(format_edge(graph),
                         overlap$positive$from, overlap$positive$to),
                     use.names = FALSE)
  negative <- unlist(Map(format_edge(graph),
                         overlap$negative$from, overlap$negative$to),
                     use.names = FALSE)

  format_list <- c()

  if (length(positive) > 0) {
    format_list <- c(paste(positive, collapse = " + "))
  }
  if (length(negative) > 0) {
    format_list <- c(format_list, " - ", paste(negative, collapse = " - "))
  }

  if (length(format_list) > 0 && !is.null(weight)) {
    format_list <- c(weight, " * ", "(", format_list, ")")
  } else if (length(format_list) == 0) {
    format_list <- c("0")
  }

  paste(format_list, collapse = "")
}

format_overlaps <- function(graph, overlaps) {
  overlaps <- vapply(overlaps, format_path_overlap(graph), character(1))
  result <- paste(Filter(function(x) x != "0", overlaps), collapse = " + ")
  if (result != "") parse(text = result) else expression(0)
}

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
#' @return A symbolic representation of the equation for the \eqn{f_4}
#'         statistics given by the admixture graph.
#'         
#' @export
sf4 <- function(graph, W, X, Y, Z) format_overlaps(graph, f4(graph, W, X, Y, Z))

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
sf3 <- function(graph, A, B, C) sf4(graph, A, B, A, C)

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
sf2 <- function(graph, A, B) sf4(graph, A, B, A, B)