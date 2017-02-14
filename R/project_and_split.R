#' Map sample names to population names.
#' 
#' Map sample names to population names.
#' 
#' This function maps the sample names in \code{W}, \code{X}, \code{Y}, and
#' \code{Z} to population names (typically what an admixture graph has for
#' leaves) and stores the original sample names so we can map them back again
#' after using the graph for making predictions.
#' 
#' @param data  The data frame to modify.
#' @param f     Function mapping sample names to population names.
#' 
#' @export
project_to_population <- function(data, f) {
  data$.W_ind <- data$W
  data$.X_ind <- data$X
  data$.Y_ind <- data$Y
  data$.Z_ind <- data$Z
  data$W <- vapply(data$.W_ind, f, "", USE.NAMES = FALSE)
  data$X <- vapply(data$.X_ind, f, "", USE.NAMES = FALSE)
  data$Y <- vapply(data$.Y_ind, f, "", USE.NAMES = FALSE)
  data$Z <- vapply(data$.Z_ind, f, "", USE.NAMES = FALSE)
  data
}

#' Reverse a projection of samples to populations.
#' 
#' Reverse a projection of samples to populations.
#' 
#' @param x  The projected data or a fitted object on projected data.
#' 
#' @export
split_population <- function(x) UseMethod("split_population")

#' Reverse a projection of samples to populations.
#' 
#' Reverse a projection of samples to populations.
#' 
#' @param x  The projected data or a fitted object on projected data.
#' 
#' @export
split_population.data.frame <- function(x) {
  x$W = x$.W_ind
  x$X = x$.X_ind
  x$Y = x$.Y_ind
  x$Z = x$.Z_ind
  x$.W_ind <- NULL
  x$.X_ind <- NULL
  x$.Y_ind <- NULL
  x$.Z_ind <- NULL
  x
}

#' Reverse a projection of samples to populations.
#' 
#' Reverse a projection of samples to populations.
#' 
#' @param x  The projected data or a fitted object on projected data.
#' 
#' @export
split_population.agraph_fit <- function(x) {
  x$data <- split_population(stats::fitted(x))
  x
}