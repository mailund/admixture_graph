
path_zero <- function(overlap) nrow(overlap$positive) == 0 && nrow(overlap$negative) == 0
path_positive <- function(overlap) nrow(overlap$positive) > 0 && nrow(overlap$negative) == 0
path_negative <- function(overlap) nrow(overlap$positive) == 0 && nrow(overlap$negative) > 0

path_non_negative <- function(overlap) path_zero(overlap) || path_positive(overlap)
path_non_positive <- function(overlap) path_zero(overlap) || path_negative(overlap)

#' All overlaps are empty
#' @export
is_zero <- function(overlaps) {
  all(unlist(Map(path_zero, overlaps)))
}

#' All overlaps are either empty or have a positive weight
#' @export
is_positive <- function(overlaps) {
  !is_zero(overlaps) && all(unlist(Map(path_non_negative, overlaps)))
}

#' All overlaps are either empty or have a negative weight
#' @export
is_negative <- function(overlaps) {
  !is_zero(overlaps) && all(unlist(Map(path_non_positive, overlaps)))
}

#' Overlapping edges have both positive and negative contributions
#' @export
is_unknown <- function(overlaps) {
  !is_zero(overlaps) && !is_positive(overlaps) && !is_negative(overlaps)
}

#' Get the sign of overlapping paths
#' @export
overlaps_sign <- function(overlaps) {
  if (is_zero(overlaps)) return('0')
  if (is_negative(overlaps)) return('-')
  if (is_positive(overlaps)) return('+')
  else return('?')
}