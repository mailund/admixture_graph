#' Filter data so all W, X, Y and Z are leaves in the graph.
#' 
#' Filter data so all \code{W}, \code{X}, \code{Y} and \code{Z} are leaves in the graph.
#' 
#' @param data   Data frame (or similar) object containing columns \code{W}, \code{X},
#'               \code{Y}, and \code{Z}.
#' @param graph  Admixture graph.
#' 
#' @return Data frame with rows where \code{W}, \code{X}, \code{Y}, or \code{Z} are
#'         not leaves are removed.
#'   
#' @export
filter_on_leaves <- function(data, graph) {
  keep <- data$W %in% graph$leaves &
    data$X %in% graph$leaves &
    data$Y %in% graph$leaves & 
    data$Z %in% graph$leaves
  data[keep, ]
}