#' Collect the postive and negative overlap between two paths.
#' 
#' Collects the postive and negative overlap between two paths.
#' 
#' @param path1  The first path.
#' @param path2  The second path.
#' 
#' @return The (admixture) probability of seeing the two paths together with the
#'         positive and negative edges in the overlap.
#'   
#' @export
path_overlap <- function(path1, path2) {
  rev_path2 <- data.frame(from = path2$to, to = path2$from, prob = path2$prob)
  list(prob = c(path_probability(path1), path_probability(path2)),
       positive = merge(path1, path2)[,1:2],
       negative = merge(path1, rev_path2)[,1:2])  
}

#' Get the list of overlaps of all paths.
#' 
#' Gets the list of overlaps of all paths.
#' 
#' @param paths1  Paths between one pair of leaves.
#' @param paths2  Paths between another pair of leaves.
#' @return A list of the overlaps of all combinations of paths from \code{paths1}
#'         and \code{paths2}.
#'   
#' @export
all_path_overlaps <- function(paths1, paths2) {
  n1 <- length(paths1)
  n2 <- length(paths2)
  
  if (n1 == 0 || n2 == 0) return(list())
  
  overlaps <- list(rep(NA, n1 * n2))

  idx <- 1
  for (i in 1:n1) {
    for (j in 1:n2) {
      overlaps[[idx]] <- path_overlap(paths1[[i]], paths2[[j]])
      idx <- idx + 1
    }
  }
  overlaps
}