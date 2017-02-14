#' List of permutations.
#' 
#' List of permutations of given elements.
#' 
#' @param populations  A vector (of populations for example) of length between 4 and 8.
#' 
#' @return A list of different permutations of the elements of \code{x}.
#'
#' @seealso \code{\link{four_leaves_graphs}}
#' @seealso \code{\link{five_leaves_graphs}}
#' @seealso \code{\link{six_leaves_graphs}}
#' @seealso \code{\link{seven_leaves_trees}}
#' @seealso \code{\link{eight_leaves_trees}}
#' @seealso \code{\link{fit_permutations_and_graphs}}
#' @seealso \code{\link{add_a_leaf}}
#' @seealso \code{\link{add_an_admixture}}
#' @seealso \code{\link{add_an_admixture2}}
#'
#' @examples
#' \donttest{
#' # The number of permutations of n elements is n!. Take 0! = 1, 1! = 1, 2! = 2
#' # and 3! = 6 for granted. Now we can estimate e:
#' FOUR <- length(make_permutations(c(1, 2, 3, 4)))
#' FIVE <- length(make_permutations(c(1, 2, 3, 4, 5)))
#' SIX <- length(make_permutations(c(1, 2, 3, 4, 5, 6)))
#' SEVEN <- length(make_permutations(c(1, 2, 3, 4, 5, 6, 7)))
#' EIGHT <- length(make_permutations(c(1, 2, 3, 4, 5, 6, 7, 8)))
#' 1/1 + 1/1 + 1/2 + 1/6 + 1/FOUR + 1/FIVE + 1/SIX + 1/SEVEN + 1/EIGHT
#' # Hey that was pretty close!
#' }
#'
#' @export
make_permutations <- function(populations) {
  if (length(populations) == 1) {
    P <- list(populations[1])
  } else {
    P <- list()
    for (j in seq(1, length(populations))) {
      endings <- make_permutations(populations[-j])
      for (end in endings) {
        P[[length(P) + 1]] <- c(populations[j], end)
      }
    }
  }
  return(P)
}

#' Four leaves graphs.
#' 
#' Kind of obsolete since the introduction of \code{\link{all_graphs}}.
#' A comprehensive listing of all the \eqn{37} admixture graphs with four leaves and
#' at most two admixture events. Our convention is that the position of the root does
#' not matter (as long as it's not after an admixture event) and that graphs that have
#' \emph{eyes}, two inner nodes with the property that all the paths between any two
#' leaves visits both or neither of them, are excluded. The reason is that the \eqn{f}
#' statistics  can't detect the exact position of the root or distinguish between an
#' eye and a simple branch. The position of the root can be moved later with the function
#' \code{\link{make_an_outgroup}}.
#' 
#' @format A list of functions on four leaves and a parameter \code{permutations} which
#'         is \code{FALSE} by default.
#'         The outputs of these functions are either single \code{\link{agraph}} objects
#'         with the input vector as leaves, or if \code{permutations} is \code{TRUE},
#'         lists of all the possible \code{\link{agraph}} objects with that leaf set up
#'         to symmetry.
#' 
#' @family graphs
#' 
#' @seealso \code{\link{all_graphs}}
#' @seealso \code{\link{make_permutations}}
#' @seealso \code{\link{fit_permutations_and_graphs}}
#' @seealso \code{\link{fit_graph_list}}
#' @seealso \code{\link{add_a_leaf}}
#' @seealso \code{\link{add_an_admixture}}
#' @seealso \code{\link{add_an_admixture2}}
#' @seealso \code{\link{make_an_outgroup}}
#'
#' @examples
#' \donttest{
#' # While the usage of this function is pretty self-explanatory, let's plot all the graphs
#' # just for browsing.
#' for (i in seq(1, length(four_leaves_graphs))) {
#'   graph <- four_leaves_graphs[[i]](c("A", "B", "C", "D"))
#'   # This is how you include quotation marks in strings by the way:
#'   title <- paste("four_leaves_graphs[[", i, "]](c(\"A\", \"B\", \"C\", \"D\"))", sep = "")
#'   plot(graph, color = "tomato3", title = title)
#' }
#' }
#'
#' @export
four_leaves_graphs <- list(
  # single tree
  tree_1 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 3
      leaf_permutations <- symmetry_4_VI(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge(leaves[1], "R"),
        edge(leaves[2], "x"),
        edge(leaves[3], "y"),
        edge(leaves[4], "y")
      ))
      admixtures <- NULL
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },

  # one admixture event
  one_admixture_1 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 12
      leaf_permutations <- symmetry_4_VII(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "x"),
        edge(leaves[1], "R"),
        edge(leaves[2], "y"), 
        edge(leaves[3], "M"),
        edge(leaves[4], "z"),
        admixture_edge("M", "y", "z")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "y", "z", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  one_admixture_2 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 12
      leaf_permutations <- symmetry_4_IV(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "y", "z", "w", "M")
      edges <- parent_edges(c(
        edge("y", "R"),
        edge("z", "R"),
        edge("w", "z"),
        edge(leaves[1], "y"), 
        edge(leaves[2], "M"), 
        edge(leaves[3], "w"),
        edge(leaves[4], "w"),
        admixture_edge("M", "y", "z")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "y", "z", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  one_admixture_3 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 6
      leaf_permutations <- symmetry_4_III(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "y", "z", "w", "M")
      edges <- parent_edges(c(
        edge("y", "R"),
        edge("z", "R"),
        edge("w", "M"),
        edge(leaves[1], "y"), 
        edge(leaves[2], "w"), 
        edge(leaves[3], "w"),
        edge(leaves[4], "z"),
        admixture_edge("M", "y", "z")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "y", "z", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  # two admixture events
  two_admixtures_1 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 12
      leaf_permutations <- symmetry_4_II(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "w", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "y"),
        edge("w", "z"),
        edge(leaves[1], "x"), 
        edge(leaves[2], "M"), 
        edge(leaves[3], "N"),
        edge(leaves[4], "w"),
        admixture_edge("M", "x", "y"),
        admixture_edge("N", "z", "w")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "y", "a"),
        admix_props("N", "z", "w", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_2 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 12
      leaf_permutations <- symmetry_4_II(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "w", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "M"),
        edge("w", "z"),
        edge(leaves[1], "x"), 
        edge(leaves[2], "w"), 
        edge(leaves[3], "N"),
        edge(leaves[4], "y"),
        admixture_edge("M", "x", "y"),
        admixture_edge("N", "z", "w")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "y", "a"),
        admix_props("N", "z", "w", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_3 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 24
      leaf_permutations <- symmetry_4_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "w", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "x"),
        edge("w", "M"),
        edge(leaves[1], "R"), 
        edge(leaves[2], "N"), 
        edge(leaves[3], "w"),
        edge(leaves[4], "z"),
        admixture_edge("M", "y", "z"),
        admixture_edge("N", "y", "w")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "y", "z", "a"),
        admix_props("N", "y", "w", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_4 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 12
      leaf_permutations <- symmetry_4_VII(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "w", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "y"),
        edge("w", "z"),
        edge(leaves[1], "x"), 
        edge(leaves[2], "M"), 
        edge(leaves[3], "w"),
        edge(leaves[4], "N"),
        admixture_edge("M", "x", "z"),
        admixture_edge("N", "w", "y")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "z", "a"),
        admix_props("N", "w", "y", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_5 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 24
      leaf_permutations <- symmetry_4_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "w", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "x"),
        edge("w", "z"),
        edge(leaves[1], "R"), 
        edge(leaves[2], "y"), 
        edge(leaves[3], "N"),
        edge(leaves[4], "w"),
        admixture_edge("M", "y", "z"),
        admixture_edge("N", "M", "w")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "y", "z", "a"),
        admix_props("N", "M", "w", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_6 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 12
      leaf_permutations <- symmetry_4_II(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "w", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "y"),
        edge("w", "y"),
        edge(leaves[1], "x"), 
        edge(leaves[2], "M"), 
        edge(leaves[3], "N"),
        edge(leaves[4], "w"),
        admixture_edge("M", "x", "z"),
        admixture_edge("N", "z", "w")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "z", "a"),
        admix_props("N", "z", "w", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_7 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 24
      leaf_permutations <- symmetry_4_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "w", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "x"),
        edge("w", "z"),
        edge(leaves[1], "x"), 
        edge(leaves[2], "M"), 
        edge(leaves[3], "N"),
        edge(leaves[4], "y"),
        admixture_edge("M", "z", "w"),
        admixture_edge("N", "w", "y")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "z", "w", "a"),
        admix_props("N", "w", "y", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_8 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 24
      leaf_permutations <- symmetry_4_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "w", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "y"),
        edge("w", "z"),
        edge(leaves[1], "x"), 
        edge(leaves[2], "M"), 
        edge(leaves[3], "w"),
        edge(leaves[4], "N"),
        admixture_edge("M", "x", "w"),
        admixture_edge("N", "z", "y")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "w", "a"),
        admix_props("N", "z", "y", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_9 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 24
      leaf_permutations <- symmetry_4_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "w", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "x"),
        edge("w", "z"),
        edge(leaves[1], "R"), 
        edge(leaves[2], "y"), 
        edge(leaves[3], "N"),
        edge(leaves[4], "z"),
        admixture_edge("M", "y", "w"),
        admixture_edge("N", "M", "w")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "y", "w", "a"),
        admix_props("N", "M", "w", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_10 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 24
      leaf_permutations <- symmetry_4_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "w", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "y"),
        edge("w", "z"),
        edge(leaves[1], "x"), 
        edge(leaves[2], "N"), 
        edge(leaves[3], "w"),
        edge(leaves[4], "z"),
        admixture_edge("M", "x", "y"),
        admixture_edge("N", "M", "w")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "y", "a"),
        admix_props("N", "M", "w", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_11 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 24
      leaf_permutations <- symmetry_4_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "w", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "M"),
        edge("w", "y"),
        edge(leaves[1], "x"), 
        edge(leaves[2], "z"), 
        edge(leaves[3], "N"),
        edge(leaves[4], "w"),
        admixture_edge("M", "x", "y"),
        admixture_edge("N", "z", "w")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "y", "a"),
        admix_props("N", "z", "w", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_12 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 24
      leaf_permutations <- symmetry_4_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "w", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "M"),
        edge("w", "z"),
        edge(leaves[1], "R"), 
        edge(leaves[2], "z"), 
        edge(leaves[3], "w"),
        edge(leaves[4], "N"),
        admixture_edge("M", "x", "y"),
        admixture_edge("N", "w", "y")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "y", "a"),
        admix_props("N", "w", "y", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_13 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 6
      leaf_permutations <- symmetry_4_III(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "w", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "x"),
        edge("w", "y"),
        edge(leaves[1], "x"), 
        edge(leaves[2], "M"), 
        edge(leaves[3], "N"),
        edge(leaves[4], "y"),
        admixture_edge("M", "z", "w"),
        admixture_edge("N", "z", "w")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "z", "w", "a"),
        admix_props("N", "z", "w", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_14 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 24
      leaf_permutations <- symmetry_4_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "w", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "x"),
        edge("w", "z"),
        edge(leaves[1], "R"), 
        edge(leaves[2], "M"), 
        edge(leaves[3], "N"),
        edge(leaves[4], "w"),
        admixture_edge("M", "y", "z"),
        admixture_edge("N", "y", "w")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "y", "z", "a"),
        admix_props("N", "y", "w", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },

  two_admixtures_15 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 24
      leaf_permutations <- symmetry_4_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "w", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "x"),
        edge("w", "y"),
        edge(leaves[1], "R"), 
        edge(leaves[2], "w"), 
        edge(leaves[3], "M"),
        edge(leaves[4], "z"),
        admixture_edge("M", "y", "N"),
        admixture_edge("N", "w", "z")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "y", "N", "a"),
        admix_props("N", "w", "z", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_16 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 12
      leaf_permutations <- symmetry_4_V(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "w", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "y"),
        edge("w", "y"),
        edge(leaves[1], "R"), 
        edge(leaves[2], "z"), 
        edge(leaves[3], "w"),
        edge(leaves[4], "N"),
        admixture_edge("M", "z", "w"),
        admixture_edge("N", "M", "x")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "z", "w", "a"),
        admix_props("N", "M", "x", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_17 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 12
      leaf_permutations <- symmetry_4_II(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "w", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "y"),
        edge("w", "M"),
        edge(leaves[1], "x"), 
        edge(leaves[2], "w"), 
        edge(leaves[3], "N"),
        edge(leaves[4], "z"),
        admixture_edge("M", "x", "z"),
        admixture_edge("N", "w", "y")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "z", "a"),
        admix_props("N", "w", "y", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_18 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 6
      leaf_permutations <- symmetry_4_III(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "w", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "y"),
        edge("w", "y"),
        edge(leaves[1], "R"), 
        edge(leaves[2], "M"), 
        edge(leaves[3], "N"),
        edge(leaves[4], "x"),
        admixture_edge("M", "z", "w"),
        admixture_edge("N", "z", "w")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "z", "w", "a"),
        admix_props("N", "z", "w", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_19 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 12
      leaf_permutations <- symmetry_4_II(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "w", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "y"),
        edge("w", "z"),
        edge(leaves[1], "R"), 
        edge(leaves[2], "w"), 
        edge(leaves[3], "N"),
        edge(leaves[4], "x"),
        admixture_edge("M", "w", "y"),
        admixture_edge("N", "z", "M")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "w", "y", "a"),
        admix_props("N", "z", "M", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_20 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 12
      leaf_permutations <- symmetry_4_IV(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "w", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "x"),
        edge("w", "N"),
        edge(leaves[1], "R"), 
        edge(leaves[2], "M"), 
        edge(leaves[3], "w"),
        edge(leaves[4], "w"),
        admixture_edge("M", "y", "z"),
        admixture_edge("N", "y", "z")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "y", "z", "a"),
        admix_props("N", "y", "z", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_21 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 6
      leaf_permutations <- symmetry_4_VIII(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "w", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "y"),
        edge("w", "N"),
        edge(leaves[1], "R"), 
        edge(leaves[2], "z"), 
        edge(leaves[3], "w"),
        edge(leaves[4], "w"),
        admixture_edge("M", "z", "x"),
        admixture_edge("N", "y", "M")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "z", "x", "a"),
        admix_props("N", "y", "M", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_22 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 12
      leaf_permutations <- symmetry_4_V(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "w", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "y"),
        edge("w", "N"),
        edge(leaves[1], "x"), 
        edge(leaves[2], "w"), 
        edge(leaves[3], "w"),
        edge(leaves[4], "y"),
        admixture_edge("M", "x", "z"),
        admixture_edge("N", "M", "z")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "z", "a"),
        admix_props("N", "M", "z", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_23 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 12
      leaf_permutations <- symmetry_4_IV(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "w", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "y"),
        edge("w", "N"),
        edge(leaves[1], "x"), 
        edge(leaves[2], "M"), 
        edge(leaves[3], "w"),
        edge(leaves[4], "w"),
        admixture_edge("M", "x", "z"),
        admixture_edge("N", "y", "z")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "z", "a"),
        admix_props("N", "y", "z", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_24 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 12
      leaf_permutations <- symmetry_4_II(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "w", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "y"),
        edge("w", "z"),
        edge(leaves[1], "R"), 
        edge(leaves[2], "N"), 
        edge(leaves[3], "w"),
        edge(leaves[4], "x"),
        admixture_edge("M", "y", "z"),
        admixture_edge("N", "M", "w")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "y", "z", "a"),
        admix_props("N", "M", "w", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_25 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 12
      leaf_permutations <- symmetry_4_II(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "w", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "y"),
        edge("w", "y"),
        edge(leaves[1], "R"), 
        edge(leaves[2], "z"), 
        edge(leaves[3], "N"),
        edge(leaves[4], "x"),
        admixture_edge("M", "z", "w"),
        admixture_edge("N", "M", "w")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "z", "w", "a"),
        admix_props("N", "M", "w", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_26 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 12
      leaf_permutations <- symmetry_4_V(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "w", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "y"),
        edge("w", "M"),
        edge(leaves[1], "x"), 
        edge(leaves[2], "w"), 
        edge(leaves[3], "w"),
        edge(leaves[4], "N"),
        admixture_edge("M", "x", "z"),
        admixture_edge("N", "z", "y")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "z", "a"),
        admix_props("N", "z", "y", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_27 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 12
      leaf_permutations <- symmetry_4_II(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "w", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "y"),
        edge("w", "z"),
        edge(leaves[1], "R"), 
        edge(leaves[2], "M"), 
        edge(leaves[3], "N"),
        edge(leaves[4], "x"),
        admixture_edge("M", "y", "w"),
        admixture_edge("N", "w", "z")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "y", "w", "a"),
        admix_props("N", "w", "z", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_28 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 12
      leaf_permutations <- symmetry_4_V(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "w", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "y"),
        edge("w", "N"),
        edge(leaves[1], "x"), 
        edge(leaves[2], "w"), 
        edge(leaves[3], "w"),
        edge(leaves[4], "z"),
        admixture_edge("M", "x", "y"),
        admixture_edge("N", "M", "z")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "y", "a"),
        admix_props("N", "M", "z", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_29 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 12
      leaf_permutations <- symmetry_4_IV(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "w", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "M"),
        edge("w", "z"),
        edge(leaves[1], "x"), 
        edge(leaves[2], "N"), 
        edge(leaves[3], "w"),
        edge(leaves[4], "w"),
        admixture_edge("M", "x", "y"),
        admixture_edge("N", "y", "z")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "y", "a"),
        admix_props("N", "y", "z", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_30 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 12
      leaf_permutations <- symmetry_4_IV(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "w", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "M"),
        edge("w", "N"),
        edge(leaves[1], "x"), 
        edge(leaves[2], "z"), 
        edge(leaves[3], "w"),
        edge(leaves[4], "w"),
        admixture_edge("M", "x", "y"),
        admixture_edge("N", "z", "y")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "y", "a"),
        admix_props("N", "z", "y", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_31 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 12
      leaf_permutations <- symmetry_4_II(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "w", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "y"),
        edge("w", "z"),
        edge(leaves[1], "R"), 
        edge(leaves[2], "w"), 
        edge(leaves[3], "N"),
        edge(leaves[4], "x"),
        admixture_edge("M", "w", "z"),
        admixture_edge("N", "M", "y")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "w", "z", "a"),
        admix_props("N", "M", "y", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_32 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 12
      leaf_permutations <- symmetry_4_II(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "w", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "y"),
        edge("w", "x"),
        edge(leaves[1], "w"), 
        edge(leaves[2], "N"), 
        edge(leaves[3], "y"),
        edge(leaves[4], "w"),
        admixture_edge("M", "x", "z"),
        admixture_edge("N", "M", "z")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "z", "a"),
        admix_props("N", "M", "z", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_33 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 12
      leaf_permutations <- symmetry_4_II(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "w", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "M"),
        edge("w", "x"),
        edge(leaves[1], "w"), 
        edge(leaves[2], "N"), 
        edge(leaves[3], "z"),
        edge(leaves[4], "w"),
        admixture_edge("M", "x", "y"),
        admixture_edge("N", "y", "z")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "y", "a"),
        admix_props("N", "y", "z", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  }
)

symmetry_4_I <- function(leaves) {
  output <- character(4)
  all_permutations <- make_permutations(leaves)
  for (candidate in all_permutations) {
    output <- rbind(output, candidate)
  }
  output <- output[-1, ]
  output
}

symmetry_4_II <- function(leaves) {
  output <- character(4)
  all_permutations <- make_permutations(leaves)
  for (candidate in all_permutations) {
    sorted <- sort(candidate)
    num <- numeric(0)
    for (j in seq(1, 4)) {
      num[j] <- match(candidate[j], sorted)
    }
    if (num[1] < num[4]) {
      output <- rbind(output, candidate)
    }
  }
  output <- output[-1, ]
  output
}

symmetry_4_III <- function(leaves) {
  output <- character(4)
  all_permutations <- make_permutations(leaves)
  for (candidate in all_permutations) {
    sorted <- sort(candidate)
    num <- numeric(0)
    for (j in seq(1, 4)) {
      num[j] <- match(candidate[j], sorted)
    }
    if (num[1] < num[4] &&
        num[2] < num[3]) {
      output <- rbind(output, candidate)
    }
  }
  output <- output[-1, ]
  output
}

symmetry_4_IV <- function(leaves) {
  output <- character(4)
  all_permutations <- make_permutations(leaves)
  for (candidate in all_permutations) {
    sorted <- sort(candidate)
    num <- numeric(0)
    for (j in seq(1, 4)) {
      num[j] <- match(candidate[j], sorted)
    }
    if (num[3] < num[4]) {
      output <- rbind(output, candidate)
    }
  }
  output <- output[-1, ]
  output
}

symmetry_4_V <- function(leaves) {
  output <- character(4)
  all_permutations <- make_permutations(leaves)
  for (candidate in all_permutations) {
    sorted <- sort(candidate)
    num <- numeric(0)
    for (j in seq(1, 4)) {
      num[j] <- match(candidate[j], sorted)
    }
    if (num[2] < num[3]) {
      output <- rbind(output, candidate)
    }
  }
  output <- output[-1, ]
  output
}

symmetry_4_VI <- function(leaves) {
  output <- character(4)
  all_permutations <- make_permutations(leaves)
  for (candidate in all_permutations) {
    sorted <- sort(candidate)
    num <- numeric(0)
    for (j in seq(1, 4)) {
      num[j] <- match(candidate[j], sorted)
    }
    if (num[1] < num[2] &&
        num[3] < num[4] &&
        num[1] + 10*num[2] < num[3] + 10*num[4]) {
      output <- rbind(output, candidate)
    }
  }
  output <- output[-1, ]
  output
}

symmetry_4_VII <- function(leaves) {
  output <- character(4)
  all_permutations <- make_permutations(leaves)
  for (candidate in all_permutations) {
    sorted <- sort(candidate)
    num <- numeric(0)
    for (j in seq(1, 4)) {
      num[j] <- match(candidate[j], sorted)
    }
    if (num[2] < num[4]) {
      output <- rbind(output, candidate)
    }
  }
  output <- output[-1, ]
  output
}

symmetry_4_VIII <- function(leaves) {
  output <- character(4)
  all_permutations <- make_permutations(leaves)
  for (candidate in all_permutations) {
    sorted <- sort(candidate)
    num <- numeric(0)
    for (j in seq(1, 4)) {
      num[j] <- match(candidate[j], sorted)
    }
    if (num[1] < num[2] &&
        num[3] < num[4]) {
      output <- rbind(output, candidate)
    }
  }
  output <- output[-1, ]
  output
}

symmetry_4_IX <- function(leaves) {
  output <- character(4)
  all_permutations <- make_permutations(leaves)
  for (candidate in all_permutations) {
    sorted <- sort(candidate)
    num <- numeric(0)
    for (j in seq(1, 4)) {
      num[j] <- match(candidate[j], sorted)
    }
    if (num[1] < num[3]) {
      output <- rbind(output, candidate)
    }
  }
  output <- output[-1, ]
  output
}

#' Five leaves graphs.
#' 
#' Kind of obsolete since the introduction of \code{\link{all_graphs}}.
#' A comprehensive listing of all the \eqn{132} admixture graphs with five leaves and
#' at most two admixture events. Our convention is that the position of the root does
#' not matter (as long as it's not after an admixture event) and that graphs that have
#' \emph{eyes}, two inner nodes with the property that all the paths between any two
#' leaves visits both or neither of them, are excluded. The reason is that the \eqn{f}
#' statistics can't detect the exact position of the root or distinguish between an eye
#' and a simple branch. The position of the root can be moved later with the function
#' \code{\link{make_an_outgroup}}.
#' 
#' @format A list of functions on five leaves and a parameter \code{permutations} which
#'         is \code{FALSE} by default.
#'         The outputs of these functions are either single \code{\link{agraph}} objects
#'         with the input vector as leaves, or if \code{permutations} is \code{TRUE},
#'         lists of all the possible \code{\link{agraph}} objects with that leaf set up
#'         to symmetry.
#'         
#' @family graphs
#' 
#' @seealso \code{\link{all_graphs}}
#' @seealso \code{\link{make_permutations}}
#' @seealso \code{\link{fit_permutations_and_graphs}}
#' @seealso \code{\link{fit_graph_list}}
#' @seealso \code{\link{add_a_leaf}}
#' @seealso \code{\link{add_an_admixture}}
#' @seealso \code{\link{add_an_admixture2}}
#' @seealso \code{\link{make_an_outgroup}}
#'
#' @examples
#' \donttest{
#' # While the usage of this function is pretty self-explanatory, let's plot all the graphs
#' # just for browsing.
#' for (i in seq(1, length(five_leaves_graphs))) {
#'   graph <- five_leaves_graphs[[i]](c("A", "B", "C", "D", "E"))
#'   # This is how you include quotation marks in strings by the way:
#'   title <- paste("five_leaves_graphs[[", i, "]](c(\"A\", \"B\", \"C\", \"D\", \"E\"))",
#'                  sep = "")
#'   plot(graph, color = "purple", title = title)
#' }
#' }
#'
#' @export
five_leaves_graphs <- list(
  # tree
  tree_1 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 15
      leaf_permutations <- symmetry_5_IV(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "y"),
        edge(leaves[1], "R"),
        edge(leaves[2], "x"),
        edge(leaves[3], "y"),
        edge(leaves[4], "z"),
        edge(leaves[5], "z")
      ))
      admixtures <- NULL
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  # one admixture
  one_admixture_1 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "w", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "x"),
        edge("w", "z"),
        edge(leaves[1], "R"),
        edge(leaves[2], "y"), 
        edge(leaves[3], "M"),
        edge(leaves[4], "w"),
        edge(leaves[5], "z"),
        admixture_edge("M", "y", "w")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "y", "w", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  one_admixture_2 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 30
      leaf_permutations <- symmetry_5_II(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "w", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "x"),
        edge("w", "M"),
        edge(leaves[1], "R"),
        edge(leaves[2], "y"), 
        edge(leaves[3], "w"),
        edge(leaves[4], "w"),
        edge(leaves[5], "z"),
        admixture_edge("M", "y", "z")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "y", "z", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  one_admixture_3 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_III(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "w", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "x"),
        edge("w", "z"),
        edge(leaves[1], "R"),
        edge(leaves[2], "y"), 
        edge(leaves[3], "M"),
        edge(leaves[4], "w"),
        edge(leaves[5], "w"),
        admixture_edge("M", "y", "z")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "y", "z", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },

  one_admixture_4 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 15
      leaf_permutations <- symmetry_5_IV(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "w", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "x"),
        edge("w", "y"),
        edge(leaves[1], "z"),
        edge(leaves[2], "z"), 
        edge(leaves[3], "M"),
        edge(leaves[4], "w"),
        edge(leaves[5], "w"),
        admixture_edge("M", "x", "y")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "y", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },

  one_admixture_5 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 30
      leaf_permutations <- symmetry_5_V(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "w", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "M"),
        edge("w", "y"),
        edge(leaves[1], "x"),
        edge(leaves[2], "z"), 
        edge(leaves[3], "z"),
        edge(leaves[4], "w"),
        edge(leaves[5], "w"),
        admixture_edge("M", "x", "y")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "y", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  one_admixture_6 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 30
      leaf_permutations <- symmetry_5_VI(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "w", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "M"),
        edge("w", "z"),
        edge(leaves[1], "x"),
        edge(leaves[2], "z"), 
        edge(leaves[3], "w"),
        edge(leaves[4], "w"),
        edge(leaves[5], "y"),
        admixture_edge("M", "x", "y")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "y", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  one_admixture_7 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_III(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "w", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "y"),
        edge("w", "z"),
        edge(leaves[1], "x"),
        edge(leaves[2], "M"), 
        edge(leaves[3], "z"),
        edge(leaves[4], "w"),
        edge(leaves[5], "w"),
        admixture_edge("M", "x", "y")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "y", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  one_admixture_8 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 30
      leaf_permutations <- symmetry_5_II(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "w", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "x"),
        edge("w", "R"),
        edge(leaves[1], "M"),
        edge(leaves[2], "y"), 
        edge(leaves[3], "w"),
        edge(leaves[4], "w"),
        edge(leaves[5], "z"),
        admixture_edge("M", "y", "z")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "y", "z", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  # two admixtures
  two_admixtures_1 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "y"),
        edge("u", "M"),
        edge("v", "u"),
        edge(leaves[1], "x"),
        edge(leaves[2], "N"), 
        edge(leaves[3], "v"),
        edge(leaves[4], "y"),
        edge(leaves[5], "z"),
        admixture_edge("M", "x", "z"),
        admixture_edge("N", "u", "v")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "z", "a"),
        admix_props("N", "u", "v", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_2 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 120
      leaf_permutations <- symmetry_5_VII(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "y"),
        edge("u", "z"),
        edge("v", "u"),
        edge(leaves[1], "x"),
        edge(leaves[2], "M"), 
        edge(leaves[3], "N"),
        edge(leaves[4], "v"),
        edge(leaves[5], "u"),
        admixture_edge("M", "x", "y"),
        admixture_edge("N", "z", "v")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "y", "a"),
        admix_props("N", "z", "v", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_3 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "y"),
        edge("u", "z"),
        edge("v", "z"),
        edge(leaves[1], "u"),
        edge(leaves[2], "x"), 
        edge(leaves[3], "M"),
        edge(leaves[4], "N"),
        edge(leaves[5], "v"),
        admixture_edge("M", "x", "y"),
        admixture_edge("N", "u", "v")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "y", "a"),
        admix_props("N", "u", "v", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_4 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "M"),
        edge("u", "z"),
        edge("v", "u"),
        edge(leaves[1], "x"),
        edge(leaves[2], "N"), 
        edge(leaves[3], "v"),
        edge(leaves[4], "u"),
        edge(leaves[5], "y"),
        admixture_edge("M", "x", "y"),
        admixture_edge("N", "z", "v")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "y", "a"),
        admix_props("N", "z", "v", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_5 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 30
      leaf_permutations <- symmetry_5_II(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "M"),
        edge("u", "z"),
        edge("v", "z"),
        edge(leaves[1], "N"),
        edge(leaves[2], "u"), 
        edge(leaves[3], "x"),
        edge(leaves[4], "y"),
        edge(leaves[5], "v"),
        admixture_edge("M", "x", "y"),
        admixture_edge("N", "u", "v")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "y", "a"),
        admix_props("N", "u", "v", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_6 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "y"),
        edge("u", "M"),
        edge("v", "u"),
        edge(leaves[1], "z"),
        edge(leaves[2], "x"), 
        edge(leaves[3], "N"),
        edge(leaves[4], "v"),
        edge(leaves[5], "z"),
        admixture_edge("M", "x", "y"),
        admixture_edge("N", "u", "v")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "y", "a"),
        admix_props("N", "u", "v", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_7 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "x"),
        edge("u", "y"),
        edge("v", "u"),
        edge(leaves[1], "z"),
        edge(leaves[2], "M"), 
        edge(leaves[3], "N"),
        edge(leaves[4], "v"),
        edge(leaves[5], "z"),
        admixture_edge("M", "x", "y"),
        admixture_edge("N", "u", "v")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "y", "a"),
        admix_props("N", "u", "v", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_8 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "M"),
        edge("u", "y"),
        edge("v", "u"),
        edge(leaves[1], "z"),
        edge(leaves[2], "x"), 
        edge(leaves[3], "N"),
        edge(leaves[4], "v"),
        edge(leaves[5], "z"),
        admixture_edge("M", "x", "y"),
        admixture_edge("N", "u", "v")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "y", "a"),
        admix_props("N", "u", "v", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_9 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 30
      leaf_permutations <- symmetry_5_II(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "v"),
        edge("u", "M"),
        edge("v", "u"),
        edge(leaves[1], "N"),
        edge(leaves[2], "z"), 
        edge(leaves[3], "x"),
        edge(leaves[4], "y"),
        edge(leaves[5], "z"),
        admixture_edge("M", "x", "y"),
        admixture_edge("N", "u", "v")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "y", "a"),
        admix_props("N", "u", "v", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_10 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 30
      leaf_permutations <- symmetry_5_II(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "N"),
        edge("u", "M"),
        edge("v", "u"),
        edge(leaves[1], "v"),
        edge(leaves[2], "z"), 
        edge(leaves[3], "x"),
        edge(leaves[4], "y"),
        edge(leaves[5], "z"),
        admixture_edge("M", "x", "y"),
        admixture_edge("N", "u", "v")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "y", "a"),
        admix_props("N", "u", "v", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_11 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "z"),
        edge("y", "x"),
        edge("z", "R"),
        edge("u", "z"),
        edge("v", "u"),
        edge(leaves[1], "y"),
        edge(leaves[2], "M"), 
        edge(leaves[3], "R"),
        edge(leaves[4], "N"),
        edge(leaves[5], "v"),
        admixture_edge("M", "x", "y"),
        admixture_edge("N", "u", "v")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "y", "a"),
        admix_props("N", "u", "v", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_12 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "M"),
        edge("u", "z"),
        edge("v", "u"),
        edge(leaves[1], "x"),
        edge(leaves[2], "z"), 
        edge(leaves[3], "N"),
        edge(leaves[4], "v"),
        edge(leaves[5], "y"),
        admixture_edge("M", "x", "y"),
        admixture_edge("N", "u", "v")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "y", "a"),
        admix_props("N", "u", "v", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_13 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 120
      leaf_permutations <- symmetry_5_VII(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "x"),
        edge("u", "z"),
        edge("v", "y"),
        edge(leaves[1], "x"),
        edge(leaves[2], "u"), 
        edge(leaves[3], "N"),
        edge(leaves[4], "v"),
        edge(leaves[5], "y"),
        admixture_edge("M", "z", "v"),
        admixture_edge("N", "u", "M")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "z", "v", "a"),
        admix_props("N", "u", "M", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_14 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 120
      leaf_permutations <- symmetry_5_VII(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "y"),
        edge("u", "z"),
        edge("v", "M"),
        edge(leaves[1], "x"),
        edge(leaves[2], "v"), 
        edge(leaves[3], "N"),
        edge(leaves[4], "z"),
        edge(leaves[5], "y"),
        admixture_edge("M", "x", "u"),
        admixture_edge("N", "v", "u")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "u", "a"),
        admix_props("N", "v", "u", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_15 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 120
      leaf_permutations <- symmetry_5_VII(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "y"),
        edge("u", "z"),
        edge("v", "u"),
        edge(leaves[1], "x"),
        edge(leaves[2], "M"), 
        edge(leaves[3], "v"),
        edge(leaves[4], "N"),
        edge(leaves[5], "y"),
        admixture_edge("M", "x", "u"),
        admixture_edge("N", "v", "z")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "u", "a"),
        admix_props("N", "v", "z", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_16 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 120
      leaf_permutations <- symmetry_5_VII(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "x"),
        edge("u", "x"),
        edge("v", "u"),
        edge(leaves[1], "z"),
        edge(leaves[2], "M"), 
        edge(leaves[3], "v"),
        edge(leaves[4], "N"),
        edge(leaves[5], "y"),
        admixture_edge("M", "z", "u"),
        admixture_edge("N", "v", "y")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "z", "u", "a"),
        admix_props("N", "v", "y", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_17 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 120
      leaf_permutations <- symmetry_5_VII(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "x"),
        edge("u", "x"),
        edge("v", "y"),
        edge(leaves[1], "z"),
        edge(leaves[2], "M"), 
        edge(leaves[3], "N"),
        edge(leaves[4], "v"),
        edge(leaves[5], "y"),
        admixture_edge("M", "z", "u"),
        admixture_edge("N", "u", "v")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "z", "u", "a"),
        admix_props("N", "u", "v", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_18 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 120
      leaf_permutations <- symmetry_5_VII(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "x"),
        edge("u", "z"),
        edge("v", "u"),
        edge(leaves[1], "x"),
        edge(leaves[2], "u"), 
        edge(leaves[3], "v"),
        edge(leaves[4], "N"),
        edge(leaves[5], "y"),
        admixture_edge("M", "z", "y"),
        admixture_edge("N", "v", "M")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "z", "y", "a"),
        admix_props("N", "v", "M", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_19 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 120
      leaf_permutations <- symmetry_5_VII(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "x"),
        edge("u", "z"),
        edge("v", "M"),
        edge(leaves[1], "x"),
        edge(leaves[2], "u"), 
        edge(leaves[3], "N"),
        edge(leaves[4], "v"),
        edge(leaves[5], "y"),
        admixture_edge("M", "z", "y"),
        admixture_edge("N", "u", "v")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "z", "y", "a"),
        admix_props("N", "u", "v", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_20 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 120
      leaf_permutations <- symmetry_5_VII(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "x"),
        edge("u", "M"),
        edge("v", "u"),
        edge(leaves[1], "x"),
        edge(leaves[2], "N"), 
        edge(leaves[3], "v"),
        edge(leaves[4], "u"),
        edge(leaves[5], "y"),
        admixture_edge("M", "z", "y"),
        admixture_edge("N", "z", "v")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "z", "y", "a"),
        admix_props("N", "z", "v", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_21 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 120
      leaf_permutations <- symmetry_5_VII(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "x"),
        edge("u", "y"),
        edge("v", "u"),
        edge(leaves[1], "x"),
        edge(leaves[2], "z"), 
        edge(leaves[3], "N"),
        edge(leaves[4], "u"),
        edge(leaves[5], "y"),
        admixture_edge("M", "z", "v"),
        admixture_edge("N", "M", "v")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "z", "v", "a"),
        admix_props("N", "M", "v", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_22 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 120
      leaf_permutations <- symmetry_5_VII(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "x"),
        edge("u", "y"),
        edge("v", "u"),
        edge(leaves[1], "x"),
        edge(leaves[2], "z"), 
        edge(leaves[3], "M"),
        edge(leaves[4], "N"),
        edge(leaves[5], "y"),
        admixture_edge("M", "z", "v"),
        admixture_edge("N", "v", "u")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "z", "v", "a"),
        admix_props("N", "v", "u", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_23 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 120
      leaf_permutations <- symmetry_5_VII(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "v"),
        edge("u", "y"),
        edge("v", "u"),
        edge(leaves[1], "x"),
        edge(leaves[2], "M"), 
        edge(leaves[3], "z"),
        edge(leaves[4], "N"),
        edge(leaves[5], "y"),
        admixture_edge("M", "x", "z"),
        admixture_edge("N", "v", "u")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "z", "a"),
        admix_props("N", "v", "u", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_24 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 120
      leaf_permutations <- symmetry_5_VII(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "x"),
        edge("u", "z"),
        edge("v", "u"),
        edge(leaves[1], "z"),
        edge(leaves[2], "u"), 
        edge(leaves[3], "v"),
        edge(leaves[4], "N"),
        edge(leaves[5], "y"),
        admixture_edge("M", "x", "y"),
        admixture_edge("N", "v", "M")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "y", "a"),
        admix_props("N", "v", "M", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_25 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 120
      leaf_permutations <- symmetry_5_VII(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "x"),
        edge("u", "z"),
        edge("v", "M"),
        edge(leaves[1], "z"),
        edge(leaves[2], "u"), 
        edge(leaves[3], "N"),
        edge(leaves[4], "v"),
        edge(leaves[5], "y"),
        admixture_edge("M", "x", "y"),
        admixture_edge("N", "u", "v")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "y", "a"),
        admix_props("N", "u", "v", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_26 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 120
      leaf_permutations <- symmetry_5_VII(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "x"),
        edge("u", "M"),
        edge("v", "u"),
        edge(leaves[1], "z"),
        edge(leaves[2], "N"), 
        edge(leaves[3], "v"),
        edge(leaves[4], "u"),
        edge(leaves[5], "y"),
        admixture_edge("M", "x", "y"),
        admixture_edge("N", "z", "v")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "y", "a"),
        admix_props("N", "z", "v", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_27 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 120
      leaf_permutations <- symmetry_5_VII(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "M"),
        edge("u", "z"),
        edge("v", "u"),
        edge(leaves[1], "N"),
        edge(leaves[2], "v"), 
        edge(leaves[3], "u"),
        edge(leaves[4], "z"),
        edge(leaves[5], "y"),
        admixture_edge("M", "x", "y"),
        admixture_edge("N", "x", "v")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "y", "a"),
        admix_props("N", "x", "v", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_28 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "y"),
        edge("u", "y"),
        edge("v", "u"),
        edge(leaves[1], "R"),
        edge(leaves[2], "z"), 
        edge(leaves[3], "N"),
        edge(leaves[4], "v"),
        edge(leaves[5], "x"),
        admixture_edge("M", "z", "u"),
        admixture_edge("N", "M", "v")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "z", "u", "a"),
        admix_props("N", "M", "v", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_29 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "y"),
        edge("u", "z"),
        edge("v", "M"),
        edge(leaves[1], "R"),
        edge(leaves[2], "z"), 
        edge(leaves[3], "N"),
        edge(leaves[4], "v"),
        edge(leaves[5], "x"),
        admixture_edge("M", "u", "y"),
        admixture_edge("N", "u", "v")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "u", "y", "a"),
        admix_props("N", "u", "v", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_30 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "y"),
        edge("u", "z"),
        edge("v", "z"),
        edge(leaves[1], "R"),
        edge(leaves[2], "v"), 
        edge(leaves[3], "N"),
        edge(leaves[4], "M"),
        edge(leaves[5], "x"),
        admixture_edge("M", "u", "y"),
        admixture_edge("N", "u", "v")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "u", "y", "a"),
        admix_props("N", "u", "v", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_31 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 120
      leaf_permutations <- symmetry_5_VII(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "x"),
        edge("u", "z"),
        edge("v", "M"),
        edge(leaves[1], "v"),
        edge(leaves[2], "u"), 
        edge(leaves[3], "N"),
        edge(leaves[4], "y"),
        edge(leaves[5], "v"),
        admixture_edge("M", "R", "u"),
        admixture_edge("N", "z", "y")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "R", "u", "a"),
        admix_props("N", "z", "y", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_32 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 120
      leaf_permutations <- symmetry_5_VII(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "x"),
        edge("u", "z"),
        edge("v", "u"),
        edge(leaves[1], "v"),
        edge(leaves[2], "x"), 
        edge(leaves[3], "N"),
        edge(leaves[4], "y"),
        edge(leaves[5], "v"),
        admixture_edge("M", "z", "y"),
        admixture_edge("N", "u", "M")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "z", "y", "a"),
        admix_props("N", "u", "M", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_33 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 120
      leaf_permutations <- symmetry_5_VII(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "x"),
        edge("u", "M"),
        edge("v", "N"),
        edge(leaves[1], "v"),
        edge(leaves[2], "x"), 
        edge(leaves[3], "u"),
        edge(leaves[4], "y"),
        edge(leaves[5], "v"),
        admixture_edge("M", "z", "y"),
        admixture_edge("N", "z", "u")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "z", "y", "a"),
        admix_props("N", "z", "u", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_34 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "y"),
        edge("u", "y"),
        edge("v", "M"),
        edge(leaves[1], "R"),
        edge(leaves[2], "N"), 
        edge(leaves[3], "v"),
        edge(leaves[4], "u"),
        edge(leaves[5], "x"),
        admixture_edge("M", "z", "u"),
        admixture_edge("N", "z", "v")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "z", "u", "a"),
        admix_props("N", "z", "v", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_35 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "y"),
        edge("u", "z"),
        edge("v", "u"),
        edge(leaves[1], "R"),
        edge(leaves[2], "z"), 
        edge(leaves[3], "v"),
        edge(leaves[4], "N"),
        edge(leaves[5], "x"),
        admixture_edge("M", "u", "y"),
        admixture_edge("N", "v", "M")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "u", "y", "a"),
        admix_props("N", "v", "M", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_36 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "y"),
        edge("u", "z"),
        edge("v", "u"),
        edge(leaves[1], "R"),
        edge(leaves[2], "N"), 
        edge(leaves[3], "v"),
        edge(leaves[4], "M"),
        edge(leaves[5], "x"),
        admixture_edge("M", "u", "y"),
        admixture_edge("N", "z", "v")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "u", "y", "a"),
        admix_props("N", "z", "v", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_37 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "y"),
        edge("u", "y"),
        edge("v", "N"),
        edge(leaves[1], "v"),
        edge(leaves[2], "x"), 
        edge(leaves[3], "M"),
        edge(leaves[4], "z"),
        edge(leaves[5], "v"),
        admixture_edge("M", "x", "u"),
        admixture_edge("N", "u", "z")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "u", "a"),
        admix_props("N", "u", "z", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_38 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "y"),
        edge("u", "M"),
        edge("v", "u"),
        edge(leaves[1], "v"),
        edge(leaves[2], "x"), 
        edge(leaves[3], "N"),
        edge(leaves[4], "y"),
        edge(leaves[5], "v"),
        admixture_edge("M", "x", "z"),
        admixture_edge("N", "u", "z")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "z", "a"),
        admix_props("N", "u", "z", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_39 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "y"),
        edge("u", "z"),
        edge("v", "N"),
        edge(leaves[1], "v"),
        edge(leaves[2], "x"), 
        edge(leaves[3], "u"),
        edge(leaves[4], "y"),
        edge(leaves[5], "v"),
        admixture_edge("M", "x", "z"),
        admixture_edge("N", "u", "M")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "z", "a"),
        admix_props("N", "u", "M", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_40 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "y"),
        edge("u", "x"),
        edge("v", "N"),
        edge(leaves[1], "v"),
        edge(leaves[2], "x"), 
        edge(leaves[3], "z"),
        edge(leaves[4], "y"),
        edge(leaves[5], "v"),
        admixture_edge("M", "u", "z"),
        admixture_edge("N", "u", "M")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "u", "z", "a"),
        admix_props("N", "u", "M", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_41 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "x"),
        edge("u", "z"),
        edge("v", "N"),
        edge(leaves[1], "v"),
        edge(leaves[2], "z"), 
        edge(leaves[3], "u"),
        edge(leaves[4], "M"),
        edge(leaves[5], "v"),
        admixture_edge("M", "u", "y"),
        admixture_edge("N", "y", "R")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "u", "y", "a"),
        admix_props("N", "y", "R", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_42 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 30
      leaf_permutations <- symmetry_5_II(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "M"),
        edge("u", "y"),
        edge("v", "R"),
        edge(leaves[1], "N"),
        edge(leaves[2], "z"), 
        edge(leaves[3], "u"),
        edge(leaves[4], "v"),
        edge(leaves[5], "z"),
        admixture_edge("M", "x", "y"),
        admixture_edge("N", "u", "v")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "y", "a"),
        admix_props("N", "u", "v", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_43 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "y"),
        edge("u", "z"),
        edge("v", "u"),
        edge(leaves[1], "x"),
        edge(leaves[2], "u"), 
        edge(leaves[3], "v"),
        edge(leaves[4], "N"),
        edge(leaves[5], "x"),
        admixture_edge("M", "z", "y"),
        admixture_edge("N", "M", "v")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "z", "y", "a"),
        admix_props("N", "M", "v", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_44 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "y"),
        edge("u", "z"),
        edge("v", "M"),
        edge(leaves[1], "x"),
        edge(leaves[2], "u"), 
        edge(leaves[3], "N"),
        edge(leaves[4], "v"),
        edge(leaves[5], "x"),
        admixture_edge("M", "z", "y"),
        admixture_edge("N", "u", "v")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "z", "y", "a"),
        admix_props("N", "u", "v", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_45 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "y"),
        edge("u", "M"),
        edge("v", "u"),
        edge(leaves[1], "x"),
        edge(leaves[2], "N"), 
        edge(leaves[3], "v"),
        edge(leaves[4], "u"),
        edge(leaves[5], "x"),
        admixture_edge("M", "z", "y"),
        admixture_edge("N", "z", "v")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "z", "y", "a"),
        admix_props("N", "z", "v", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_46 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "y"),
        edge("u", "z"),
        edge("v", "u"),
        edge(leaves[1], "x"),
        edge(leaves[2], "N"), 
        edge(leaves[3], "u"),
        edge(leaves[4], "z"),
        edge(leaves[5], "x"),
        admixture_edge("M", "v", "y"),
        admixture_edge("N", "M", "v")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "v", "y", "a"),
        admix_props("N", "M", "v", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_47 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "y"),
        edge("u", "z"),
        edge("v", "y"),
        edge(leaves[1], "x"),
        edge(leaves[2], "N"), 
        edge(leaves[3], "u"),
        edge(leaves[4], "z"),
        edge(leaves[5], "x"),
        admixture_edge("M", "v", "u"),
        admixture_edge("N", "M", "v")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "v", "u", "a"),
        admix_props("N", "M", "v", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_48 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "M"),
        edge("u", "x"),
        edge("v", "u"),
        edge(leaves[1], "z"),
        edge(leaves[2], "x"), 
        edge(leaves[3], "N"),
        edge(leaves[4], "y"),
        edge(leaves[5], "z"),
        admixture_edge("M", "v", "y"),
        admixture_edge("N", "u", "v")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "v", "y", "a"),
        admix_props("N", "u", "v", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_49 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "y"),
        edge("u", "y"),
        edge("v", "u"),
        edge(leaves[1], "R"),
        edge(leaves[2], "N"), 
        edge(leaves[3], "M"),
        edge(leaves[4], "z"),
        edge(leaves[5], "x"),
        admixture_edge("M", "v", "z"),
        admixture_edge("N", "u", "v")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "v", "z", "a"),
        admix_props("N", "u", "v", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_50 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "v"),
        edge("u", "y"),
        edge("v", "u"),
        edge(leaves[1], "R"),
        edge(leaves[2], "N"), 
        edge(leaves[3], "z"),
        edge(leaves[4], "M"),
        edge(leaves[5], "x"),
        admixture_edge("M", "y", "z"),
        admixture_edge("N", "u", "v")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "y", "z", "a"),
        admix_props("N", "u", "v", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_51 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "x"),
        edge("u", "z"),
        edge("v", "N"),
        edge(leaves[1], "v"),
        edge(leaves[2], "z"), 
        edge(leaves[3], "u"),
        edge(leaves[4], "y"),
        edge(leaves[5], "v"),
        admixture_edge("M", "y", "x"),
        admixture_edge("N", "u", "M")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "y", "x", "a"),
        admix_props("N", "u", "M", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_52 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "x"),
        edge("u", "M"),
        edge("v", "u"),
        edge(leaves[1], "v"),
        edge(leaves[2], "z"), 
        edge(leaves[3], "N"),
        edge(leaves[4], "y"),
        edge(leaves[5], "v"),
        admixture_edge("M", "y", "x"),
        admixture_edge("N", "u", "z")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "y", "x", "a"),
        admix_props("N", "u", "z", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_53 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "M"),
        edge("u", "z"),
        edge("v", "z"),
        edge(leaves[1], "v"),
        edge(leaves[2], "u"), 
        edge(leaves[3], "N"),
        edge(leaves[4], "y"),
        edge(leaves[5], "v"),
        admixture_edge("M", "y", "x"),
        admixture_edge("N", "u", "x")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "y", "x", "a"),
        admix_props("N", "u", "x", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_54 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "M"),
        edge("u", "z"),
        edge("v", "N"),
        edge(leaves[1], "v"),
        edge(leaves[2], "u"), 
        edge(leaves[3], "z"),
        edge(leaves[4], "y"),
        edge(leaves[5], "v"),
        admixture_edge("M", "y", "x"),
        admixture_edge("N", "u", "x")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "y", "x", "a"),
        admix_props("N", "u", "x", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_55 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "M"),
        edge("u", "y"),
        edge("v", "u"),
        edge(leaves[1], "v"),
        edge(leaves[2], "x"), 
        edge(leaves[3], "z"),
        edge(leaves[4], "N"),
        edge(leaves[5], "v"),
        admixture_edge("M", "y", "x"),
        admixture_edge("N", "u", "z")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "y", "x", "a"),
        admix_props("N", "u", "z", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_56 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "y"),
        edge("u", "z"),
        edge("v", "z"),
        edge(leaves[1], "v"),
        edge(leaves[2], "x"), 
        edge(leaves[3], "u"),
        edge(leaves[4], "N"),
        edge(leaves[5], "v"),
        admixture_edge("M", "y", "x"),
        admixture_edge("N", "u", "M")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "y", "x", "a"),
        admix_props("N", "u", "M", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_57 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "y"),
        edge("u", "z"),
        edge("v", "y"),
        edge(leaves[1], "v"),
        edge(leaves[2], "x"), 
        edge(leaves[3], "z"),
        edge(leaves[4], "N"),
        edge(leaves[5], "v"),
        admixture_edge("M", "u", "x"),
        admixture_edge("N", "u", "M")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "u", "x", "a"),
        admix_props("N", "u", "M", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_58 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "y"),
        edge("u", "z"),
        edge("v", "x"),
        edge(leaves[1], "v"),
        edge(leaves[2], "y"), 
        edge(leaves[3], "M"),
        edge(leaves[4], "N"),
        edge(leaves[5], "v"),
        admixture_edge("M", "u", "x"),
        admixture_edge("N", "u", "z")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "u", "x", "a"),
        admix_props("N", "u", "z", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_59 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 30
      leaf_permutations <- symmetry_5_II(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "x"),
        edge("u", "z"),
        edge("v", "N"),
        edge(leaves[1], "M"),
        edge(leaves[2], "u"), 
        edge(leaves[3], "v"),
        edge(leaves[4], "v"),
        edge(leaves[5], "y"),
        admixture_edge("M", "z", "x"),
        admixture_edge("N", "u", "y")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "z", "x", "a"),
        admix_props("N", "u", "y", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_60 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "x"),
        edge("u", "z"),
        edge("v", "u"),
        edge(leaves[1], "v"),
        edge(leaves[2], "z"), 
        edge(leaves[3], "N"),
        edge(leaves[4], "y"),
        edge(leaves[5], "v"),
        admixture_edge("M", "y", "x"),
        admixture_edge("N", "u", "M")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "y", "x", "a"),
        admix_props("N", "u", "M", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_61 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "x"),
        edge("u", "M"),
        edge("v", "N"),
        edge(leaves[1], "v"),
        edge(leaves[2], "z"), 
        edge(leaves[3], "u"),
        edge(leaves[4], "y"),
        edge(leaves[5], "v"),
        admixture_edge("M", "y", "x"),
        admixture_edge("N", "u", "z")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "y", "x", "a"),
        admix_props("N", "u", "z", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_62 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "M"),
        edge("u", "z"),
        edge("v", "u"),
        edge(leaves[1], "v"),
        edge(leaves[2], "x"), 
        edge(leaves[3], "N"),
        edge(leaves[4], "z"),
        edge(leaves[5], "v"),
        admixture_edge("M", "y", "R"),
        admixture_edge("N", "u", "y")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "y", "R", "a"),
        admix_props("N", "u", "y", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_63 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 30
      leaf_permutations <- symmetry_5_II(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "y"),
        edge("u", "N"),
        edge("v", "z"),
        edge(leaves[1], "x"),
        edge(leaves[2], "u"), 
        edge(leaves[3], "v"),
        edge(leaves[4], "v"),
        edge(leaves[5], "u"),
        admixture_edge("M", "y", "x"),
        admixture_edge("N", "M", "z")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "y", "x", "a"),
        admix_props("N", "M", "z", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_64 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 30
      leaf_permutations <- symmetry_5_II(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "M"),
        edge("u", "x"),
        edge("v", "u"),
        edge(leaves[1], "N"),
        edge(leaves[2], "z"), 
        edge(leaves[3], "y"),
        edge(leaves[4], "y"),
        edge(leaves[5], "z"),
        admixture_edge("M", "v", "R"),
        admixture_edge("N", "v", "u")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "v", "R", "a"),
        admix_props("N", "v", "u", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_65 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 30
      leaf_permutations <- symmetry_5_II(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "M"),
        edge("u", "z"),
        edge("v", "N"),
        edge(leaves[1], "x"),
        edge(leaves[2], "v"), 
        edge(leaves[3], "u"),
        edge(leaves[4], "u"),
        edge(leaves[5], "v"),
        admixture_edge("M", "x", "y"),
        admixture_edge("N", "z", "y")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "y", "a"),
        admix_props("N", "z", "y", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_66 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 30
      leaf_permutations <- symmetry_5_II(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "x"),
        edge("u", "x"),
        edge("v", "y"),
        edge(leaves[1], "N"),
        edge(leaves[2], "v"), 
        edge(leaves[3], "u"),
        edge(leaves[4], "u"),
        edge(leaves[5], "v"),
        admixture_edge("M", "z", "y"),
        admixture_edge("N", "z", "M")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "z", "y", "a"),
        admix_props("N", "z", "M", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_67 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 30
      leaf_permutations <- symmetry_5_II(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "y"),
        edge("u", "x"),
        edge("v", "z"),
        edge(leaves[1], "N"),
        edge(leaves[2], "v"), 
        edge(leaves[3], "u"),
        edge(leaves[4], "u"),
        edge(leaves[5], "v"),
        admixture_edge("M", "x", "y"),
        admixture_edge("N", "z", "M")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "y", "a"),
        admix_props("N", "z", "M", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_68 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 30
      leaf_permutations <- symmetry_5_II(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "y"),
        edge("u", "z"),
        edge("v", "M"),
        edge(leaves[1], "N"),
        edge(leaves[2], "u"), 
        edge(leaves[3], "v"),
        edge(leaves[4], "v"),
        edge(leaves[5], "u"),
        admixture_edge("M", "y", "R"),
        admixture_edge("N", "x", "z")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "y", "R", "a"),
        admix_props("N", "x", "z", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_69 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 30
      leaf_permutations <- symmetry_5_II(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "y"),
        edge("u", "M"),
        edge("v", "u"),
        edge(leaves[1], "N"),
        edge(leaves[2], "v"), 
        edge(leaves[3], "R"),
        edge(leaves[4], "x"),
        edge(leaves[5], "v"),
        admixture_edge("M", "y", "z"),
        admixture_edge("N", "u", "z")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "y", "z", "a"),
        admix_props("N", "u", "z", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_70 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 30
      leaf_permutations <- symmetry_5_II(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "y"),
        edge("u", "M"),
        edge("v", "N"),
        edge(leaves[1], "u"),
        edge(leaves[2], "v"), 
        edge(leaves[3], "R"),
        edge(leaves[4], "x"),
        edge(leaves[5], "v"),
        admixture_edge("M", "y", "z"),
        admixture_edge("N", "u", "z")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "y", "z", "a"),
        admix_props("N", "u", "z", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_71 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 30
      leaf_permutations <- symmetry_5_II(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "y"),
        edge("u", "M"),
        edge("v", "N"),
        edge(leaves[1], "x"),
        edge(leaves[2], "u"), 
        edge(leaves[3], "v"),
        edge(leaves[4], "v"),
        edge(leaves[5], "u"),
        admixture_edge("M", "x", "z"),
        admixture_edge("N", "y", "z")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "z", "a"),
        admix_props("N", "y", "z", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_72 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 30
      leaf_permutations <- symmetry_5_II(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "x"),
        edge("u", "N"),
        edge("v", "y"),
        edge(leaves[1], "z"),
        edge(leaves[2], "u"), 
        edge(leaves[3], "v"),
        edge(leaves[4], "v"),
        edge(leaves[5], "u"),
        admixture_edge("M", "y", "x"),
        admixture_edge("N", "M", "z")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "y", "x", "a"),
        admix_props("N", "M", "z", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_73 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 30
      leaf_permutations <- symmetry_5_II(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "x"),
        edge("u", "x"),
        edge("v", "N"),
        edge(leaves[1], "y"),
        edge(leaves[2], "u"), 
        edge(leaves[3], "v"),
        edge(leaves[4], "v"),
        edge(leaves[5], "u"),
        admixture_edge("M", "y", "z"),
        admixture_edge("N", "M", "z")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "y", "z", "a"),
        admix_props("N", "M", "z", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_74 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 30
      leaf_permutations <- symmetry_5_II(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "y"),
        edge("u", "x"),
        edge("v", "N"),
        edge(leaves[1], "y"),
        edge(leaves[2], "u"), 
        edge(leaves[3], "v"),
        edge(leaves[4], "v"),
        edge(leaves[5], "u"),
        admixture_edge("M", "x", "z"),
        admixture_edge("N", "M", "z")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "z", "a"),
        admix_props("N", "M", "z", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_75 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "y"),
        edge("u", "N"),
        edge("v", "u"),
        edge(leaves[1], "v"),
        edge(leaves[2], "x"), 
        edge(leaves[3], "y"),
        edge(leaves[4], "u"),
        edge(leaves[5], "v"),
        admixture_edge("M", "x", "z"),
        admixture_edge("N", "M", "z")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "z", "a"),
        admix_props("N", "M", "z", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_76 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "y"),
        edge("u", "N"),
        edge("v", "u"),
        edge(leaves[1], "v"),
        edge(leaves[2], "x"), 
        edge(leaves[3], "M"),
        edge(leaves[4], "u"),
        edge(leaves[5], "v"),
        admixture_edge("M", "R", "z"),
        admixture_edge("N", "y", "z")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "R", "z", "a"),
        admix_props("N", "y", "z", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_77 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "y"),
        edge("u", "x"),
        edge("v", "u"),
        edge(leaves[1], "v"),
        edge(leaves[2], "u"), 
        edge(leaves[3], "N"),
        edge(leaves[4], "z"),
        edge(leaves[5], "v"),
        admixture_edge("M", "x", "y"),
        admixture_edge("N", "M", "z")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "y", "a"),
        admix_props("N", "M", "z", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_78 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "x"),
        edge("u", "x"),
        edge("v", "u"),
        edge(leaves[1], "v"),
        edge(leaves[2], "u"), 
        edge(leaves[3], "N"),
        edge(leaves[4], "y"),
        edge(leaves[5], "v"),
        admixture_edge("M", "z", "y"),
        admixture_edge("N", "M", "z")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "z", "y", "a"),
        admix_props("N", "M", "z", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_79 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "y"),
        edge("u", "x"),
        edge("v", "u"),
        edge(leaves[1], "v"),
        edge(leaves[2], "u"), 
        edge(leaves[3], "N"),
        edge(leaves[4], "y"),
        edge(leaves[5], "v"),
        admixture_edge("M", "z", "x"),
        admixture_edge("N", "M", "z")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "z", "x", "a"),
        admix_props("N", "M", "z", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_80 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "y"),
        edge("u", "x"),
        edge("v", "u"),
        edge(leaves[1], "v"),
        edge(leaves[2], "u"), 
        edge(leaves[3], "N"),
        edge(leaves[4], "M"),
        edge(leaves[5], "v"),
        admixture_edge("M", "z", "x"),
        admixture_edge("N", "y", "z")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "z", "x", "a"),
        admix_props("N", "y", "z", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_81 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "y"),
        edge("u", "M"),
        edge("v", "u"),
        edge(leaves[1], "v"),
        edge(leaves[2], "u"), 
        edge(leaves[3], "N"),
        edge(leaves[4], "x"),
        edge(leaves[5], "v"),
        admixture_edge("M", "z", "x"),
        admixture_edge("N", "y", "z")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "z", "x", "a"),
        admix_props("N", "y", "z", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_82 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "x"),
        edge("u", "z"),
        edge("v", "u"),
        edge(leaves[1], "v"),
        edge(leaves[2], "u"), 
        edge(leaves[3], "N"),
        edge(leaves[4], "y"),
        edge(leaves[5], "v"),
        admixture_edge("M", "y", "x"),
        admixture_edge("N", "M", "z")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "y", "x", "a"),
        admix_props("N", "M", "z", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_83 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "M"),
        edge("u", "N"),
        edge("v", "u"),
        edge(leaves[1], "v"),
        edge(leaves[2], "u"), 
        edge(leaves[3], "x"),
        edge(leaves[4], "z"),
        edge(leaves[5], "v"),
        admixture_edge("M", "y", "R"),
        admixture_edge("N", "y", "z")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "y", "R", "a"),
        admix_props("N", "y", "z", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_84 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "M"),
        edge("u", "z"),
        edge("v", "u"),
        edge(leaves[1], "v"),
        edge(leaves[2], "u"), 
        edge(leaves[3], "x"),
        edge(leaves[4], "N"),
        edge(leaves[5], "v"),
        admixture_edge("M", "y", "R"),
        admixture_edge("N", "y", "z")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "y", "R", "a"),
        admix_props("N", "y", "z", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_85 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "y"),
        edge("u", "N"),
        edge("v", "u"),
        edge(leaves[1], "v"),
        edge(leaves[2], "u"), 
        edge(leaves[3], "x"),
        edge(leaves[4], "z"),
        edge(leaves[5], "v"),
        admixture_edge("M", "y", "x"),
        admixture_edge("N", "M", "z")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "y", "x", "a"),
        admix_props("N", "M", "z", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_86 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 120
      leaf_permutations <- symmetry_5_VII(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "y"),
        edge("u", "z"),
        edge("v", "z"),
        edge(leaves[1], "x"),
        edge(leaves[2], "N"), 
        edge(leaves[3], "v"),
        edge(leaves[4], "u"),
        edge(leaves[5], "y"),
        admixture_edge("M", "u", "x"),
        admixture_edge("N", "M", "v")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "u", "x", "a"),
        admix_props("N", "M", "v", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_87 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 120
      leaf_permutations <- symmetry_5_VII(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "y"),
        edge("u", "z"),
        edge("v", "M"),
        edge(leaves[1], "x"),
        edge(leaves[2], "N"), 
        edge(leaves[3], "v"),
        edge(leaves[4], "u"),
        edge(leaves[5], "y"),
        admixture_edge("M", "u", "x"),
        admixture_edge("N", "z", "v")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "u", "x", "a"),
        admix_props("N", "z", "v", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_88 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "y"),
        edge("u", "z"),
        edge("v", "y"),
        edge(leaves[1], "x"),
        edge(leaves[2], "N"), 
        edge(leaves[3], "M"),
        edge(leaves[4], "z"),
        edge(leaves[5], "v"),
        admixture_edge("M", "u", "x"),
        admixture_edge("N", "u", "v")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "u", "x", "a"),
        admix_props("N", "u", "v", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_89 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "y"),
        edge("u", "z"),
        edge("v", "u"),
        edge(leaves[1], "x"),
        edge(leaves[2], "N"), 
        edge(leaves[3], "M"),
        edge(leaves[4], "z"),
        edge(leaves[5], "v"),
        admixture_edge("M", "u", "x"),
        admixture_edge("N", "y", "v")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "u", "x", "a"),
        admix_props("N", "y", "v", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_90 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 120
      leaf_permutations <- symmetry_5_VII(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "v"),
        edge("u", "x"),
        edge("v", "y"),
        edge(leaves[1], "x"),
        edge(leaves[2], "N"), 
        edge(leaves[3], "M"),
        edge(leaves[4], "z"),
        edge(leaves[5], "y"),
        admixture_edge("M", "u", "v"),
        admixture_edge("N", "z", "u")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "u", "v", "a"),
        admix_props("N", "z", "u", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_91 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "y"),
        edge("u", "x"),
        edge("v", "z"),
        edge(leaves[1], "x"),
        edge(leaves[2], "N"), 
        edge(leaves[3], "u"),
        edge(leaves[4], "v"),
        edge(leaves[5], "z"),
        admixture_edge("M", "u", "v"),
        admixture_edge("N", "y", "M")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "u", "v", "a"),
        admix_props("N", "y", "M", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_92 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 120
      leaf_permutations <- symmetry_5_VII(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "x"),
        edge("u", "y"),
        edge("v", "u"),
        edge(leaves[1], "x"),
        edge(leaves[2], "N"), 
        edge(leaves[3], "z"),
        edge(leaves[4], "v"),
        edge(leaves[5], "y"),
        admixture_edge("M", "z", "v"),
        admixture_edge("N", "u", "M")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "z", "v", "a"),
        admix_props("N", "u", "M", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_93 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 120
      leaf_permutations <- symmetry_5_VII(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "x"),
        edge("u", "z"),
        edge("v", "u"),
        edge(leaves[1], "x"),
        edge(leaves[2], "N"), 
        edge(leaves[3], "M"),
        edge(leaves[4], "u"),
        edge(leaves[5], "y"),
        admixture_edge("M", "z", "v"),
        admixture_edge("N", "v", "y")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "z", "v", "a"),
        admix_props("N", "v", "y", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_94 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "x"),
        edge("u", "y"),
        edge("v", "z"),
        edge(leaves[1], "u"),
        edge(leaves[2], "N"), 
        edge(leaves[3], "M"),
        edge(leaves[4], "R"),
        edge(leaves[5], "v"),
        admixture_edge("M", "u", "v"),
        admixture_edge("N", "z", "y")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "u", "v", "a"),
        admix_props("N", "z", "y", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_95 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 30
      leaf_permutations <- symmetry_5_II(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "x"),
        edge("u", "y"),
        edge("v", "z"),
        edge(leaves[1], "R"),
        edge(leaves[2], "M"), 
        edge(leaves[3], "y"),
        edge(leaves[4], "z"),
        edge(leaves[5], "N"),
        admixture_edge("M", "u", "v"),
        admixture_edge("N", "u", "v")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "u", "v", "a"),
        admix_props("N", "u", "v", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_96 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "y"),
        edge("u", "z"),
        edge("v", "z"),
        edge(leaves[1], "u"),
        edge(leaves[2], "x"), 
        edge(leaves[3], "N"),
        edge(leaves[4], "y"),
        edge(leaves[5], "v"),
        admixture_edge("M", "v", "u"),
        admixture_edge("N", "M", "x")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "v", "u", "a"),
        admix_props("N", "M", "x", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_97 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "M"),
        edge("u", "y"),
        edge("v", "y"),
        edge(leaves[1], "u"),
        edge(leaves[2], "x"), 
        edge(leaves[3], "N"),
        edge(leaves[4], "z"),
        edge(leaves[5], "v"),
        admixture_edge("M", "u", "v"),
        admixture_edge("N", "x", "z")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "u", "v", "a"),
        admix_props("N", "x", "z", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_98 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "M"),
        edge("z", "y"),
        edge("u", "x"),
        edge("v", "x"),
        edge(leaves[1], "u"),
        edge(leaves[2], "y"), 
        edge(leaves[3], "N"),
        edge(leaves[4], "z"),
        edge(leaves[5], "v"),
        admixture_edge("M", "u", "v"),
        admixture_edge("N", "R", "z")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "u", "v", "a"),
        admix_props("N", "R", "z", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_99 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "x"),
        edge("u", "z"),
        edge("v", "N"),
        edge(leaves[1], "v"),
        edge(leaves[2], "y"), 
        edge(leaves[3], "x"),
        edge(leaves[4], "u"),
        edge(leaves[5], "v"),
        admixture_edge("M", "u", "y"),
        admixture_edge("N", "M", "z")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "u", "y", "a"),
        admix_props("N", "M", "z", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_100 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 30
      leaf_permutations <- symmetry_5_II(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "x"),
        edge("u", "y"),
        edge("v", "N"),
        edge(leaves[1], "M"),
        edge(leaves[2], "v"), 
        edge(leaves[3], "x"),
        edge(leaves[4], "y"),
        edge(leaves[5], "v"),
        admixture_edge("M", "u", "z"),
        admixture_edge("N", "u", "z")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "u", "z", "a"),
        admix_props("N", "u", "z", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_101 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "y"),
        edge("u", "z"),
        edge("v", "z"),
        edge(leaves[1], "v"),
        edge(leaves[2], "M"), 
        edge(leaves[3], "x"),
        edge(leaves[4], "N"),
        edge(leaves[5], "v"),
        admixture_edge("M", "u", "x"),
        admixture_edge("N", "u", "y")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "u", "x", "a"),
        admix_props("N", "u", "y", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_102 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "y"),
        edge("u", "z"),
        edge("v", "N"),
        edge(leaves[1], "v"),
        edge(leaves[2], "M"), 
        edge(leaves[3], "x"),
        edge(leaves[4], "z"),
        edge(leaves[5], "v"),
        admixture_edge("M", "u", "x"),
        admixture_edge("N", "u", "y")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "u", "x", "a"),
        admix_props("N", "u", "y", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_103 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "y"),
        edge("u", "z"),
        edge("v", "u"),
        edge(leaves[1], "v"),
        edge(leaves[2], "N"), 
        edge(leaves[3], "x"),
        edge(leaves[4], "y"),
        edge(leaves[5], "v"),
        admixture_edge("M", "u", "x"),
        admixture_edge("N", "M", "z")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "u", "x", "a"),
        admix_props("N", "M", "z", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_104 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "x"),
        edge("u", "x"),
        edge("v", "z"),
        edge(leaves[1], "v"),
        edge(leaves[2], "N"), 
        edge(leaves[3], "u"),
        edge(leaves[4], "y"),
        edge(leaves[5], "v"),
        admixture_edge("M", "z", "y"),
        admixture_edge("N", "M", "u")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "z", "y", "a"),
        admix_props("N", "M", "u", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_105 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "x"),
        edge("u", "M"),
        edge("v", "z"),
        edge(leaves[1], "v"),
        edge(leaves[2], "N"), 
        edge(leaves[3], "u"),
        edge(leaves[4], "y"),
        edge(leaves[5], "v"),
        admixture_edge("M", "z", "y"),
        admixture_edge("N", "x", "u")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "z", "y", "a"),
        admix_props("N", "x", "u", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_106 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "y"),
        edge("u", "z"),
        edge("v", "x"),
        edge(leaves[1], "v"),
        edge(leaves[2], "N"), 
        edge(leaves[3], "u"),
        edge(leaves[4], "y"),
        edge(leaves[5], "v"),
        admixture_edge("M", "x", "u"),
        admixture_edge("N", "M", "z")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "u", "a"),
        admix_props("N", "M", "z", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_107 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 30
      leaf_permutations <- symmetry_5_II(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "x"),
        edge("u", "x"),
        edge("v", "y"),
        edge(leaves[1], "y"),
        edge(leaves[2], "M"), 
        edge(leaves[3], "z"),
        edge(leaves[4], "z"),
        edge(leaves[5], "N"),
        admixture_edge("M", "v", "u"),
        admixture_edge("N", "v", "u")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "v", "u", "a"),
        admix_props("N", "v", "u", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_108 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "y"),
        edge("u", "z"),
        edge("v", "y"),
        edge(leaves[1], "v"),
        edge(leaves[2], "x"), 
        edge(leaves[3], "u"),
        edge(leaves[4], "N"),
        edge(leaves[5], "v"),
        admixture_edge("M", "x", "u"),
        admixture_edge("N", "M", "z")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "u", "a"),
        admix_props("N", "M", "z", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_109 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "x"),
        edge("u", "y"),
        edge("v", "N"),
        edge(leaves[1], "v"),
        edge(leaves[2], "M"), 
        edge(leaves[3], "z"),
        edge(leaves[4], "y"),
        edge(leaves[5], "v"),
        admixture_edge("M", "x", "u"),
        admixture_edge("N", "u", "z")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "u", "a"),
        admix_props("N", "u", "z", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_110 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "x"),
        edge("u", "y"),
        edge("v", "z"),
        edge(leaves[1], "v"),
        edge(leaves[2], "M"), 
        edge(leaves[3], "N"),
        edge(leaves[4], "y"),
        edge(leaves[5], "v"),
        admixture_edge("M", "x", "u"),
        admixture_edge("N", "u", "z")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "u", "a"),
        admix_props("N", "u", "z", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_111 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 30
      leaf_permutations <- symmetry_5_II(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "x"),
        edge("u", "x"),
        edge("v", "N"),
        edge(leaves[1], "u"),
        edge(leaves[2], "z"), 
        edge(leaves[3], "v"),
        edge(leaves[4], "v"),
        edge(leaves[5], "y"),
        admixture_edge("M", "z", "y"),
        admixture_edge("N", "u", "M")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "z", "y", "a"),
        admix_props("N", "u", "M", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_112 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 30
      leaf_permutations <- symmetry_5_II(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "y"),
        edge("u", "M"),
        edge("v", "u"),
        edge(leaves[1], "N"),
        edge(leaves[2], "x"), 
        edge(leaves[3], "v"),
        edge(leaves[4], "v"),
        edge(leaves[5], "z"),
        admixture_edge("M", "z", "R"),
        admixture_edge("N", "u", "y")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "z", "R", "a"),
        admix_props("N", "u", "y", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_113 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 30
      leaf_permutations <- symmetry_5_II(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "y"),
        edge("u", "M"),
        edge("v", "N"),
        edge(leaves[1], "u"),
        edge(leaves[2], "x"), 
        edge(leaves[3], "v"),
        edge(leaves[4], "v"),
        edge(leaves[5], "z"),
        admixture_edge("M", "z", "R"),
        admixture_edge("N", "u", "y")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "z", "R", "a"),
        admix_props("N", "u", "y", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_114 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 30
      leaf_permutations <- symmetry_5_II(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "y"),
        edge("u", "y"),
        edge("v", "u"),
        edge(leaves[1], "N"),
        edge(leaves[2], "x"), 
        edge(leaves[3], "v"),
        edge(leaves[4], "v"),
        edge(leaves[5], "z"),
        admixture_edge("M", "z", "R"),
        admixture_edge("N", "u", "M")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "z", "R", "a"),
        admix_props("N", "u", "M", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_115 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 30
      leaf_permutations <- symmetry_5_II(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "y"),
        edge("u", "z"),
        edge("v", "N"),
        edge(leaves[1], "u"),
        edge(leaves[2], "v"), 
        edge(leaves[3], "R"),
        edge(leaves[4], "x"),
        edge(leaves[5], "v"),
        admixture_edge("M", "u", "y"),
        admixture_edge("N", "z", "M")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "u", "y", "a"),
        admix_props("N", "z", "M", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_116 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 30
      leaf_permutations <- symmetry_5_II(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "y"),
        edge("u", "y"),
        edge("v", "N"),
        edge(leaves[1], "M"),
        edge(leaves[2], "v"), 
        edge(leaves[3], "R"),
        edge(leaves[4], "x"),
        edge(leaves[5], "v"),
        admixture_edge("M", "z", "u"),
        admixture_edge("N", "z", "u")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "z", "u", "a"),
        admix_props("N", "z", "u", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_117 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 15
      leaf_permutations <- symmetry_5_IV(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "x"),
        edge("u", "M"),
        edge("v", "N"),
        edge(leaves[1], "u"),
        edge(leaves[2], "u"), 
        edge(leaves[3], "R"),
        edge(leaves[4], "v"),
        edge(leaves[5], "v"),
        admixture_edge("M", "z", "y"),
        admixture_edge("N", "z", "y")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "z", "y", "a"),
        admix_props("N", "z", "y", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_118 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 30
      leaf_permutations <- symmetry_5_II(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "x"),
        edge("u", "N"),
        edge("v", "u"),
        edge(leaves[1], "u"),
        edge(leaves[2], "y"), 
        edge(leaves[3], "v"),
        edge(leaves[4], "v"),
        edge(leaves[5], "z"),
        admixture_edge("M", "z", "y"),
        admixture_edge("N", "x", "M")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "z", "y", "a"),
        admix_props("N", "x", "M", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_119 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "x"),
        edge("u", "N"),
        edge("v", "u"),
        edge(leaves[1], "v"),
        edge(leaves[2], "R"), 
        edge(leaves[3], "M"),
        edge(leaves[4], "u"),
        edge(leaves[5], "v"),
        admixture_edge("M", "z", "y"),
        admixture_edge("N", "z", "y")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "z", "y", "a"),
        admix_props("N", "z", "y", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_120 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "y"),
        edge("u", "R"),
        edge("v", "u"),
        edge(leaves[1], "v"),
        edge(leaves[2], "z"), 
        edge(leaves[3], "N"),
        edge(leaves[4], "u"),
        edge(leaves[5], "v"),
        admixture_edge("M", "z", "x"),
        admixture_edge("N", "M", "y")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "z", "x", "a"),
        admix_props("N", "M", "y", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_121 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 30
      leaf_permutations <- symmetry_5_II(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "x"),
        edge("u", "R"),
        edge("v", "u"),
        edge(leaves[1], "u"),
        edge(leaves[2], "M"), 
        edge(leaves[3], "v"),
        edge(leaves[4], "v"),
        edge(leaves[5], "N"),
        admixture_edge("M", "z", "y"),
        admixture_edge("N", "z", "y")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "z", "y", "a"),
        admix_props("N", "z", "y", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_122 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 15
      leaf_permutations <- symmetry_5_IV(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "z"),
        edge("z", "R"),
        edge("u", "x"),
        edge("v", "y"),
        edge(leaves[1], "u"),
        edge(leaves[2], "u"), 
        edge(leaves[3], "N"),
        edge(leaves[4], "v"),
        edge(leaves[5], "v"),
        admixture_edge("M", "x", "y"),
        admixture_edge("N", "z", "M")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "y", "a"),
        admix_props("N", "z", "M", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  two_admixtures_123 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 60
      leaf_permutations <- symmetry_5_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "M", "N")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "x"),
        edge("u", "z"),
        edge("v", "M"),
        edge(leaves[1], "u"),
        edge(leaves[2], "z"), 
        edge(leaves[3], "N"),
        edge(leaves[4], "v"),
        edge(leaves[5], "u"),
        admixture_edge("M", "R", "y"),
        admixture_edge("N", "v", "y")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "R", "y", "a"),
        admix_props("N", "v", "y", "b")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  }
)

symmetry_5_I <- function(leaves) {
  output <- character(5)
  all_permutations <- make_permutations(leaves)
  for (candidate in all_permutations) {
    sorted <- sort(candidate)
    num <- numeric(0)
    for (j in seq(1, 5)) {
      num[j] <- match(candidate[j], sorted)
    }
    if (num[1] < num[5]) {
      output <- rbind(output, candidate)
    }
  }
  output <- output[-1, ]
  output
}

symmetry_5_II <- function(leaves) {
  output <- character(5)
  all_permutations <- make_permutations(leaves)
  for (candidate in all_permutations) {
    sorted <- sort(candidate)
    num <- numeric(0)
    for (j in seq(1, 5)) {
      num[j] <- match(candidate[j], sorted)
    }
    if (num[2] < num[5] &&
        num[3] < num[4]) {
      output <- rbind(output, candidate)
    }
  }
  output <- output[-1, ]
  output
}

symmetry_5_III <- function(leaves) {
  output <- character(5)
  all_permutations <- make_permutations(leaves)
  for (candidate in all_permutations) {
    sorted <- sort(candidate)
    num <- numeric(0)
    for (j in seq(1, 5)) {
      num[j] <- match(candidate[j], sorted)
    }
    if (num[4] < num[5]) {
      output <- rbind(output, candidate)
    }
  }
  output <- output[-1, ]
  output
}

symmetry_5_IV <- function(leaves) {
  output <- character(5)
  all_permutations <- make_permutations(leaves)
  for (candidate in all_permutations) {
    sorted <- sort(candidate)
    num <- numeric(0)
    for (j in seq(1, 5)) {
      num[j] <- match(candidate[j], sorted)
    }
    if (num[1] < num[2] &&
        num[4] < num[5] &&
        num[1] + 10*num[2] < num[4] + 10*num[5]) {
      output <- rbind(output, candidate)
    }
  }
  output <- output[-1, ]
  output
}

symmetry_5_V <- function(leaves) {
  output <- character(5)
  all_permutations <- make_permutations(leaves)
  for (candidate in all_permutations) {
    sorted <- sort(candidate)
    num <- numeric(0)
    for (j in seq(1, 5)) {
      num[j] <- match(candidate[j], sorted)
    }
    if (num[2] < num[3] &&
        num[4] < num[5]) {
      output <- rbind(output, candidate)
    }
  }
  output <- output[-1, ]
  output
}

symmetry_5_VI <- function(leaves) {
  output <- character(5)
  all_permutations <- make_permutations(leaves)
  for (candidate in all_permutations) {
    sorted <- sort(candidate)
    num <- numeric(0)
    for (j in seq(1, 5)) {
      num[j] <- match(candidate[j], sorted)
    }
    if (num[3] < num[4] &&
        num[1] < num[5]) {
      output <- rbind(output, candidate)
    }
  }
  output <- output[-1, ]
  output
}

symmetry_5_VII <- function(leaves) {
  output <- character(5)
  all_permutations <- make_permutations(leaves)
  for (candidate in all_permutations) {
    output <- rbind(output, candidate)
  }
  output <- output[-1, ]
  output
}

#' Six leaves graphs.
#' 
#' Kind of obsolete since the introduction of \code{\link{all_graphs}}.
#' A comprehensive listing of all the \eqn{21} admixture graphs with six leaves and
#' at most one admixture event. Our convention is that the position of the root does
#' not matter (as long as it's not after an admixture event) and that graphs that have
#' \emph{eyes}, two inner nodes with the property that all the paths between any two
#' leaves visits both or neither of them, are excluded. The reason  is that the \eqn{f}
#' statistics can't detect the exact position of the root or distinguish between an
#' eye and a simple branch. The position of the root can be moved later with the function
#' \code{\link{make_an_outgroup}}.
#' 
#' @format A list of functions on six leaves and a parameter \code{permutations} which
#'         is \code{FALSE} by default.
#'         The outputs of these functions are either single \code{\link{agraph}} objects
#'         with the input vector as leaves, or if \code{permutations} is \code{TRUE},
#'         lists of all the possible \code{\link{agraph}} objects with that leaf set up
#'         to symmetry.
#'
#' @family graphs
#'
#' @seealso \code{\link{all_graphs}}
#' @seealso \code{\link{make_permutations}}
#' @seealso \code{\link{fit_permutations_and_graphs}}
#' @seealso \code{\link{fit_graph_list}}
#' @seealso \code{\link{add_a_leaf}}
#' @seealso \code{\link{add_an_admixture}}
#' @seealso \code{\link{add_an_admixture2}}
#' @seealso \code{\link{make_an_outgroup}}
#'
#' @examples
#' \donttest{
#' # While the usage of this function is pretty self-explanatory, let's plot all the graphs
#' # just for browsing.
#' for (i in seq(1, length(six_leaves_graphs))) {
#'   graph <- six_leaves_graphs[[i]](c("A", "B", "C", "D", "E", "F"))
#'   # This is how you include quotation marks in strings by the way:
#'   title <- paste("six_leaves_graphs[[", i,
#'                  "]](c(\"A\", \"B\", \"C\", \"D\", \"E\", \"F\"))", sep = "")
#'   plot(graph, color = "yellow4", title = title)
#' }
#' }
#'
#' @export
six_leaves_graphs <- list(
  # tree
  tree_1 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 15
      leaf_permutations <- symmetry_6_XVII(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "w")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "y"),
        edge("w", "y"),
        edge(leaves[1], "R"),
        edge(leaves[2], "x"),
        edge(leaves[3], "z"),
        edge(leaves[4], "z"),
        edge(leaves[5], "w"),
        edge(leaves[6], "w")
      ))
      admixtures <- NULL
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },

  tree_2 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 90
      leaf_permutations <- symmetry_6_XVIII(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "w")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "y"),
        edge("w", "z"),
        edge(leaves[1], "R"),
        edge(leaves[2], "x"),
        edge(leaves[3], "y"),
        edge(leaves[4], "z"),
        edge(leaves[5], "w"),
        edge(leaves[6], "w")
      ))
      admixtures <- NULL
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  # one admixture
  one_admixture_1 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 360
      leaf_permutations <- symmetry_6_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "w", "u", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "x"),
        edge("w", "y"),
        edge("u", "z"),
        edge(leaves[1], "R"),
        edge(leaves[2], "y"), 
        edge(leaves[3], "w"),
        edge(leaves[4], "M"),
        edge(leaves[5], "u"),
        edge(leaves[6], "z"),
        admixture_edge("M", "w", "u")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "w", "u", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  one_admixture_2 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 180
      leaf_permutations <- symmetry_6_II(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "w", "u", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "x"),
        edge("w", "y"),
        edge("u", "M"),
        edge(leaves[1], "x"),
        edge(leaves[2], "z"), 
        edge(leaves[3], "u"),
        edge(leaves[4], "u"),
        edge(leaves[5], "w"),
        edge(leaves[6], "y"),
        admixture_edge("M", "z", "w")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "z", "w", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  one_admixture_3 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 360
      leaf_permutations <- symmetry_6_III(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "w", "u", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "x"),
        edge("w", "y"),
        edge("u", "w"),
        edge(leaves[1], "x"),
        edge(leaves[2], "z"), 
        edge(leaves[3], "M"),
        edge(leaves[4], "u"),
        edge(leaves[5], "u"),
        edge(leaves[6], "y"),
        admixture_edge("M", "z", "w")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "z", "w", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },

  one_admixture_4 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 360
      leaf_permutations <- symmetry_6_IV(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "w", "u", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "x"),
        edge("w", "y"),
        edge("u", "y"),
        edge(leaves[1], "x"),
        edge(leaves[2], "z"), 
        edge(leaves[3], "M"),
        edge(leaves[4], "w"),
        edge(leaves[5], "u"),
        edge(leaves[6], "u"),
        admixture_edge("M", "z", "w")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "z", "w", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  one_admixture_5 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 180
      leaf_permutations <- symmetry_6_V(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "w", "u", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "x"),
        edge("w", "M"),
        edge("u", "w"),
        edge(leaves[1], "R"),
        edge(leaves[2], "y"), 
        edge(leaves[3], "w"),
        edge(leaves[4], "u"),
        edge(leaves[5], "u"),
        edge(leaves[6], "z"),
        admixture_edge("M", "y", "z")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "y", "z", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  one_admixture_6 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 360
      leaf_permutations <- symmetry_6_VI(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "w", "u", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "x"),
        edge("w", "z"),
        edge("u", "w"),
        edge(leaves[1], "R"),
        edge(leaves[2], "y"), 
        edge(leaves[3], "M"),
        edge(leaves[4], "w"),
        edge(leaves[5], "u"),
        edge(leaves[6], "u"),
        admixture_edge("M", "y", "z")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "y", "z", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  one_admixture_7 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 180
      leaf_permutations <- symmetry_6_VII(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "w", "u", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "x"),
        edge("w", "R"),
        edge("u", "w"),
        edge(leaves[1], "u"),
        edge(leaves[2], "u"), 
        edge(leaves[3], "w"),
        edge(leaves[4], "y"),
        edge(leaves[5], "M"),
        edge(leaves[6], "z"),
        admixture_edge("M", "y", "z")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "y", "z", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  one_admixture_8 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 180
      leaf_permutations <- symmetry_6_VIII(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "w", "u", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "x"),
        edge("w", "M"),
        edge("u", "z"),
        edge(leaves[1], "R"),
        edge(leaves[2], "y"), 
        edge(leaves[3], "w"),
        edge(leaves[4], "w"),
        edge(leaves[5], "u"),
        edge(leaves[6], "u"),
        admixture_edge("M", "y", "z")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "y", "z", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  one_admixture_9 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 180
      leaf_permutations <- symmetry_6_IX(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "w", "u", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "x"),
        edge("w", "R"),
        edge("u", "z"),
        edge(leaves[1], "w"),
        edge(leaves[2], "w"), 
        edge(leaves[3], "y"),
        edge(leaves[4], "M"),
        edge(leaves[5], "u"),
        edge(leaves[6], "u"),
        admixture_edge("M", "y", "z")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "y", "z", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  one_admixture_10 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 90
      leaf_permutations <- symmetry_6_X(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "w", "u", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "x"),
        edge("w", "R"),
        edge("u", "M"),
        edge(leaves[1], "w"),
        edge(leaves[2], "w"), 
        edge(leaves[3], "y"),
        edge(leaves[4], "u"),
        edge(leaves[5], "u"),
        edge(leaves[6], "z"),
        admixture_edge("M", "y", "z")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "y", "z", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  one_admixture_11 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 90
      leaf_permutations <- symmetry_6_XI(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "w", "u", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "x"),
        edge("w", "y"),
        edge("u", "z"),
        edge(leaves[1], "R"),
        edge(leaves[2], "w"), 
        edge(leaves[3], "w"),
        edge(leaves[4], "M"),
        edge(leaves[5], "u"),
        edge(leaves[6], "u"),
        admixture_edge("M", "y", "z")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "y", "z", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  one_admixture_12 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 180
      leaf_permutations <- symmetry_6_XII(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "w", "u", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "M"),
        edge("w", "z"),
        edge("u", "w"),
        edge(leaves[1], "x"),
        edge(leaves[2], "z"), 
        edge(leaves[3], "w"),
        edge(leaves[4], "u"),
        edge(leaves[5], "u"),
        edge(leaves[6], "y"),
        admixture_edge("M", "x", "y")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "y", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  one_admixture_13 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 360
      leaf_permutations <- symmetry_6_VI(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "w", "u", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "y"),
        edge("w", "z"),
        edge("u", "w"),
        edge(leaves[1], "x"),
        edge(leaves[2], "M"), 
        edge(leaves[3], "z"),
        edge(leaves[4], "w"),
        edge(leaves[5], "u"),
        edge(leaves[6], "u"),
        admixture_edge("M", "x", "y")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "y", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  one_admixture_14 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 45
      leaf_permutations <- symmetry_6_XIII(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "w", "u", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "M"),
        edge("w", "z"),
        edge("u", "z"),
        edge(leaves[1], "x"),
        edge(leaves[2], "w"), 
        edge(leaves[3], "w"),
        edge(leaves[4], "u"),
        edge(leaves[5], "u"),
        edge(leaves[6], "y"),
        admixture_edge("M", "x", "y")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "y", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  one_admixture_15 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 90
      leaf_permutations <- symmetry_6_XIV(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "w", "u", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "y"),
        edge("w", "z"),
        edge("u", "z"),
        edge(leaves[1], "x"),
        edge(leaves[2], "M"), 
        edge(leaves[3], "w"),
        edge(leaves[4], "w"),
        edge(leaves[5], "u"),
        edge(leaves[6], "u"),
        admixture_edge("M", "x", "y")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "y", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  one_admixture_16 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 180
      leaf_permutations <- symmetry_6_VIII(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "w", "u", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "M"),
        edge("w", "z"),
        edge("u", "y"),
        edge(leaves[1], "x"),
        edge(leaves[2], "z"), 
        edge(leaves[3], "w"),
        edge(leaves[4], "w"),
        edge(leaves[5], "u"),
        edge(leaves[6], "u"),
        admixture_edge("M", "x", "y")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "y", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  one_admixture_17 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 180
      leaf_permutations <- symmetry_6_XV(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "w", "u", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "M"),
        edge("w", "y"),
        edge("u", "w"),
        edge(leaves[1], "x"),
        edge(leaves[2], "z"), 
        edge(leaves[3], "z"),
        edge(leaves[4], "w"),
        edge(leaves[5], "u"),
        edge(leaves[6], "u"),
        admixture_edge("M", "x", "y")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "y", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  one_admixture_18 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 180
      leaf_permutations <- symmetry_6_IX(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "w", "u", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "x"),
        edge("w", "y"),
        edge("u", "w"),
        edge(leaves[1], "z"),
        edge(leaves[2], "z"), 
        edge(leaves[3], "M"),
        edge(leaves[4], "w"),
        edge(leaves[5], "u"),
        edge(leaves[6], "u"),
        admixture_edge("M", "x", "y")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "y", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  one_admixture_19 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 45
      leaf_permutations <- symmetry_6_XVI(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "w", "u", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "x"),
        edge("w", "M"),
        edge("u", "y"),
        edge(leaves[1], "z"),
        edge(leaves[2], "z"), 
        edge(leaves[3], "w"),
        edge(leaves[4], "w"),
        edge(leaves[5], "u"),
        edge(leaves[6], "u"),
        admixture_edge("M", "x", "y")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "y", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  }
)

symmetry_6_I <- function(leaves) {
  output <- character(6)
  all_permutations <- make_permutations(leaves)
  for (candidate in all_permutations) {
    sorted <- sort(candidate)
    num <- numeric(0)
    for (j in seq(1, 6)) {
      num[j] <- match(candidate[j], sorted)
    }
    if (num[3] < num[5]) {
      output <- rbind(output, candidate)
    }
  }
  output <- output[-1, ]
  output
}

symmetry_6_II <- function(leaves) {
  output <- character(6)
  all_permutations <- make_permutations(leaves)
  for (candidate in all_permutations) {
    sorted <- sort(candidate)
    num <- numeric(0)
    for (j in seq(1, 6)) {
      num[j] <- match(candidate[j], sorted)
    }
    if (num[3] < num[4] &&
        num[1] < num[6]) {
      output <- rbind(output, candidate)
    }
  }
  output <- output[-1, ]
  output
}

symmetry_6_III <- function(leaves) {
  output <- character(6)
  all_permutations <- make_permutations(leaves)
  for (candidate in all_permutations) {
    sorted <- sort(candidate)
    num <- numeric(0)
    for (j in seq(1, 6)) {
      num[j] <- match(candidate[j], sorted)
    }
    if (num[4] < num[5]) {
      output <- rbind(output, candidate)
    }
  }
  output <- output[-1, ]
  output
}

symmetry_6_IV <- function(leaves) {
  output <- character(6)
  all_permutations <- make_permutations(leaves)
  for (candidate in all_permutations) {
    sorted <- sort(candidate)
    num <- numeric(0)
    for (j in seq(1, 6)) {
      num[j] <- match(candidate[j], sorted)
    }
    if (num[5] < num[6]) {
      output <- rbind(output, candidate)
    }
  }
  output <- output[-1, ]
  output
}

symmetry_6_V <- function(leaves) {
  output <- character(6)
  all_permutations <- make_permutations(leaves)
  for (candidate in all_permutations) {
    sorted <- sort(candidate)
    num <- numeric(0)
    for (j in seq(1, 6)) {
      num[j] <- match(candidate[j], sorted)
    }
    if (num[4] < num[5] &&
        num[2] < num[6]) {
      output <- rbind(output, candidate)
    }
  }
  output <- output[-1, ]
  output
}

symmetry_6_VI <- function(leaves) {
  output <- character(6)
  all_permutations <- make_permutations(leaves)
  for (candidate in all_permutations) {
    sorted <- sort(candidate)
    num <- numeric(0)
    for (j in seq(1, 6)) {
      num[j] <- match(candidate[j], sorted)
    }
    if (num[5] < num[6]) {
      output <- rbind(output, candidate)
    }
  }
  output <- output[-1, ]
  output
}

symmetry_6_VII <- function(leaves) {
  output <- character(6)
  all_permutations <- make_permutations(leaves)
  for (candidate in all_permutations) {
    sorted <- sort(candidate)
    num <- numeric(0)
    for (j in seq(1, 6)) {
      num[j] <- match(candidate[j], sorted)
    }
    if (num[1] < num[2] &&
        num[4] < num[6]) {
      output <- rbind(output, candidate)
    }
  }
  output <- output[-1, ]
  output
}

symmetry_6_VIII <- function(leaves) {
  output <- character(6)
  all_permutations <- make_permutations(leaves)
  for (candidate in all_permutations) {
    sorted <- sort(candidate)
    num <- numeric(0)
    for (j in seq(1, 6)) {
      num[j] <- match(candidate[j], sorted)
    }
    if (num[3] < num[4] &&
        num[5] < num[6]) {
      output <- rbind(output, candidate)
    }
  }
  output <- output[-1, ]
  output
}

symmetry_6_IX <- function(leaves) {
  output <- character(6)
  all_permutations <- make_permutations(leaves)
  for (candidate in all_permutations) {
    sorted <- sort(candidate)
    num <- numeric(0)
    for (j in seq(1, 6)) {
      num[j] <- match(candidate[j], sorted)
    }
    if (num[1] < num[2] &&
        num[5] < num[6]) {
      output <- rbind(output, candidate)
    }
  }
  output <- output[-1, ]
  output
}

symmetry_6_X <- function(leaves) {
  output <- character(6)
  all_permutations <- make_permutations(leaves)
  for (candidate in all_permutations) {
    sorted <- sort(candidate)
    num <- numeric(0)
    for (j in seq(1, 6)) {
      num[j] <- match(candidate[j], sorted)
    }
    if (num[1] < num[2] &&
        num[3] < num[6] &&
        num[4] < num[5]) {
      output <- rbind(output, candidate)
    }
  }
  output <- output[-1, ]
  output
}

symmetry_6_XI <- function(leaves) {
  output <- character(6)
  all_permutations <- make_permutations(leaves)
  for (candidate in all_permutations) {
    sorted <- sort(candidate)
    num <- numeric(0)
    for (j in seq(1, 6)) {
      num[j] <- match(candidate[j], sorted)
    }
    if (num[2] < num[3] &&
        num[5] < num[6] &&
        num[2] + 10*num[3] < num[5] + 10*num[6]) {
      output <- rbind(output, candidate)
    }
  }
  output <- output[-1, ]
  output
}

symmetry_6_XII <- function(leaves) {
  output <- character(6)
  all_permutations <- make_permutations(leaves)
  for (candidate in all_permutations) {
    sorted <- sort(candidate)
    num <- numeric(0)
    for (j in seq(1, 6)) {
      num[j] <- match(candidate[j], sorted)
    }
    if (num[1] < num[6] &&
        num[4] < num[5]) {
      output <- rbind(output, candidate)
    }
  }
  output <- output[-1, ]
  output
}

symmetry_6_XIII <- function(leaves) {
  output <- character(6)
  all_permutations <- make_permutations(leaves)
  for (candidate in all_permutations) {
    sorted <- sort(candidate)
    num <- numeric(0)
    for (j in seq(1, 6)) {
      num[j] <- match(candidate[j], sorted)
    }
    if (num[1] < num[6] &&
        num[2] < num[3] &&
        num[4] < num[5] &&
        num[2] + 10*num[3] < num[4] + 10*num[5]) {
      output <- rbind(output, candidate)
    }
  }
  output <- output[-1, ]
  output
}

symmetry_6_XIV <- function(leaves) {
  output <- character(6)
  all_permutations <- make_permutations(leaves)
  for (candidate in all_permutations) {
    sorted <- sort(candidate)
    num <- numeric(0)
    for (j in seq(1, 6)) {
      num[j] <- match(candidate[j], sorted)
    }
    if (num[3] < num[4] &&
        num[5] < num[6] &&
        num[3] + 10*num[4] < num[5] + 10*num[6]) {
      output <- rbind(output, candidate)
    }
  }
  output <- output[-1, ]
  output
}

symmetry_6_XV <- function(leaves) {
  output <- character(6)
  all_permutations <- make_permutations(leaves)
  for (candidate in all_permutations) {
    sorted <- sort(candidate)
    num <- numeric(0)
    for (j in seq(1, 6)) {
      num[j] <- match(candidate[j], sorted)
    }
    if (num[2] < num[3] &&
        num[5] < num[6]) {
      output <- rbind(output, candidate)
    }
  }
  output <- output[-1, ]
  output
}

symmetry_6_XVI <- function(leaves) {
  output <- character(6)
  all_permutations <- make_permutations(leaves)
  for (candidate in all_permutations) {
    sorted <- sort(candidate)
    num <- numeric(0)
    for (j in seq(1, 6)) {
      num[j] <- match(candidate[j], sorted)
    }
    if (num[3] < num[4] &&
        num[1] < num[2] &&
        num[5] < num[6] &&
        num[1] + 10*num[2] < num[5] + 10*num[6]) {
      output <- rbind(output, candidate)
    }
  }
  output <- output[-1, ]
  output
}

symmetry_6_XVII <- function(leaves) {
  output <- character(6)
  all_permutations <- make_permutations(leaves)
  for (candidate in all_permutations) {
    sorted <- sort(candidate)
    num <- numeric(0)
    for (j in seq(1, 6)) {
      num[j] <- match(candidate[j], sorted)
    }
    if (num[1] < num[2] &&
        num[3] < num[4] &&
        num[5] < num[6] &&
        num[1] + 10*num[2] < num[3] + 10*num[4] &&
        num[3] + 10*num[4] < num[5] + 10*num[6]) {
      output <- rbind(output, candidate)
    }
  }
  output <- output[-1, ]
  output
}

symmetry_6_XVIII <- function(leaves) {
  output <- character(6)
  all_permutations <- make_permutations(leaves)
  for (candidate in all_permutations) {
    sorted <- sort(candidate)
    num <- numeric(0)
    for (j in seq(1, 6)) {
      num[j] <- match(candidate[j], sorted)
    }
    if (num[1] < num[2] &&
        num[3] < num[4] &&
        num[5] < num[6]) {
      output <- rbind(output, candidate)
    }
  }
  output <- output[-1, ]
  output
}

#' Seven leaves graphs.
#' 
#' Kind of obsolete since the introduction of \code{\link{all_graphs}}.
#' A comprehensive listing of all the \eqn{48} admixture graphs with seven leaves and
#' at most one admixture event. Our convention is that the position of the root does
#' not matter (as long as it's not after an admixture event) and that graphs that have
#' \emph{eyes}, two inner nodes with the property that all the paths between any two
#' leaves visits both or neither of them, are excluded. The reason  is that the \eqn{f}
#' statistics can't detect the exact position of the root or distinguish between an
#' eye and a simple branch. The position of the root can be moved later with the function
#' \code{\link{make_an_outgroup}}.
#' 
#' @format A list of functions on seven leaves and a parameter \code{permutations} which
#'         is \code{FALSE} by default.
#'         The outputs of these functions are either single \code{\link{agraph}} objects
#'         with the input vector as leaves, or if \code{permutations} is \code{TRUE},
#'         lists of all the possible \code{\link{agraph}} objects with that leaf set up
#'         to symmetry.
#'
#' @family graphs
#'   
#' @seealso \code{\link{all_graphs}}
#' @seealso \code{\link{make_permutations}}
#' @seealso \code{\link{fit_permutations_and_graphs}}
#' @seealso \code{\link{fit_graph_list}}
#' @seealso \code{\link{add_a_leaf}}
#' @seealso \code{\link{add_an_admixture}}
#' @seealso \code{\link{add_an_admixture2}}
#' @seealso \code{\link{make_an_outgroup}}
#'
#' @examples
#' \donttest{
#' # While the usage of this function is pretty self-explanatory, let's plot all the graphs
#' # just for browsing.
#' for (i in seq(1, length(seven_leaves_graphs))) {
#'   graph <- seven_leaves_graphs[[i]](c("A", "B", "C", "D", "E", "F", "G"))
#'   # This is how you include quotation marks in strings by the way:
#'   title <- paste("seven_leaves_graphs[[", i,
#'                  "]](c(\"A\", \"B\", \"C\", \"D\", \"E\", \"F\", \"G\"))", sep = "")
#'   plot(graph, color = "seagreen", title = title)
#' }
#' }
#'
#' @export
seven_leaves_graphs <- list(
  # tree
  tree_1 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 630
      leaf_permutations <- symmetry_7_III(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "y"),
        edge("u", "z"),
        edge("v", "u"),
        edge(leaves[1], "R"),
        edge(leaves[2], "x"),
        edge(leaves[3], "y"),
        edge(leaves[4], "u"),
        edge(leaves[5], "v"),
        edge(leaves[6], "v"),
        edge(leaves[7], "z")
      ))
      admixtures <- NULL
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },

  tree_2 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 315
      leaf_permutations <- symmetry_7_IV(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "x"),
        edge("u", "R"),
        edge("v", "u"),
        edge(leaves[1], "y"),
        edge(leaves[2], "y"),
        edge(leaves[3], "z"),
        edge(leaves[4], "z"),
        edge(leaves[5], "v"),
        edge(leaves[6], "v"),
        edge(leaves[7], "u")
      ))
      admixtures <- NULL
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  # one admixture
  one_admixture_1 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 2520
      leaf_permutations <- symmetry_7_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "w", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "y"),
        edge("u", "R"),
        edge("v", "u"),
        edge("w", "v"),
        edge(leaves[1], "u"),
        edge(leaves[2], "x"),
        edge(leaves[3], "y"),
        edge(leaves[4], "z"),
        edge(leaves[5], "M"),
        edge(leaves[6], "w"),
        edge(leaves[7], "v"),
        admixture_edge("M", "z", "w")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "z", "w", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  one_admixture_2 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 1260
      leaf_permutations <- symmetry_7_II(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "w", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "x"),
        edge("u", "y"),
        edge("v", "z"),
        edge("w", "M"),
        edge(leaves[1], "w"),
        edge(leaves[2], "w"),
        edge(leaves[3], "u"),
        edge(leaves[4], "v"),
        edge(leaves[5], "z"),
        edge(leaves[6], "R"),
        edge(leaves[7], "y"),
        admixture_edge("M", "u", "v")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "u", "v", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  one_admixture_3 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 2520
      leaf_permutations <- symmetry_7_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "w", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "x"),
        edge("u", "y"),
        edge("v", "z"),
        edge("w", "v"),
        edge(leaves[1], "w"),
        edge(leaves[2], "w"),
        edge(leaves[3], "z"),
        edge(leaves[4], "R"),
        edge(leaves[5], "y"),
        edge(leaves[6], "u"),
        edge(leaves[7], "M"),
        admixture_edge("M", "u", "v")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "u", "v", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  one_admixture_4 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 2520
      leaf_permutations <- symmetry_7_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "w", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "x"),
        edge("u", "y"),
        edge("v", "z"),
        edge("w", "z"),
        edge(leaves[1], "w"),
        edge(leaves[2], "w"),
        edge(leaves[3], "R"),
        edge(leaves[4], "y"),
        edge(leaves[5], "u"),
        edge(leaves[6], "M"),
        edge(leaves[7], "v"),
        admixture_edge("M", "u", "v")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "u", "v", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  one_admixture_5 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 1260
      leaf_permutations <- symmetry_7_II(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "w", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "y"),
        edge("u", "y"),
        edge("v", "z"),
        edge("w", "u"),
        edge(leaves[1], "R"),
        edge(leaves[2], "x"),
        edge(leaves[3], "z"),
        edge(leaves[4], "u"),
        edge(leaves[5], "w"),
        edge(leaves[6], "M"),
        edge(leaves[7], "v"),
        admixture_edge("M", "v", "w")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "v", "w", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  one_admixture_6 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 1260
      leaf_permutations <- symmetry_7_II(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "w", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "x"),
        edge("u", "y"),
        edge("v", "M"),
        edge("w", "u"),
        edge(leaves[1], "w"),
        edge(leaves[2], "w"),
        edge(leaves[3], "v"),
        edge(leaves[4], "v"),
        edge(leaves[5], "z"),
        edge(leaves[6], "x"),
        edge(leaves[7], "y"),
        admixture_edge("M", "z", "u")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "z", "u", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  one_admixture_7 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 1260
      leaf_permutations <- symmetry_7_II(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "w", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "x"),
        edge("u", "y"),
        edge("v", "u"),
        edge("w", "y"),
        edge(leaves[1], "w"),
        edge(leaves[2], "w"),
        edge(leaves[3], "v"),
        edge(leaves[4], "v"),
        edge(leaves[5], "M"),
        edge(leaves[6], "z"),
        edge(leaves[7], "x"),
        admixture_edge("M", "z", "u")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "z", "u", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  one_admixture_8 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 630
      leaf_permutations <- symmetry_7_III(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "w", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "x"),
        edge("u", "y"),
        edge("v", "x"),
        edge("w", "y"),
        edge(leaves[1], "w"),
        edge(leaves[2], "w"),
        edge(leaves[3], "v"),
        edge(leaves[4], "v"),
        edge(leaves[5], "z"),
        edge(leaves[6], "M"),
        edge(leaves[7], "u"),
        admixture_edge("M", "z", "u")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "z", "u", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  one_admixture_9 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 630
      leaf_permutations <- symmetry_7_III(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "w", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "x"),
        edge("u", "y"),
        edge("v", "z"),
        edge("w", "u"),
        edge(leaves[1], "w"),
        edge(leaves[2], "w"),
        edge(leaves[3], "y"),
        edge(leaves[4], "x"),
        edge(leaves[5], "v"),
        edge(leaves[6], "v"),
        edge(leaves[7], "M"),
        admixture_edge("M", "z", "u")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "z", "u", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  one_admixture_10 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 1260
      leaf_permutations <- symmetry_7_II(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "w", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "x"),
        edge("u", "y"),
        edge("v", "M"),
        edge("w", "y"),
        edge(leaves[1], "v"),
        edge(leaves[2], "v"),
        edge(leaves[3], "w"),
        edge(leaves[4], "w"),
        edge(leaves[5], "x"),
        edge(leaves[6], "z"),
        edge(leaves[7], "u"),
        admixture_edge("M", "z", "u")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "z", "u", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  one_admixture_11 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 1260
      leaf_permutations <- symmetry_7_II(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "w", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "x"),
        edge("u", "y"),
        edge("v", "z"),
        edge("w", "y"),
        edge(leaves[1], "w"),
        edge(leaves[2], "w"),
        edge(leaves[3], "v"),
        edge(leaves[4], "v"),
        edge(leaves[5], "M"),
        edge(leaves[6], "u"),
        edge(leaves[7], "x"),
        admixture_edge("M", "z", "u")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "z", "u", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  one_admixture_12 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 2520
      leaf_permutations <- symmetry_7_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "w", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "y"),
        edge("u", "z"),
        edge("v", "u"),
        edge("w", "z"),
        edge(leaves[1], "R"),
        edge(leaves[2], "x"),
        edge(leaves[3], "y"),
        edge(leaves[4], "u"),
        edge(leaves[5], "v"),
        edge(leaves[6], "M"),
        edge(leaves[7], "w"),
        admixture_edge("M", "v", "w")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "v", "w", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  one_admixture_13 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 2520
      leaf_permutations <- symmetry_7_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "w", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "x"),
        edge("u", "z"),
        edge("v", "y"),
        edge("w", "v"),
        edge(leaves[1], "w"),
        edge(leaves[2], "w"),
        edge(leaves[3], "v"),
        edge(leaves[4], "x"),
        edge(leaves[5], "z"),
        edge(leaves[6], "u"),
        edge(leaves[7], "M"),
        admixture_edge("M", "u", "y")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "u", "y", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  one_admixture_14 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 1260
      leaf_permutations <- symmetry_7_II(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "w", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "x"),
        edge("u", "y"),
        edge("v", "M"),
        edge("w", "v"),
        edge(leaves[1], "w"),
        edge(leaves[2], "w"),
        edge(leaves[3], "x"),
        edge(leaves[4], "y"),
        edge(leaves[5], "v"),
        edge(leaves[6], "z"),
        edge(leaves[7], "u"),
        admixture_edge("M", "z", "u")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "z", "u", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  one_admixture_15 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 315
      leaf_permutations <- symmetry_7_IV(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "w", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "y"),
        edge("u", "y"),
        edge("v", "z"),
        edge("w", "u"),
        edge(leaves[1], "v"),
        edge(leaves[2], "v"),
        edge(leaves[3], "w"),
        edge(leaves[4], "w"),
        edge(leaves[5], "x"),
        edge(leaves[6], "x"),
        edge(leaves[7], "M"),
        admixture_edge("M", "z", "u")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "z", "u", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  one_admixture_16 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 630
      leaf_permutations <- symmetry_7_III(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "w", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "y"),
        edge("u", "y"),
        edge("v", "M"),
        edge("w", "u"),
        edge(leaves[1], "x"),
        edge(leaves[2], "R"),
        edge(leaves[3], "w"),
        edge(leaves[4], "w"),
        edge(leaves[5], "v"),
        edge(leaves[6], "v"),
        edge(leaves[7], "z"),
        admixture_edge("M", "z", "u")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "z", "u", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  one_admixture_17 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 315
      leaf_permutations <- symmetry_7_IV(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "w", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "x"),
        edge("u", "y"),
        edge("v", "M"),
        edge("w", "z"),
        edge(leaves[1], "u"),
        edge(leaves[2], "u"),
        edge(leaves[3], "w"),
        edge(leaves[4], "w"),
        edge(leaves[5], "v"),
        edge(leaves[6], "v"),
        edge(leaves[7], "R"),
        admixture_edge("M", "y", "z")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "y", "z", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  one_admixture_18 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 1260
      leaf_permutations <- symmetry_7_II(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "w", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "x"),
        edge("u", "y"),
        edge("v", "M"),
        edge("w", "v"),
        edge(leaves[1], "u"),
        edge(leaves[2], "u"),
        edge(leaves[3], "w"),
        edge(leaves[4], "w"),
        edge(leaves[5], "v"),
        edge(leaves[6], "z"),
        edge(leaves[7], "R"),
        admixture_edge("M", "y", "z")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "y", "z", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  one_admixture_19 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 1260
      leaf_permutations <- symmetry_7_II(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "w", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "x"),
        edge("u", "M"),
        edge("v", "z"),
        edge("w", "v"),
        edge(leaves[1], "u"),
        edge(leaves[2], "u"),
        edge(leaves[3], "w"),
        edge(leaves[4], "w"),
        edge(leaves[5], "v"),
        edge(leaves[6], "R"),
        edge(leaves[7], "y"),
        admixture_edge("M", "y", "z")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "y", "z", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  one_admixture_20 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 1260
      leaf_permutations <- symmetry_7_II(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "w", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "y"),
        edge("u", "z"),
        edge("v", "z"),
        edge("w", "v"),
        edge(leaves[1], "w"),
        edge(leaves[2], "w"),
        edge(leaves[3], "R"),
        edge(leaves[4], "x"),
        edge(leaves[5], "y"),
        edge(leaves[6], "u"),
        edge(leaves[7], "M"),
        admixture_edge("M", "u", "v")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "u", "v", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  one_admixture_21 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 1260
      leaf_permutations <- symmetry_7_II(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "w", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "y"),
        edge("u", "y"),
        edge("v", "z"),
        edge("w", "v"),
        edge(leaves[1], "R"),
        edge(leaves[2], "x"),
        edge(leaves[3], "w"),
        edge(leaves[4], "w"),
        edge(leaves[5], "v"),
        edge(leaves[6], "M"),
        edge(leaves[7], "u"),
        admixture_edge("M", "z", "u")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "z", "u", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  one_admixture_22 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 1260
      leaf_permutations <- symmetry_7_II(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "w", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "x"),
        edge("u", "y"),
        edge("v", "z"),
        edge("w", "v"),
        edge(leaves[1], "u"),
        edge(leaves[2], "u"),
        edge(leaves[3], "w"),
        edge(leaves[4], "w"),
        edge(leaves[5], "v"),
        edge(leaves[6], "R"),
        edge(leaves[7], "M"),
        admixture_edge("M", "y", "z")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "y", "z", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  one_admixture_23 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 630
      leaf_permutations <- symmetry_7_III(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "w", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "y"),
        edge("u", "z"),
        edge("v", "z"),
        edge("w", "M"),
        edge(leaves[1], "w"),
        edge(leaves[2], "w"),
        edge(leaves[3], "u"),
        edge(leaves[4], "v"),
        edge(leaves[5], "R"),
        edge(leaves[6], "x"),
        edge(leaves[7], "y"),
        admixture_edge("M", "u", "v")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "u", "v", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  one_admixture_24 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 630
      leaf_permutations <- symmetry_7_III(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "w", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "y"),
        edge("u", "y"),
        edge("v", "M"),
        edge("w", "v"),
        edge(leaves[1], "R"),
        edge(leaves[2], "x"),
        edge(leaves[3], "z"),
        edge(leaves[4], "u"),
        edge(leaves[5], "w"),
        edge(leaves[6], "w"),
        edge(leaves[7], "v"),
        admixture_edge("M", "z", "u")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "z", "u", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  one_admixture_25 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 630
      leaf_permutations <- symmetry_7_V(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "w", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "x"),
        edge("u", "z"),
        edge("v", "u"),
        edge("w", "u"),
        edge(leaves[1], "w"),
        edge(leaves[2], "w"),
        edge(leaves[3], "v"),
        edge(leaves[4], "v"),
        edge(leaves[5], "M"),
        edge(leaves[6], "y"),
        edge(leaves[7], "R"),
        admixture_edge("M", "y", "z")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "y", "z", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  one_admixture_26 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 315
      leaf_permutations <- symmetry_7_IV(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "w", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "x"),
        edge("u", "M"),
        edge("v", "u"),
        edge("w", "u"),
        edge(leaves[1], "v"),
        edge(leaves[2], "v"),
        edge(leaves[3], "w"),
        edge(leaves[4], "w"),
        edge(leaves[5], "y"),
        edge(leaves[6], "z"),
        edge(leaves[7], "R"),
        admixture_edge("M", "y", "z")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "y", "z", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  one_admixture_27 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 315
      leaf_permutations <- symmetry_7_IV(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "w", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "x"),
        edge("u", "R"),
        edge("v", "u"),
        edge("w", "u"),
        edge(leaves[1], "y"),
        edge(leaves[2], "y"),
        edge(leaves[3], "z"),
        edge(leaves[4], "z"),
        edge(leaves[5], "v"),
        edge(leaves[6], "w"),
        edge(leaves[7], "M"),
        admixture_edge("M", "v", "w")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "v", "w", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  one_admixture_28 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 1260
      leaf_permutations <- symmetry_7_II(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "w", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "x"),
        edge("u", "M"),
        edge("v", "u"),
        edge("w", "v"),
        edge(leaves[1], "w"),
        edge(leaves[2], "w"),
        edge(leaves[3], "y"),
        edge(leaves[4], "z"),
        edge(leaves[5], "v"),
        edge(leaves[6], "u"),
        edge(leaves[7], "R"),
        admixture_edge("M", "y", "z")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "y", "z", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  one_admixture_29 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 2520
      leaf_permutations <- symmetry_7_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "w", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "x"),
        edge("u", "z"),
        edge("v", "u"),
        edge("w", "v"),
        edge(leaves[1], "w"),
        edge(leaves[2], "w"),
        edge(leaves[3], "v"),
        edge(leaves[4], "u"),
        edge(leaves[5], "M"),
        edge(leaves[6], "y"),
        edge(leaves[7], "R"),
        admixture_edge("M", "y", "z")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "y", "z", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  one_admixture_30 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 1260
      leaf_permutations <- symmetry_7_II(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "w", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "x"),
        edge("z", "y"),
        edge("u", "z"),
        edge("v", "u"),
        edge("w", "u"),
        edge(leaves[1], "R"),
        edge(leaves[2], "x"),
        edge(leaves[3], "v"),
        edge(leaves[4], "w"),
        edge(leaves[5], "y"),
        edge(leaves[6], "z"),
        edge(leaves[7], "M"),
        admixture_edge("M", "v", "w")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "v", "w", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  one_admixture_31 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 315
      leaf_permutations <- symmetry_7_IV(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "w", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "x"),
        edge("u", "y"),
        edge("v", "M"),
        edge("w", "v"),
        edge(leaves[1], "u"),
        edge(leaves[2], "u"),
        edge(leaves[3], "z"),
        edge(leaves[4], "z"),
        edge(leaves[5], "w"),
        edge(leaves[6], "w"),
        edge(leaves[7], "v"),
        admixture_edge("M", "x", "y")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "y", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  one_admixture_32 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 630
      leaf_permutations <- symmetry_7_III(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "w", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "x"),
        edge("u", "z"),
        edge("v", "M"),
        edge("w", "y"),
        edge(leaves[1], "v"),
        edge(leaves[2], "v"),
        edge(leaves[3], "w"),
        edge(leaves[4], "w"),
        edge(leaves[5], "u"),
        edge(leaves[6], "u"),
        edge(leaves[7], "z"),
        admixture_edge("M", "x", "y")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "y", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  one_admixture_33 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 630
      leaf_permutations <- symmetry_7_III(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "w", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "x"),
        edge("u", "y"),
        edge("v", "z"),
        edge("w", "u"),
        edge(leaves[1], "v"),
        edge(leaves[2], "v"),
        edge(leaves[3], "w"),
        edge(leaves[4], "w"),
        edge(leaves[5], "z"),
        edge(leaves[6], "u"),
        edge(leaves[7], "M"),
        admixture_edge("M", "x", "y")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "y", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  one_admixture_34 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 1260
      leaf_permutations <- symmetry_7_II(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "w", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "x"),
        edge("u", "z"),
        edge("v", "M"),
        edge("w", "v"),
        edge(leaves[1], "u"),
        edge(leaves[2], "u"),
        edge(leaves[3], "w"),
        edge(leaves[4], "w"),
        edge(leaves[5], "v"),
        edge(leaves[6], "y"),
        edge(leaves[7], "z"),
        admixture_edge("M", "x", "y")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "y", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  one_admixture_35 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 315
      leaf_permutations <- symmetry_7_IV(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "w", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "x"),
        edge("u", "z"),
        edge("v", "z"),
        edge("w", "y"),
        edge(leaves[1], "u"),
        edge(leaves[2], "u"),
        edge(leaves[3], "v"),
        edge(leaves[4], "v"),
        edge(leaves[5], "w"),
        edge(leaves[6], "w"),
        edge(leaves[7], "M"),
        admixture_edge("M", "x", "y")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "y", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  one_admixture_36 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 315
      leaf_permutations <- symmetry_7_IV(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "w", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "x"),
        edge("u", "z"),
        edge("v", "z"),
        edge("w", "M"),
        edge(leaves[1], "u"),
        edge(leaves[2], "u"),
        edge(leaves[3], "v"),
        edge(leaves[4], "v"),
        edge(leaves[5], "w"),
        edge(leaves[6], "w"),
        edge(leaves[7], "y"),
        admixture_edge("M", "x", "y")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "y", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  one_admixture_37 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 315
      leaf_permutations <- symmetry_7_IV(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "w", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "M"),
        edge("u", "z"),
        edge("v", "z"),
        edge("w", "y"),
        edge(leaves[1], "u"),
        edge(leaves[2], "u"),
        edge(leaves[3], "v"),
        edge(leaves[4], "v"),
        edge(leaves[5], "w"),
        edge(leaves[6], "w"),
        edge(leaves[7], "x"),
        admixture_edge("M", "x", "y")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "y", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  one_admixture_38 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 1260
      leaf_permutations <- symmetry_7_II(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "w", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "x"),
        edge("u", "y"),
        edge("v", "u"),
        edge("w", "v"),
        edge(leaves[1], "z"),
        edge(leaves[2], "z"),
        edge(leaves[3], "u"),
        edge(leaves[4], "v"),
        edge(leaves[5], "w"),
        edge(leaves[6], "w"),
        edge(leaves[7], "M"),
        admixture_edge("M", "x", "y")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "y", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  one_admixture_39 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 1260
      leaf_permutations <- symmetry_7_II(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "w", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "x"),
        edge("u", "z"),
        edge("v", "u"),
        edge("w", "M"),
        edge(leaves[1], "w"),
        edge(leaves[2], "w"),
        edge(leaves[3], "v"),
        edge(leaves[4], "v"),
        edge(leaves[5], "u"),
        edge(leaves[6], "z"),
        edge(leaves[7], "y"),
        admixture_edge("M", "x", "y")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "y", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  one_admixture_40 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 1260
      leaf_permutations <- symmetry_7_II(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "w", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "x"),
        edge("u", "M"),
        edge("v", "u"),
        edge("w", "v"),
        edge(leaves[1], "z"),
        edge(leaves[2], "z"),
        edge(leaves[3], "w"),
        edge(leaves[4], "w"),
        edge(leaves[5], "v"),
        edge(leaves[6], "u"),
        edge(leaves[7], "y"),
        admixture_edge("M", "x", "y")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "y", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  one_admixture_41 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 630
      leaf_permutations <- symmetry_7_V(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "w", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "x"),
        edge("u", "z"),
        edge("v", "u"),
        edge("w", "u"),
        edge(leaves[1], "v"),
        edge(leaves[2], "v"),
        edge(leaves[3], "w"),
        edge(leaves[4], "w"),
        edge(leaves[5], "z"),
        edge(leaves[6], "M"),
        edge(leaves[7], "y"),
        admixture_edge("M", "x", "y")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "y", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  one_admixture_42 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 315
      leaf_permutations <- symmetry_7_IV(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "w", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "M"),
        edge("u", "z"),
        edge("v", "u"),
        edge("w", "u"),
        edge(leaves[1], "v"),
        edge(leaves[2], "v"),
        edge(leaves[3], "w"),
        edge(leaves[4], "w"),
        edge(leaves[5], "x"),
        edge(leaves[6], "y"),
        edge(leaves[7], "z"),
        admixture_edge("M", "x", "y")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "y", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  one_admixture_43 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 1260
      leaf_permutations <- symmetry_7_II(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "w", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "M"),
        edge("u", "z"),
        edge("v", "u"),
        edge("w", "v"),
        edge(leaves[1], "x"),
        edge(leaves[2], "y"),
        edge(leaves[3], "w"),
        edge(leaves[4], "w"),
        edge(leaves[5], "v"),
        edge(leaves[6], "u"),
        edge(leaves[7], "z"),
        admixture_edge("M", "x", "y")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "y", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  one_admixture_44 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 2520
      leaf_permutations <- symmetry_7_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "w", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "x"),
        edge("u", "z"),
        edge("v", "u"),
        edge("w", "v"),
        edge(leaves[1], "w"),
        edge(leaves[2], "w"),
        edge(leaves[3], "v"),
        edge(leaves[4], "u"),
        edge(leaves[5], "z"),
        edge(leaves[6], "M"),
        edge(leaves[7], "y"),
        admixture_edge("M", "x", "y")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "y", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  one_admixture_45 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 630
      leaf_permutations <- symmetry_7_III(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "w", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "M"),
        edge("u", "z"),
        edge("v", "z"),
        edge("w", "v"),
        edge(leaves[1], "x"),
        edge(leaves[2], "y"),
        edge(leaves[3], "u"),
        edge(leaves[4], "u"),
        edge(leaves[5], "w"),
        edge(leaves[6], "w"),
        edge(leaves[7], "v"),
        admixture_edge("M", "x", "y")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "y", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  one_admixture_46 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 1260
      leaf_permutations <- symmetry_7_II(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "w", "M")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "x"),
        edge("u", "z"),
        edge("v", "z"),
        edge("w", "v"),
        edge(leaves[1], "u"),
        edge(leaves[2], "u"),
        edge(leaves[3], "w"),
        edge(leaves[4], "w"),
        edge(leaves[5], "v"),
        edge(leaves[6], "M"),
        edge(leaves[7], "y"),
        admixture_edge("M", "x", "y")
      ))
      admixtures <- admixture_proportions(c(
        admix_props("M", "x", "y", "a")
      ))
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  }
)

symmetry_7_I <- function(leaves) {
  output <- character(7)
  all_permutations <- make_permutations(leaves)
  for (candidate in all_permutations) {
    sorted <- sort(candidate)
    num <- numeric(0)
    for (j in seq(1, 7)) {
      num[j] <- match(candidate[j], sorted)
    }
    if (num[1] < num[2]) {
      output <- rbind(output, candidate)
    }
  }
  output <- output[-1, ]
  output
}

symmetry_7_II <- function(leaves) {
  output <- character(7)
  all_permutations <- make_permutations(leaves)
  for (candidate in all_permutations) {
    sorted <- sort(candidate)
    num <- numeric(0)
    for (j in seq(1, 7)) {
      num[j] <- match(candidate[j], sorted)
    }
    if (num[1] < num[2] &&
        num[3] < num[4]) {
      output <- rbind(output, candidate)
    }
  }
  output <- output[-1, ]
  output
}

symmetry_7_III <- function(leaves) {
  output <- character(7)
  all_permutations <- make_permutations(leaves)
  for (candidate in all_permutations) {
    sorted <- sort(candidate)
    num <- numeric(0)
    for (j in seq(1, 7)) {
      num[j] <- match(candidate[j], sorted)
    }
    if (num[1] < num[2] &&
        num[3] < num[4] &&
        num[5] < num[6]) {
      output <- rbind(output, candidate)
    }
  }
  output <- output[-1, ]
  output
}

symmetry_7_IV <- function(leaves) {
  output <- character(7)
  all_permutations <- make_permutations(leaves)
  for (candidate in all_permutations) {
    sorted <- sort(candidate)
    num <- numeric(0)
    for (j in seq(1, 7)) {
      num[j] <- match(candidate[j], sorted)
    }
    if (num[1] < num[2] &&
        num[3] < num[4] &&
        num[5] < num[6] &&
        num[1] < num[2] &&
        num[1] + 10*num[2] < num[3] + 10*num[4]) {
      output <- rbind(output, candidate)
    }
  }
  output <- output[-1, ]
  output
}

symmetry_7_V <- function(leaves) {
  output <- character(7)
  all_permutations <- make_permutations(leaves)
  for (candidate in all_permutations) {
    sorted <- sort(candidate)
    num <- numeric(0)
    for (j in seq(1, 7)) {
      num[j] <- match(candidate[j], sorted)
    }
    if (num[1] < num[2] &&
        num[3] < num[4] &&
        num[1] + 10*num[2] < num[3] + 10*num[4]) {
      output <- rbind(output, candidate)
    }
  }
  output <- output[-1, ]
  output
}

#' Eight leaves trees.
#' 
#' Kind of obsolete since the introduction of \code{\link{all_graphs}}.
#' A comprehensive listing of all the four unrooted trees with eight leaves.
#' The position of the root can be moved later with the function
#' \code{\link{make_an_outgroup}}.
#' 
#' @format A list of functions on eight leaves.
#'         The outputs of these functions are \code{\link{agraph}} objects.
#'
#' @family graphs
#'
#' @seealso \code{\link{all_graphs}}
#' @seealso \code{\link{make_permutations}}
#' @seealso \code{\link{fit_permutations_and_graphs}}
#' @seealso \code{\link{fit_graph_list}}
#' @seealso \code{\link{add_a_leaf}}
#' @seealso \code{\link{add_an_admixture}}
#' @seealso \code{\link{add_an_admixture2}}
#' @seealso \code{\link{make_an_outgroup}}
#'
#' @examples
#' \donttest{
#' # While the usage of this function is pretty self-explanatory, let's plot all the graphs
#' # just for browsing.
#' for (i in seq(1, length(eight_leaves_trees))) {
#'   graph <- eight_leaves_trees[[i]](c("A", "B", "C", "D", "E", "F", "G", "H"))
#'   # This is how you include quotation marks in strings by the way:
#'   title <- paste("eight_leaves_trees[[", i,
#'                  "]](c(\"A\", \"B\", \"C\", \"D\", \"E\", \"F\", \"G\", \"H\"))", sep = "")
#'   plot(graph, color = "brown", title = title)
#' }
#' }
#'
#' @export
eight_leaves_trees <- list(
  tree_1 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 5040
      leaf_permutations <- symmetry_8_I(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "w")
      edges <- parent_edges(c(
        edge("x", "y"),
        edge("y", "z"),
        edge("z", "R"),
        edge("u", "R"),
        edge("v", "u"),
        edge("w", "v"),
        edge(leaves[1], "x"),
        edge(leaves[2], "x"),
        edge(leaves[3], "y"),
        edge(leaves[4], "z"),
        edge(leaves[5], "u"),
        edge(leaves[6], "v"),
        edge(leaves[7], "w"),
        edge(leaves[8], "w")
      ))
      admixtures <- NULL
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  tree_2 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 2520
      leaf_permutations <- symmetry_8_II(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "w")
      edges <- parent_edges(c(
        edge("x", "y"),
        edge("y", "z"),
        edge("z", "R"),
        edge("u", "v"),
        edge("v", "R"),
        edge("w", "v"),
        edge(leaves[1], "x"),
        edge(leaves[2], "x"),
        edge(leaves[3], "y"),
        edge(leaves[4], "z"),
        edge(leaves[5], "u"),
        edge(leaves[6], "u"),
        edge(leaves[7], "w"),
        edge(leaves[8], "w")
      ))
      admixtures <- NULL
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  tree_3 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 2520
      leaf_permutations <- symmetry_8_III(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "w")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "z"),
        edge("z", "u"),
        edge("u", "R"),
        edge("v", "u"),
        edge("w", "v"),
        edge(leaves[1], "x"),
        edge(leaves[2], "x"),
        edge(leaves[3], "y"),
        edge(leaves[4], "y"),
        edge(leaves[5], "z"),
        edge(leaves[6], "v"),
        edge(leaves[7], "w"),
        edge(leaves[8], "w")
      ))
      admixtures <- NULL
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  },
  
  tree_3 <- function(leaves, permutations = FALSE) {
    if (permutations == TRUE) {
      P <- 315
      leaf_permutations <- symmetry_8_IV(leaves)
    } else {
      P <- 1
      leaf_permutations <- rbind(leaves)
    }
    result <- vector(mode = "list", length = P)
    for (j in seq(1, P)) {
      leaves <- leaf_permutations[j, ]
      inner_nodes <- c("R", "x", "y", "z", "u", "v", "w")
      edges <- parent_edges(c(
        edge("x", "R"),
        edge("y", "R"),
        edge("z", "x"),
        edge("u", "x"),
        edge("v", "y"),
        edge("w", "y"),
        edge(leaves[1], "z"),
        edge(leaves[2], "z"),
        edge(leaves[3], "u"),
        edge(leaves[4], "u"),
        edge(leaves[5], "v"),
        edge(leaves[6], "v"),
        edge(leaves[7], "w"),
        edge(leaves[8], "w")
      ))
      admixtures <- NULL
      result[[j]] <- agraph(leaves, inner_nodes, edges, admixtures)  
    }
    if (permutations == FALSE) {result <- result[[1]]}
    result
  }
)

symmetry_8_I <- function(leaves) {
  output <- character(8)
  all_permutations <- make_permutations(leaves)
  for (candidate in all_permutations) {
    sorted <- sort(candidate)
    num <- numeric(0)
    for (j in seq(1, 8)) {
      num[j] <- match(candidate[j], sorted)
    }
    if (num[1] < num[2] &&
        num[4] < num[5] &&
        num[7] < num[8]) {
      output <- rbind(output, candidate)
    }
  }
  output <- output[-1, ]
  output
}

symmetry_8_II <- function(leaves) {
  output <- character(8)
  all_permutations <- make_permutations(leaves)
  for (candidate in all_permutations) {
    sorted <- sort(candidate)
    num <- numeric(0)
    for (j in seq(1, 8)) {
      num[j] <- match(candidate[j], sorted)
    }
    if (num[1] < num[2] &&
        num[5] < num[6] &&
        num[7] < num[8] &&
        num[5] + 10*num[6] < num[7] + 10*num[8]) {
      output <- rbind(output, candidate)
    }
  }
  output <- output[-1, ]
  output
}

symmetry_8_III <- function(leaves) {
  output <- character(8)
  all_permutations <- make_permutations(leaves)
  for (candidate in all_permutations) {
    sorted <- sort(candidate)
    num <- numeric(0)
    for (j in seq(1, 8)) {
      num[j] <- match(candidate[j], sorted)
    }
    if (num[1] < num[2] &&
        num[3] < num[4] &&
        num[5] < num[6] &&
        num[7] < num[8]) {
      output <- rbind(output, candidate)
    }
  }
  output <- output[-1, ]
  output
}

symmetry_8_IV <- function(leaves) {
  output <- character(8)
  all_permutations <- make_permutations(leaves)
  for (candidate in all_permutations) {
    sorted <- sort(candidate)
    num <- numeric(0)
    for (j in seq(1, 8)) {
      num[j] <- match(candidate[j], sorted)
    }
    if (num[1] < num[2] &&
        num[3] < num[4] &&
        num[1] + 10*num[2] < num[3] + 10*num[4] &&
        num[5] < num[6] &&
        num[7] < num[8] &&
        num[5] + 10*num[6] < num[7] + 10*num[8] &&
        num[1] + 10*num[2] + 100*num[3] + 1000*num[4] < num[5] + 10*num[6] + 100*num[7] + 1000*num[8]) {
      output <- rbind(output, candidate)
    }
  }
  output <- output[-1, ]
  output
}

#' Fit lots of graphs to data.
#'
#' Combines a list of (population) permutations and a list of graph topologies
#' to a big list of graphs, then fits those graphs to given data using parallel
#' computation. This function needs \code{doParallel}, \code{foreach} and
#' \code{parallel} installed.
#'   
#' @param data          The data table.
#' @param permutations  List of population permutations.
#' @param graphs        List of functions for producing graphs.
#' @param cores         Number of cores used.
#' 
#' @return A list of \code{\link{fast_fit}} results.
#'
#' @seealso \code{\link{make_permutations}}
#' @seealso \code{\link{four_leaves_graphs}}
#' @seealso \code{\link{five_leaves_graphs}}
#' @seealso \code{\link{six_leaves_graphs}}
#' @seealso \code{\link{seven_leaves_graphs}}
#' @seealso \code{\link{eight_leaves_trees}}
#' @seealso \code{\link{fit_graph_list}}
#' 
#' @examples
#' \donttest{
#' # Let's experiment by fitting all the graphs with five leaves and at most one admixture
#' # event to a five population subset of the bear data. Note that with three data rows only
#' # we do wisely by not concluding too much about the actual bear family tree; this is to
#' # illustrate the function usage only!
#' 
#' data(bears)
#' data <- bears[16:18, ]
#' print(data)
#' permutations <- make_permutations(c("PB", "BLK", "Sweden", "Denali", "Kenai"))
#' graphs <- five_leaves_graphs
#' 
#' # We go with one core only as I don't know what kind of machine you are using.
#' 
#' fitted_graphs <- fit_permutations_and_graphs(data, permutations, graphs, 1)
#' 
#' # Now sort the fitted objects by best_error and see how the best graph looks like.
#' 
#' errors <- sapply(fitted_graphs, function(x) x$best_error)
#' best_graphs <- fitted_graphs[order(errors)]
#' plot(best_graphs[[1]]$graph, color = "goldenrod", title = best_graphs[[1]]$best_error)
#'
#' # The same value for best_error actually occurs in the list 152 times because of our
#' # unsufficient data.
#' }
#'
#' @export
fit_permutations_and_graphs <- function(data, permutations, graphs, cores) {
  if (!requireNamespace("doParallel", quietly = TRUE)) {
    stop("This function needs the package 'doParallel' to be installed.")
  }
  if (!requireNamespace("foreach", quietly = TRUE)) {
    stop("This function needs the package 'foreach' to be installed.")
  }
  if (!requireNamespace("parallel", quietly = TRUE)) {
    stop("This function needs the package 'parallel' to be installed.")
  }
  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)
  # The following is a mess because I don't want to attach packages within my own, so I have to use
  # the binary operators as they were normal functions.
  i <- seq(1, length(permutations))
  j <- seq(1, length(graphs))
  foreach::'%dopar%'(foreach::'%:%'(foreach::foreach(i = i, .combine = c, .packages = "admixturegraph"),
                                    foreach::foreach(j = j, .packages = "admixturegraph")), {
    permutation <- permutations[[i]]
    graph_function <- graphs[[j]]
    graph <- graph_function(permutation)
    result <- fast_fit(filter_on_leaves(data, graph), graph)
    return(result)
  })
}

#' Fit lots of graphs to data.
#'
#' Fits a list of graphs to given data using parallel computation. This function
#' needs \code{doParallel}, \code{foreach} and \code{parallel} installed.
#'   
#' @param data          The data table.
#' @param graphs        List of graphs.
#' @param cores         Number of cores used.
#' 
#' @return A list of \code{\link{fast_fit}} results.
#'
#' @seealso \code{\link{all_graphs}}
#' @seealso \code{\link{fit_permutations_and_graphs}}
#' 
#' @export
fit_graph_list <- function(data, graphs, cores) {
  if (!requireNamespace("doParallel", quietly = TRUE)) {
    stop("This function needs the package 'doParallel' to be installed.")
  }
  if (!requireNamespace("foreach", quietly = TRUE)) {
    stop("This function needs the package 'foreach' to be installed.")
  }
  if (!requireNamespace("parallel", quietly = TRUE)) {
    stop("This function needs the package 'parallel' to be installed.")
  }
  cl <- parallel::makeCluster(cores)
  doParallel::registerDoParallel(cl)
  # The following is a mess because I don't want to attach packages within my own, so I have to use
  # the binary operators as they were normal functions.
  j <- seq(1, length(graphs))
  foreach::'%dopar%'(foreach::foreach(j = j, .packages = "admixturegraph"), {
    graph <- graphs[[j]]
    result <- fast_fit(filter_on_leaves(data, graph), graph)
    return(result)
  })
}

#' Adds a new leaf to a graph.
#' 
#' Given an admixture graph, selects an edge and branches off a new edge ending at a new leaf. 
#' 
#' @param graph      An admixture graph.
#' @param leaf_name  A name for the new leaf.
#' @param outgroup   An optional parameter for the preferred outgroup, which can be the new leaf.
#' 
#' @return A list of graphs made by adding a new leaf to the input graph. The list has no
#'         duplicate elements.
#'         
#' @seealso \code{\link{all_graphs}}
#' @seealso \code{\link{fit_permutations_and_graphs}}
#' @seealso \code{\link{fit_graph_list}}
#' @seealso \code{\link{add_an_admixture}}
#' @seealso \code{\link{add_an_admixture2}}
#' @seealso \code{\link{make_an_outgroup}}
#' 
#' @examples
#' \donttest{
#' # Take a look at how much trees there are: 
#' leaves <- c("1", "2")
#' inner_nodes <- c("R")
#' edges <- parent_edges(c(edge("1", "R"), edge("2", "R")))
#' admixtures <- NULL
#' Lambda <- agraph(leaves, inner_nodes, edges, admixtures)
#' set <- list(Lambda)
#' for (i in seq(1, 6)) {
#'   new_set <- list()
#'   for (tree in set) {
#'     new_set <- c(new_set, add_a_leaf(tree, paste(i + 2)))
#'   }
#'   set <- new_set
#'   cat("There are ")
#'   cat(length(set))
#'   cat(" different trees with ")
#'   cat(i + 2)
#'   cat(" labeled leaves.")
#'   cat("\n")
#' }
#' # In general, there are 1*3*5*...*(2n - 5) different trees with n labeled leaves
#' # (A001147 in the online encyclopedia of integer sequences).
#' # Exhaustive search through the graph space is hard!
#' }
#' 
#' @export
add_a_leaf <- function(graph, leaf_name, outgroup = "") {
  graph_list <- list()
  broken_graph <- break_graph(graph)
  inner_name <- paste("inner", leaf_name, sep = "_")
  root <- broken_graph$root
  forgive <- TRUE
  for (i in seq(1, length(broken_graph$edges))) {
    leaves <- c(broken_graph$leaves, leaf_name)
    inner_nodes <- c(broken_graph$inner_nodes, inner_name)
    new_edges <- broken_graph$edges
    chosen <- new_edges[[i]]
    edge_argument <- character(0)
    new_admixtures <- broken_graph$admixtures
    admix_argument <- character(0)
    new_edges[[i]] <- NULL
    new_edges[[length(new_edges) + 1]] <- c(leaf_name, inner_name)
    new_edges[[length(new_edges) + 1]] <- c(chosen[1], inner_name)
    new_edges[[length(new_edges) + 1]] <- c(inner_name, chosen[2])
    for (j in seq(1, length(new_edges))) {
      edge_argument <- c(edge_argument, edge(new_edges[[j]][1], new_edges[[j]][2]))
    }
    if (length(broken_graph$admixtures) > 0) {
      for (k in seq(1, length(new_admixtures))) {
        edge_argument <- c(edge_argument, admixture_edge(new_admixtures[[k]][1], new_admixtures[[k]][2],
                                                         new_admixtures[[k]][3]))
        admix_argument <- c(admix_argument, admix_props(new_admixtures[[k]][1], new_admixtures[[k]][2],
                                                        new_admixtures[[k]][3], new_admixtures[[k]][4]))
      }
    }
    edges <- parent_edges(edge_argument)
    if (length(broken_graph$admixtures) > 0) {
      admixtures <- admixture_proportions(admix_argument)
    } else {
      admixtures <- NULL
    }
    skip <- FALSE
    if (chosen[2] == root) {
      if (forgive == TRUE) {
        forgive <- FALSE
      } else {
        skip <- TRUE
      }
    }
    if (skip == FALSE) {
      new_graph <- agraph(leaves, inner_nodes, edges, admixtures)
      graph_list[[length(graph_list) + 1]] <- make_an_outgroup(new_graph, outgroup)
    }
  }
  if (length(broken_graph$admixtures) > 0) {
    for (i in seq(1, length(broken_graph$admixtures))) {
      leaves <- c(broken_graph$leaves, leaf_name)
      inner_nodes <- c(broken_graph$inner_nodes, inner_name)
      new_edges <- broken_graph$edges
      edge_argument <- character(0)
      new_admixtures <- broken_graph$admixtures
      chosen <- new_admixtures[[i]]
      admix_argument <- character(0)
      new_admixtures[[i]] <- NULL
      new_edges[[length(new_edges) + 1]] <- c(leaf_name, inner_name)
      new_edges[[length(new_edges) + 1]] <- c(inner_name, chosen[2])
      new_admixtures[[length(new_admixtures) + 1]] <- c(chosen[1], inner_name, chosen[3], chosen[4])
      for (j in seq(1, length(new_edges))) {
        edge_argument <- c(edge_argument, edge(new_edges[[j]][1], new_edges[[j]][2]))
      }
      for (k in seq(1, length(new_admixtures))) {
        edge_argument <- c(edge_argument, admixture_edge(new_admixtures[[k]][1], new_admixtures[[k]][2],
                                                         new_admixtures[[k]][3]))
        admix_argument <- c(admix_argument, admix_props(new_admixtures[[k]][1], new_admixtures[[k]][2],
                                                        new_admixtures[[k]][3], new_admixtures[[k]][4]))
      }
      edges <- parent_edges(edge_argument)
      admixtures <- admixture_proportions(admix_argument)
      new_graph <- agraph(leaves, inner_nodes, edges, admixtures)
      graph_list[[length(graph_list) + 1]] <- make_an_outgroup(new_graph, outgroup)
      new_edges <- broken_graph$edges
      edge_argument <- character(0)
      new_admixtures <- broken_graph$admixtures
      chosen <- new_admixtures[[i]]
      admix_argument <- character(0)
      new_admixtures[[i]] <- NULL
      new_edges[[length(new_edges) + 1]] <- c(leaf_name, inner_name)
      new_edges[[length(new_edges) + 1]] <- c(inner_name, chosen[3])
      new_admixtures[[length(new_admixtures) + 1]] <- c(chosen[1], chosen[2], inner_name, chosen[4])
      for (j in seq(1, length(new_edges))) {
        edge_argument <- c(edge_argument, edge(new_edges[[j]][1], new_edges[[j]][2]))
      }
      for (k in seq(1, length(new_admixtures))) {
        edge_argument <- c(edge_argument, admixture_edge(new_admixtures[[k]][1], new_admixtures[[k]][2],
                                                         new_admixtures[[k]][3]))
        admix_argument <- c(admix_argument, admix_props(new_admixtures[[k]][1], new_admixtures[[k]][2],
                                                        new_admixtures[[k]][3], new_admixtures[[k]][4]))
      }
      edges <- parent_edges(edge_argument)
      admixtures <- admixture_proportions(admix_argument)
      new_graph <- agraph(leaves, inner_nodes, edges, admixtures)
      graph_list[[length(graph_list) + 1]] <- make_an_outgroup(new_graph, outgroup)
    }  
  }
  return(graph_list)
}

#' Adds a new admixture event to a graph.
#' 
#' Given an admixture graph, selects a child edge and a parent edges and adds a new edge from the 
#' parent edge to the child edge with an admixture event, if possible.
#' Thus, the resulting graph is an extension of the input graph in the sense that erasing one of
#' the admixture edges (the new one) we get the original admixture graph. However, we know that 
#' in practice when fitting data to admixture graphs, the best graph with \eqn{k} admixture events
#' is not always an extension like that from the best graph with \eqn{k - 1} admixture events.
#' For a more relaxed way of adding a new admixture event, try \code{\link{add_an_admixture2}}.
#' 
#' @param graph                    An admixture graph.
#' @param admixture_variable_name  A name for the new admixture proportion.
#' @param labels_matter            When \code{FALSE} (the default value), we consider two admixture
#'                                 graphs similar when they have the same topology but permuted
#'                                 admixture proportion names. When \code{TRUE}, the already existing
#'                                 admixture events and the edges leading to them are considered
#'                                 labeled.
#' @param outgroup                 An optional parameter for the preferred outgroup.
#' 
#' @return A list of graphs made by adding a new admixture event to the input graph. The list has
#'         no duplicate elements (what that means depends on the value of \code{labels_matter}).
#' 
#' @seealso \code{\link{all_graphs}}
#' @seealso \code{\link{fit_permutations_and_graphs}}
#' @seealso \code{\link{fit_graph_list}}
#' @seealso \code{\link{add_a_leaf}}
#' @seealso \code{\link{add_an_admixture2}}
#' @seealso \code{\link{make_an_outgroup}}
#' 
#' @examples
#' \donttest{
#' # To illustrate what the parameter labels_matter does, consider the following graph:
#' 
#' leaves <- c("A", "B", "C")
#' inner_nodes <- c("R", "x", "y", "M")
#' edges <- parent_edges(c(edge("x", "R"),
#'                         edge("y", "R"),
#'                         edge("A", "x"),
#'                         edge("B", "M"),
#'                         edge("C", "y"),
#'                         admixture_edge("M", "x", "y")))
#' admixtures <- admixture_proportions(c(admix_props("M", "x", "y", "p")))
#' graph <- agraph(leaves, inner_nodes, edges, admixtures)
#' plot(graph, show_admixture_labels = TRUE, title = "graph")
#' 
#' # There are 15 ways this graph can be extended to a graph with two admixture events by
#' # adding an admixture edge, as can be seeing by having the program explicitly construct
#' # all the cases:
#' 
#' short_list <- add_an_admixture(graph, "q")
#' print(length(short_list))
#' 
#' # However, maybe we already have a specific historical event in mind corresponding to the
#' # original admixture event in graph, or a fixed value for the admixture proportion p.
#' # Then, for example, it makes a difference to us whether we consider the possibility of
#' # another admixture event occurring before that event,
#' 
#' leaves <- c("A", "B", "C")
#' inner_nodes <- c("R", "x", "y", "z", "M", "N")
#' edges <- parent_edges(c(edge("x", "R"),
#'                         edge("z", "R"),
#'                         edge("y", "z"),
#'                         edge("A", "x"),
#'                         edge("B", "M"),
#'                         edge("C", "y"),
#'                         admixture_edge("N", "x", "z"),
#'                         admixture_edge("M", "N", "y")))
#' admixtures <- admixture_proportions(c(admix_props("N", "x", "z", "q"),
#'                                       admix_props("M", "N", "y", "p")))
#' example1 <- agraph(leaves, inner_nodes, edges, admixtures)
#' plot(example1, show_admixture_labels = TRUE, title = "example 1")
#' 
#' # or after that event,
#' 
#' leaves <- c("A", "B", "C")
#' inner_nodes <- c("R", "x", "y", "z", "M", "N")
#' edges <- parent_edges(c(edge("x", "R"),
#'                         edge("y", "R"),
#'                         edge("z", "y"),
#'                         edge("A", "x"),
#'                         edge("B", "N"),
#'                         edge("C", "z"),
#'                         admixture_edge("M", "x", "y"),
#'                         admixture_edge("N", "M", "z")))
#' admixtures <- admixture_proportions(c(admix_props("M", "x", "y", "p"),
#'                                       admix_props("N", "M", "z", "q")))
#' example2 <- agraph(leaves, inner_nodes, edges, admixtures)
#' plot(example2, show_admixture_labels = TRUE, title = "example 2")
#' 
#' # even though as (acyclic) directed graphs with labeled leaves example 1
#' # and example 2 are isomorphic.
#' # Counting cases like that dirrerent, we get 21 possible extensions:
#'
#' long_list <- add_an_admixture(graph, "q", labels_matter = TRUE)
#' print(length(long_list))
#' }
#' 
#' @export
add_an_admixture <- function(graph, admixture_variable_name, labels_matter = FALSE, outgroup = "") {
  graph_list <- list()
  broken_graph <- break_graph(graph)
  # We might have to choose a different root after adding the admixture, so we start by removing
  # the original root and treating only the edges colliding in an admix event as directed and the
  # rest as undirected.
  leaves <- broken_graph$leaves
  original_inner_nodes <- broken_graph$inner_nodes
  original_edges <- broken_graph$edges
  original_admixtures <- broken_graph$admixtures
  root <- broken_graph$root
  for (i in seq(1, length(original_inner_nodes))) {
    node <- original_inner_nodes[i]
    if (node == root) {
      original_inner_nodes <- original_inner_nodes[-i]
      break
    }
  }
  memory <- ""
  to_be_deleted <- nchar(0)
  for (i in seq(1, length(original_edges))) {
    edge <- original_edges[[i]]
    if (edge[2] == root) {
      if (nchar(memory) == 0) {
        memory <- edge[1]
        to_be_deleted <- i
      } else {
        original_edges[[i]] <- c(memory, edge[1])
      }
    }
  }
  original_edges[[to_be_deleted]] <- NULL
  original_directed_edges <- list()
  for (admixture in original_admixtures) {
    original_directed_edges[[length(original_directed_edges) + 1]] <- c(admixture[2], admixture[1])
    original_directed_edges[[length(original_directed_edges) + 1]] <- c(admixture[3], admixture[1])
  }
  inner_name <- paste("inner", admixture_variable_name, sep = "_")
  admix_name <- paste("admix", admixture_variable_name, sep = "_")
  inner_nodes <- c(original_inner_nodes, inner_name, admix_name)
  for (i in seq(1, length(original_edges))) {
    for (j in seq(1, length(original_edges))) {
      if (i != j) {
        # From edge to edge:
        # Weeding out half of the duplicated cases.
        default_problem <- FALSE
        if (length(intersect(original_edges[[i]], original_edges[[j]])) > 0) {
          common_node <- intersect(original_edges[[i]], original_edges[[j]])
          for (k in seq(1, i)) {
            if (length(intersect(common_node, original_edges[[k]])) > 0 && k < i && k != j) {
              # All three undirected -> compare indices of i and k.
              default_problem <- TRUE
            }
          }
          if (length(original_directed_edges) > 0) {
            for (k in seq(1, length(original_directed_edges))) {
              if (original_directed_edges[[k]][1] == common_node) {
                # The only directed not j -> always choose the directed i.
                default_problem <- TRUE
              }
            }
          }
        }
        # Optionally weeding out half of the cases that are indistinguishable from existing admixtures.
        if (labels_matter == FALSE) {
          if (length(original_directed_edges) > 0) {
            for (k in seq(1, length(original_directed_edges))) {
              if (length(intersect(original_directed_edges[[k]][1], original_edges[[i]])) > 0 &&
                  length(intersect(original_directed_edges[[k]][2], original_edges[[j]])) > 0) {
                default_problem <- TRUE
              }
            }
          }
        }
        # Direction [2]:
        admixtures <- original_admixtures
        admixtures[[length(admixtures) + 1]] <- c(admix_name, inner_name, original_edges[[j]][1], admixture_variable_name)
        directed_edges <- original_directed_edges
        directed_edges[[length(directed_edges) + 1]] <- c(inner_name, admix_name)
        directed_edges[[length(directed_edges) + 1]] <- c(original_edges[[j]][1], admix_name)
        edges <- original_edges
        edges[[i]] <- c(inner_name, original_edges[[i]][1])
        edges[[j]] <- c(inner_name, original_edges[[i]][2])
        edges[[length(edges) + 1]] <- c(admix_name, original_edges[[j]][2])
        starting_directed_edge <- c(admix_name, original_edges[[j]][2])
        flow_result <- flow(leaves, edges, directed_edges, starting_directed_edge, default_problem)
        if (flow_result$problem == FALSE) {
          graph <- root_graph(leaves, inner_nodes, flow_result$edges, flow_result$directed_edges, admixtures, outgroup)
          graph_list[[length(graph_list) + 1]] <- graph
        }
        # Direction [1]:
        admixtures <- original_admixtures
        admixtures[[length(admixtures) + 1]] <- c(admix_name, inner_name, original_edges[[j]][2], admixture_variable_name)
        directed_edges <- original_directed_edges
        directed_edges[[length(directed_edges) + 1]] <- c(inner_name, admix_name)
        directed_edges[[length(directed_edges) + 1]] <- c(original_edges[[j]][2], admix_name)
        edges <- original_edges
        edges[[i]] <- c(inner_name, original_edges[[i]][1])
        edges[[j]] <- c(inner_name, original_edges[[i]][2])
        edges[[length(edges) + 1]] <- c(admix_name, original_edges[[j]][1])
        starting_directed_edge <- c(admix_name, original_edges[[j]][1])
        flow_result <- flow(leaves, edges, directed_edges, starting_directed_edge, default_problem)
        if (flow_result$problem == FALSE) {
          graph <- root_graph(leaves, inner_nodes, flow_result$edges, flow_result$directed_edges, admixtures, outgroup)
          graph_list[[length(graph_list) + 1]] <- graph
        }
      }
    }
    if (length(original_directed_edges) > 0) {
      for (j in seq(1, length(original_directed_edges))) {
        # From edge to directed edge:
        # Weeding out half of the duplicated cases.
        default_problem <- FALSE
        if (length(intersect(original_edges[[i]], original_directed_edges[[j]])) > 0) {
          common_node <- intersect(original_edges[[i]], original_directed_edges[[j]])
          for (k in seq(1, i)) {
            if (length(intersect(common_node, original_edges[[k]])) > 0 && k < i) {
              # The only directed j -> compare indices of i and k.
              default_problem <- TRUE
            }
          }
          for (k in seq(1, length(original_directed_edges))) {
            if (length(intersect(common_node, original_directed_edges[[k]])) > 0 && k != j) {
              # The only undirected not j, diverging at common node -> always choose the directed i.
              default_problem <- TRUE
            }
          }
        }
        # Optionally weeding out half of the cases that are indistinguishable from existing admixtures.
        if (labels_matter == FALSE) {
          for (k in seq(1, length(original_directed_edges))) {
            if (length(intersect(original_directed_edges[[k]][1], original_edges[[i]])) > 0 &&
                original_directed_edges[[k]][2] == original_directed_edges[[j]][1]) {
              default_problem <- TRUE
            }
          }
          if (length(intersect(original_edges[[i]], original_directed_edges[[j]])) == 0) {
            for (k in seq(1, j)) {
              if (length(intersect(original_directed_edges[[k]][1], original_edges[[i]])) > 0 &&
                  original_directed_edges[[k]][2] == original_directed_edges[[j]][2] && k < j) {
                for (l in seq(1, length(original_edges))) {
                  if (length(intersect(original_edges[[l]], original_directed_edges[[j]][1])) &&
                      length(intersect(original_edges[[l]], original_directed_edges[[k]][1]))) {
                    default_problem <- TRUE
                  }
                }
              }
            }
          }
        }
        admixtures <- original_admixtures
        for (k in seq(1, length(admixtures))) {
          admixture <- admixtures[[k]]
          if (admixture[1] == original_directed_edges[[j]][2]) {
            if (admixture[2] == original_directed_edges[[j]][1]) {
              admixture[2] <- admix_name
            }
            if (admixture[3] == original_directed_edges[[j]][1]) {
              admixture[3] <- admix_name
            }
            admixtures[[k]] <- admixture
          }
        }
        admixtures[[length(admixtures) + 1]] <- c(admix_name, inner_name, original_directed_edges[[j]][1], admixture_variable_name)
        directed_edges <- original_directed_edges
        directed_edges[[j]] <- c(inner_name, admix_name)
        directed_edges[[length(directed_edges) + 1]] <- c(original_directed_edges[[j]][1], admix_name)
        directed_edges[[length(directed_edges) + 1]] <- c(admix_name, original_directed_edges[[j]][2])
        edges <- original_edges
        edges[[i]] <- c(inner_name, original_edges[[i]][1])
        edges[[length(edges) + 1]] <- c(inner_name, original_edges[[i]][2])
        starting_directed_edge <- c(admix_name, original_directed_edges[[j]][2])
        flow_result <- flow(leaves, edges, directed_edges, starting_directed_edge, default_problem)
        if (flow_result$problem == FALSE) {
          graph <- root_graph(leaves, inner_nodes, flow_result$edges, flow_result$directed_edges, admixtures, outgroup)
          graph_list[[length(graph_list) + 1]] <- graph
        }
      }
    }
  }
  if (length(original_directed_edges) > 0) {
    for (i in seq(1, length(original_directed_edges))) {
      for (j in seq(1, length(original_edges))) {
        # From directed edge to edge:
        # Weeding out half of the duplicated cases.
        default_problem <- FALSE
        if (length(intersect(original_directed_edges[[i]], original_edges[[j]])) > 0) {
          common_node <- intersect(original_directed_edges[[i]], original_edges[[j]])
          for (k in seq(1, length(original_directed_edges))) {
            if (common_node == original_directed_edges[[k]][1] && k < i) {
              # The only undirected j, diverging at common node -> compare indices of i and k.
              default_problem <- TRUE
            }
            if (common_node == original_directed_edges[[k]][2]) {
              # The only undirected not i, meeting at common node -> always choose the directed j.
              # Actually this is a matter of whether the labels matter.
              if (labels_matter == FALSE) {
                default_problem <- TRUE
              }
            }
          }
        }
        # Optionally weeding out half of the cases that are indistinguishable from existing admixtures.
        if (labels_matter == FALSE) {
          for (k in seq(1, length(original_directed_edges))) {
            if (original_directed_edges[[k]][1] == original_directed_edges[[i]][1] &&
              length(intersect(original_directed_edges[[k]][2], original_edges[[j]])) > 0) {
              default_problem <- TRUE
            }
          }
        }
        # Direction [2]:
        admixtures <- original_admixtures
        for (k in seq(1, length(admixtures))) {
          admixture <- admixtures[[k]]
          if (admixture[1] == original_directed_edges[[i]][2]) {
            if (admixture[2] == original_directed_edges[[i]][1]) {
              admixture[2] <- inner_name
            }
            if (admixture[3] == original_directed_edges[[i]][1]) {
              admixture[3] <- inner_name
            }
            admixtures[[k]] <- admixture
          }
        }
        admixtures[[length(admixtures) + 1]] <- c(admix_name, inner_name, original_edges[[j]][1], admixture_variable_name)
        directed_edges <- original_directed_edges
        directed_edges[[i]] <- c(inner_name, original_directed_edges[[i]][2])
        directed_edges[[length(directed_edges) + 1]] <- c(inner_name, admix_name)
        directed_edges[[length(directed_edges) + 1]] <- c(original_edges[[j]][1], admix_name)
        edges <- original_edges
        edges[[j]] <- c(original_directed_edges[[i]][1], inner_name)
        edges[[length(edges) + 1]] <- c(admix_name, original_edges[[j]][2])
        starting_directed_edge <- c(admix_name, original_edges[[j]][2])
        flow_result <- flow(leaves, edges, directed_edges, starting_directed_edge, default_problem)
        if (flow_result$problem == FALSE) {
          graph <- root_graph(leaves, inner_nodes, flow_result$edges, flow_result$directed_edges, admixtures, outgroup)
          graph_list[[length(graph_list) + 1]] <- graph
        }
        # Direction [1]:
        admixtures <- original_admixtures
        for (k in seq(1, length(admixtures))) {
          admixture <- admixtures[[k]]
          if (admixture[1] == original_directed_edges[[i]][2]) {
            if (admixture[2] == original_directed_edges[[i]][1]) {
              admixture[2] <- inner_name
            }
            if (admixture[3] == original_directed_edges[[i]][1]) {
              admixture[3] <- inner_name
            }
            admixtures[[k]] <- admixture
          }
        }
        admixtures[[length(admixtures) + 1]] <- c(admix_name, inner_name, original_edges[[j]][2], admixture_variable_name)
        directed_edges <- original_directed_edges
        directed_edges[[i]] <- c(inner_name, original_directed_edges[[i]][2])
        directed_edges[[length(directed_edges) + 1]] <- c(inner_name, admix_name)
        directed_edges[[length(directed_edges) + 1]] <- c(original_edges[[j]][2], admix_name)
        edges <- original_edges
        edges[[j]] <- c(original_directed_edges[[i]][1], inner_name)
        edges[[length(edges) + 1]] <- c(admix_name, original_edges[[j]][1])
        starting_directed_edge <- c(admix_name, original_edges[[j]][1])
        flow_result <- flow(leaves, edges, directed_edges, starting_directed_edge, default_problem)
        if (flow_result$problem == FALSE) {
          graph <- root_graph(leaves, inner_nodes, flow_result$edges, flow_result$directed_edges, admixtures, outgroup)
          graph_list[[length(graph_list) + 1]] <- graph
        }
      }
      for (j in seq(1, length(original_directed_edges))) {
        if (i != j) {
          # From directed edge to directed edge:
          # Weeding out half of the duplicated cases.
          default_problem <- FALSE
          if (length(intersect(original_directed_edges[[i]], original_directed_edges[[j]])) > 0) {
            common_node <- intersect(original_directed_edges[[i]], original_directed_edges[[j]])
            for (k in seq(1, length(original_directed_edges))) {
              if (original_directed_edges[[k]][2] == common_node && k != i && k != j) {
                # All three directed -> choose i and j meeting at common node.
                # Actually this is a matter of whether the labels matter.
                if (labels_matter == FALSE) {
                  default_problem <- TRUE
                }
              }
            }
          }
          # Optionally weeding out half of the cases that are indistinguishable from existing admixtures.
          if (labels_matter == FALSE) {
            for (k in seq(1, length(original_directed_edges))) {
              if (original_directed_edges[[k]][1] == original_directed_edges[[i]][1] &&
                  original_directed_edges[[k]][2] == original_directed_edges[[j]][1]) {
                default_problem <- TRUE
              }
            }
            if (length(intersect(original_directed_edges[[i]], original_directed_edges[[j]])) == 0) {
              for (k in seq(1, j)) {
                if (original_directed_edges[[k]][1] == original_directed_edges[[i]][1] &&  k != i &&
                    original_directed_edges[[k]][2] == original_directed_edges[[j]][2] && k < j) {
                  for (l in seq(1, length(original_edges))) {
                    if (length(intersect(original_edges[[l]], original_directed_edges[[j]][1])) &&
                        length(intersect(original_edges[[l]], original_directed_edges[[k]][1]))) {
                      default_problem <- TRUE
                    }
                  }
                }
              }
            }
          }
          admixtures <- original_admixtures
          if (original_directed_edges[[i]][2] == original_directed_edges[[j]][2]) {
            for (k in seq(1,length(admixtures))) {
              admixture <- admixtures[[k]]
              if (admixture[1] == original_directed_edges[[i]][2]) {
                if (admixture[2] == original_directed_edges[[i]][1]) {
                  admixture[2] <- inner_name
                  admixture[3] <- admix_name
                } else {
                  admixture[3] <- inner_name
                  admixture[2] <- admix_name
                }
                admixtures[[k]] <- admixture
              }
            }
          } else {
            for (k in seq(1,length(admixtures))) {
              admixture <- admixtures[[k]]
              if (admixture[1] == original_directed_edges[[i]][2]) {
                if (admixture[2] == original_directed_edges[[i]][1]) {
                  admixture[2] <- inner_name
                }
                if (admixture[3] == original_directed_edges[[i]][1]) {
                  admixture[3] <- inner_name
                }
                admixtures[[k]] <- admixture
              }
              if (admixture[1] == original_directed_edges[[j]][2]) {
                if (admixture[2] == original_directed_edges[[j]][1]) {
                  admixture[2] <- admix_name
                }
                if (admixture[3] == original_directed_edges[[j]][1]) {
                  admixture[3] <- admix_name
                }
                admixtures[[k]] <- admixture
              }
            }
          }
          admixtures[[length(admixtures) + 1]] <- c(admix_name, inner_name, original_directed_edges[[j]][1], admixture_variable_name)
          directed_edges <- original_directed_edges
          directed_edges[[i]] <- c(inner_name, original_directed_edges[[i]][2])
          directed_edges[[j]] <- c(inner_name, admix_name)
          directed_edges[[length(directed_edges) + 1]] <- c(admix_name, original_directed_edges[[j]][2])
          directed_edges[[length(directed_edges) + 1]] <- c(original_directed_edges[[j]][1], admix_name)
          edges <- original_edges
          edges[[length(edges) + 1]] <- c(original_directed_edges[[i]][1], inner_name)
          starting_directed_edge <- c(admix_name, original_directed_edges[[j]][2])
          flow_result <- flow(leaves, edges, directed_edges, starting_directed_edge, default_problem)
          if (flow_result$problem == FALSE) {
            graph <- root_graph(leaves, inner_nodes, flow_result$edges, flow_result$directed_edges, admixtures, outgroup)
            graph_list[[length(graph_list) + 1]] <- graph
          }
        }
      }
    }
  }
  return(graph_list)
}

#' Adds a new admixture event to a graph.
#' 
#' Given an admixture graph, selects a child edge and two parent edges, disconnects the child edge
#' from its original parent node and connects it to the two parent edges with an admixture event,
#' if possible. Thus, contrary to \code{\link{add_an_admixture}}, the resulting graph need not be an
#' extension of the input graph in the sense that erasing one of the admixture edges we get
#' the original admixture graph. In practice, we know that when fitting data to admixture graphs,
#' the best graph with \eqn{k} admixture events is not always an extension like that from the best
#' graph with \eqn{k - 1} admixture events. Most likely it doesn't need to be an extension like
#' this (the two new admixture edges can both go where ever) either.
#' 
#' @param graph                    An admixture graph.
#' @param admixture_variable_name  A name for the new admixture proportion.
#' @param outgroup                 An optional parameter for the preferred outgroup.
#' 
#' @return A list of graphs made by adding a new admixture event to the input graph. The list contains
#'         duplicate elements, and may even contain graphs with \emph{eyes} (two inner nodes with the
#'         property that all the paths between any two leaves visits both or neither of them).
#' 
#' @seealso \code{\link{all_graphs}}
#' @seealso \code{\link{fit_permutations_and_graphs}}
#' @seealso \code{\link{fit_graph_list}}
#' @seealso \code{\link{add_a_leaf}}
#' @seealso \code{\link{add_an_admixture}}
#' @seealso \code{\link{make_an_outgroup}}
#' 
#' @export
add_an_admixture2 <- function(graph, admixture_variable_name, outgroup = "") {
  graph_list <- list()
  broken_graph <- break_graph(graph)
  # We might have to choose a different root after adding the admixture, so we start by removing
  # the original root and treating only the edges colliding in an admix event as directed and the
  # rest as undirected.
  leaves <- broken_graph$leaves
  original_inner_nodes <- broken_graph$inner_nodes
  original_edges <- broken_graph$edges
  original_admixtures <- broken_graph$admixtures
  root <- broken_graph$root
  for (i in seq(1, length(original_inner_nodes))) {
    node <- original_inner_nodes[i]
    if (node == root) {
      original_inner_nodes <- original_inner_nodes[-i]
      break
    }
  }
  memory <- ""
  to_be_deleted <- nchar(0)
  for (i in seq(1, length(original_edges))) {
    edge <- original_edges[[i]]
    if (edge[2] == root) {
      if (nchar(memory) == 0) {
        memory <- edge[1]
        to_be_deleted <- i
      } else {
        original_edges[[i]] <- c(memory, edge[1])
      }
    }
  }
  original_edges[[to_be_deleted]] <- NULL
  original_directed_edges <- list()
  for (admixture in original_admixtures) {
    original_directed_edges[[length(original_directed_edges) + 1]] <- c(admixture[2], admixture[1])
    original_directed_edges[[length(original_directed_edges) + 1]] <- c(admixture[3], admixture[1])
  }
  first_name <- paste("first", admixture_variable_name, sep = "_")
  second_name <- paste("second", admixture_variable_name, sep = "_")
  admix_name <- paste("admix", admixture_variable_name, sep = "_")
  # Lose the admixture events, let directed edges carry all the information.
  original_directed_edges <- save_admixture_information(original_admixtures, original_directed_edges)
  # Start:
  for (i in seq(1, length(original_edges))) {
    # Admixture event added in a normal edge:
    # Direction [2]:
    # Detach the selected edge from its old parent and merge some edges together.
    we_can_continue <- TRUE
    unnecessary_node <- original_edges[[i]][1]
    temporary_directed_edges <- original_directed_edges
    temporary_edges <- original_edges
    starting_directed_edge <- c(admix_name, temporary_edges[[i]][2])
    temporary_edges[[i]] <- NULL
    delete_indices <- numeric(0)
    endpoints <- character(0)
    for (l in seq(1, length(temporary_edges))) {
      if (temporary_edges[[l]][1] == unnecessary_node) {
        delete_indices <- c(delete_indices, l)
        endpoints <- c(endpoints, temporary_edges[[l]][2])
      }
      if (temporary_edges[[l]][2] == unnecessary_node) {
        delete_indices <- c(delete_indices, l)
        endpoints <- c(endpoints, temporary_edges[[l]][1])
      }
    }
    if (length(delete_indices) == 0) {
      we_can_continue <- FALSE
      temporary_inner_nodes <- original_inner_nodes
      temporary_inner_nodes <- c(temporary_inner_nodes, admix_name, second_name)
      temporary_edges[[length(temporary_edges) + 1]] <- temporary_edges[[1]]
      temporary_edges[[1]] <- starting_directed_edge
      for (k in seq(2, length(temporary_edges))) {
        # SPECIAL_N2-N
        edges <- temporary_edges
        directed_edges <- temporary_directed_edges
        inner_nodes <- temporary_inner_nodes
        edges[[length(edges) + 1]] <- c(temporary_edges[[k]][1], second_name)
        edges[[k]][1] <- second_name
        directed_edges[[length(directed_edges) + 1]] <- c(unnecessary_node, admix_name, admixture_variable_name, "")
        directed_edges[[length(directed_edges) + 1]] <- c(second_name, admix_name, "", admixture_variable_name)
        flow_result <- flow(leaves, edges, directed_edges, starting_directed_edge)
        admixtures <- load_admixture_information(directed_edges)
        for (l in seq(1, length(admixtures))) {
          if (admixtures[[l]][2] == admixtures[[l]][3]) {
            flow_result$problem <- TRUE
          }
        }
        if (flow_result$problem == FALSE) {
          graph <- root_graph(leaves, inner_nodes, flow_result$edges, flow_result$directed_edges, admixtures, outgroup)
          graph_list[[length(graph_list) + 1]] <- graph
        }
      }
      if (length(temporary_directed_edges) > 0) {
        for (k in seq(1, length(temporary_directed_edges))) {
          # SPECIAL_N2-D
          edges <- temporary_edges
          directed_edges <- temporary_directed_edges
          inner_nodes <- temporary_inner_nodes
          edges[[length(edges) + 1]] <- c(temporary_directed_edges[[k]][1], second_name)
          directed_edges[[k]][1] <- second_name
          directed_edges[[length(directed_edges) + 1]] <- c(unnecessary_node, admix_name, admixture_variable_name, "")
          directed_edges[[length(directed_edges) + 1]] <- c(second_name, admix_name, "", admixture_variable_name)
          flow_result <- flow(leaves, edges, directed_edges, starting_directed_edge)
          admixtures <- load_admixture_information(directed_edges)
          for (l in seq(1, length(admixtures))) {
            if (admixtures[[l]][2] == admixtures[[l]][3]) {
              flow_result$problem <- TRUE
            }
          }
            if (flow_result$problem == FALSE) {
            graph <- root_graph(leaves, inner_nodes, flow_result$edges, flow_result$directed_edges, admixtures, outgroup)
            graph_list[[length(graph_list) + 1]] <- graph
          }
        }
      }
    } else if (length(delete_indices) == 1) {
      for (l in seq(1, length(temporary_directed_edges))) {  
        if (temporary_directed_edges[[l]][1] == unnecessary_node) {
          temporary_edges[[delete_indices[1]]] <- NULL
          temporary_directed_edges[[l]][1] <- endpoints[1]
        }
      }
    } else if (length(delete_indices) == 2) {
      temporary_edges[[delete_indices[2]]] <- NULL
      temporary_edges[[delete_indices[1]]] <- NULL
      temporary_edges[[length(temporary_edges) + 1]] <- c(endpoints[1], endpoints[2])
    }
    temporary_inner_nodes <- original_inner_nodes[original_inner_nodes != unnecessary_node]
    temporary_inner_nodes <- c(temporary_inner_nodes, admix_name, first_name, second_name)
    temporary_edges[[length(temporary_edges) + 1]] <- temporary_edges[[1]]
    temporary_edges[[1]] <- starting_directed_edge
    # Attach the selected egde to two new parents.
    if (we_can_continue == TRUE) {
      for (j in seq(2, length(temporary_edges))) {
        for (k in seq(j, length(temporary_edges))) {
          # N2-NN
          edges <- temporary_edges
          directed_edges <- temporary_directed_edges
          inner_nodes <- temporary_inner_nodes
          edges[[j]][2] <- first_name
          edges[[length(edges) + 1]] <- c(first_name, temporary_edges[[j]][2])
          edges[[length(edges) + 1]] <- c(temporary_edges[[k]][1], second_name)
          edges[[k]][1] <- second_name
          # Looks inconsistent but this way it works even when j = k.
          directed_edges[[length(directed_edges) + 1]] <- c(first_name, admix_name, admixture_variable_name, "")
          directed_edges[[length(directed_edges) + 1]] <- c(second_name, admix_name, "", admixture_variable_name)
          flow_result <- flow(leaves, edges, directed_edges, starting_directed_edge)
          admixtures <- load_admixture_information(directed_edges)
          for (l in seq(1, length(admixtures))) {
            if (admixtures[[l]][2] == admixtures[[l]][3]) {
              flow_result$problem <- TRUE
            }
          }
          if (flow_result$problem == FALSE) {
            graph <- root_graph(leaves, inner_nodes, flow_result$edges, flow_result$directed_edges, admixtures, outgroup)
            graph_list[[length(graph_list) + 1]] <- graph
          }
        }
        if (length(temporary_directed_edges) > 0) {
          for (k in seq(1, length(temporary_directed_edges))) {
            # N2-ND
            edges <- temporary_edges
            directed_edges <- temporary_directed_edges
            inner_nodes <- temporary_inner_nodes
            edges[[j]] <- c(temporary_edges[[j]][1], first_name)
            edges[[length(edges) + 1]] <- c(first_name, temporary_edges[[j]][2])
            edges[[length(edges) + 1]] <- c(temporary_directed_edges[[k]][1], second_name)
            directed_edges[[k]][1] <- second_name
            directed_edges[[length(directed_edges) + 1]] <- c(first_name, admix_name, admixture_variable_name, "")
            directed_edges[[length(directed_edges) + 1]] <- c(second_name, admix_name, "", admixture_variable_name)
            flow_result <- flow(leaves, edges, directed_edges, starting_directed_edge)
            admixtures <- load_admixture_information(directed_edges)
            for (l in seq(1, length(admixtures))) {
              if (admixtures[[l]][2] == admixtures[[l]][3]) {
                flow_result$problem <- TRUE
              }
            }
            if (flow_result$problem == FALSE) {
              graph <- root_graph(leaves, inner_nodes, flow_result$edges, flow_result$directed_edges, admixtures, outgroup)
              graph_list[[length(graph_list) + 1]] <- graph
            }
          }
        }
      }
      if (length(temporary_directed_edges) > 0) {
        for (j in seq(1, length(temporary_directed_edges))) {
          for (k in seq(j, length(temporary_directed_edges))) {
            # N2-DD
            edges <- temporary_edges
            directed_edges <- temporary_directed_edges
            inner_nodes <- temporary_inner_nodes
            edges[[length(edges) + 1]] <- c(directed_edges[[j]][1], first_name)
            directed_edges[[j]][1] <- first_name
            edges[[length(edges) + 1]] <- c(directed_edges[[k]][1], second_name)
            directed_edges[[k]][1] <- second_name
            # Looks inconsistent but this way it works even when j = k.
            directed_edges[[length(directed_edges) + 1]] <- c(first_name, admix_name, admixture_variable_name, "")
            directed_edges[[length(directed_edges) + 1]] <- c(second_name, admix_name, "", admixture_variable_name)
            flow_result <- flow(leaves, edges, directed_edges, starting_directed_edge)
            admixtures <- load_admixture_information(directed_edges)
            for (l in seq(1, length(admixtures))) {
              if (admixtures[[l]][2] == admixtures[[l]][3]) {
                flow_result$problem <- TRUE
              }
            }
            if (flow_result$problem == FALSE) {
              graph <- root_graph(leaves, inner_nodes, flow_result$edges, flow_result$directed_edges, admixtures, outgroup)
              graph_list[[length(graph_list) + 1]] <- graph
            }
          }
        }
      }
    }
    # Direction [1]:
    # Detach the selected edge from its old parent and merge some edges together.
    we_can_continue <- TRUE
    unnecessary_node <- original_edges[[i]][2]
    temporary_directed_edges <- original_directed_edges
    temporary_edges <- original_edges
    starting_directed_edge <- c(admix_name, temporary_edges[[i]][1])
    temporary_edges[[i]] <- NULL
    delete_indices <- numeric(0)
    endpoints <- character(0)
    for (l in seq(1, length(temporary_edges))) {
      if (temporary_edges[[l]][1] == unnecessary_node) {
        delete_indices <- c(delete_indices, l)
        endpoints <- c(endpoints, temporary_edges[[l]][2])
      }
      if (temporary_edges[[l]][2] == unnecessary_node) {
        delete_indices <- c(delete_indices, l)
        endpoints <- c(endpoints, temporary_edges[[l]][1])
      }
    }
    if (length(delete_indices) == 0) {
      we_can_continue <- FALSE
      temporary_inner_nodes <- original_inner_nodes
      temporary_inner_nodes <- c(temporary_inner_nodes, admix_name, second_name)
      temporary_edges[[length(temporary_edges) + 1]] <- temporary_edges[[1]]
      temporary_edges[[1]] <- starting_directed_edge
      for (k in seq(2, length(temporary_edges))) {
        # SPECIAL_N1-N
        edges <- temporary_edges
        directed_edges <- temporary_directed_edges
        inner_nodes <- temporary_inner_nodes
        edges[[length(edges) + 1]] <- c(temporary_edges[[k]][1], second_name)
        edges[[k]][1] <- second_name
        directed_edges[[length(directed_edges) + 1]] <- c(unnecessary_node, admix_name, admixture_variable_name, "")
        directed_edges[[length(directed_edges) + 1]] <- c(second_name, admix_name, "", admixture_variable_name)
        flow_result <- flow(leaves, edges, directed_edges, starting_directed_edge)
        admixtures <- load_admixture_information(directed_edges)
        for (l in seq(1, length(admixtures))) {
          if (admixtures[[l]][2] == admixtures[[l]][3]) {
            flow_result$problem <- TRUE
          }
        }
        if (flow_result$problem == FALSE) {
          graph <- root_graph(leaves, inner_nodes, flow_result$edges, flow_result$directed_edges, admixtures, outgroup)
          graph_list[[length(graph_list) + 1]] <- graph
        }
      }
      if (length(temporary_directed_edges) > 0) {
        for (k in seq(1, length(temporary_directed_edges))) {
          # SPECIAL_N1-D
          edges <- temporary_edges
          directed_edges <- temporary_directed_edges
          inner_nodes <- temporary_inner_nodes
          edges[[length(edges) + 1]] <- c(temporary_directed_edges[[k]][1], second_name)
          directed_edges[[k]][1] <- second_name
          directed_edges[[length(directed_edges) + 1]] <- c(unnecessary_node, admix_name, admixture_variable_name, "")
          directed_edges[[length(directed_edges) + 1]] <- c(second_name, admix_name, "", admixture_variable_name)
          flow_result <- flow(leaves, edges, directed_edges, starting_directed_edge)
          admixtures <- load_admixture_information(directed_edges)
          for (l in seq(1, length(admixtures))) {
            if (admixtures[[l]][2] == admixtures[[l]][3]) {
              flow_result$problem <- TRUE
            }
          }
          if (flow_result$problem == FALSE) {
            graph <- root_graph(leaves, inner_nodes, flow_result$edges, flow_result$directed_edges, admixtures, outgroup)
            graph_list[[length(graph_list) + 1]] <- graph
          }
        }
      }
    } else if (length(delete_indices) == 1) {
      for (l in seq(1, length(temporary_directed_edges))) {  
        if (temporary_directed_edges[[l]][1] == unnecessary_node) {
          temporary_edges[[delete_indices[1]]] <- NULL
          temporary_directed_edges[[l]][1] <- endpoints[1]
        }
      }
    } else if (length(delete_indices) == 2) {
      temporary_edges[[delete_indices[2]]] <- NULL
      temporary_edges[[delete_indices[1]]] <- NULL
      temporary_edges[[length(temporary_edges) + 1]] <- c(endpoints[1], endpoints[2])
    }
    temporary_inner_nodes <- original_inner_nodes[original_inner_nodes != unnecessary_node]
    temporary_inner_nodes <- c(temporary_inner_nodes, admix_name, first_name, second_name)
    temporary_edges[[length(temporary_edges) + 1]] <- temporary_edges[[1]]
    temporary_edges[[1]] <- starting_directed_edge
    # Attach the selected egde to two new parents.
    if (we_can_continue == TRUE) {
      for (j in seq(2, length(temporary_edges))) {
        for (k in seq(j, length(temporary_edges))) {
          # N1-NN
          edges <- temporary_edges
          directed_edges <- temporary_directed_edges
          inner_nodes <- temporary_inner_nodes
          edges[[j]][2] <- first_name
          edges[[length(edges) + 1]] <- c(first_name, temporary_edges[[j]][2])
          edges[[length(edges) + 1]] <- c(temporary_edges[[k]][1], second_name)
          edges[[k]][1] <- second_name
          # Looks inconsistent but this way it works even when j = k.
          directed_edges[[length(directed_edges) + 1]] <- c(first_name, admix_name, admixture_variable_name, "")
          directed_edges[[length(directed_edges) + 1]] <- c(second_name, admix_name, "", admixture_variable_name)
          flow_result <- flow(leaves, edges, directed_edges, starting_directed_edge)
          admixtures <- load_admixture_information(directed_edges)
          for (l in seq(1, length(admixtures))) {
            if (admixtures[[l]][2] == admixtures[[l]][3]) {
              flow_result$problem <- TRUE
            }
          }
          if (flow_result$problem == FALSE) {
            graph <- root_graph(leaves, inner_nodes, flow_result$edges, flow_result$directed_edges, admixtures, outgroup)
            graph_list[[length(graph_list) + 1]] <- graph
          }
        }
        if (length(temporary_directed_edges) > 0) {
          for (k in seq(1, length(temporary_directed_edges))) {
            # N1-ND
            edges <- temporary_edges
            directed_edges <- temporary_directed_edges
            inner_nodes <- temporary_inner_nodes
            edges[[j]] <- c(temporary_edges[[j]][1], first_name)
            edges[[length(edges) + 1]] <- c(first_name, temporary_edges[[j]][2])
            edges[[length(edges) + 1]] <- c(temporary_directed_edges[[k]][1], second_name)
            directed_edges[[k]][1] <- second_name
            directed_edges[[length(directed_edges) + 1]] <- c(first_name, admix_name, admixture_variable_name,"")
            directed_edges[[length(directed_edges) + 1]] <- c(second_name, admix_name, "", admixture_variable_name)
            flow_result <- flow(leaves, edges, directed_edges, starting_directed_edge)
            admixtures <- load_admixture_information(directed_edges)
            for (l in seq(1, length(admixtures))) {
              if (admixtures[[l]][2] == admixtures[[l]][3]) {
                flow_result$problem <- TRUE
              }
            }
            if (flow_result$problem == FALSE) {
              graph <- root_graph(leaves, inner_nodes, flow_result$edges, flow_result$directed_edges, admixtures, outgroup)
              graph_list[[length(graph_list) + 1]] <- graph
            }
          }
        }
      }
      if (length(temporary_directed_edges) > 0) {
        for (j in seq(1, length(temporary_directed_edges))) {
          for (k in seq(j, length(temporary_directed_edges))) {
            # N1-DD
            edges <- temporary_edges
            directed_edges <- temporary_directed_edges
            inner_nodes <- temporary_inner_nodes
            edges[[length(edges) + 1]] <- c(directed_edges[[j]][1], first_name)
            directed_edges[[j]][1] <- first_name
            edges[[length(edges) + 1]] <- c(directed_edges[[k]][1], second_name)
            directed_edges[[k]][1] <- second_name
            # Looks inconsistent but this way it works even when j = k.
            directed_edges[[length(directed_edges) + 1]] <- c(first_name, admix_name, admixture_variable_name, "")
            directed_edges[[length(directed_edges) + 1]] <- c(second_name, admix_name, "", admixture_variable_name)
            flow_result <- flow(leaves, edges, directed_edges, starting_directed_edge)
            admixtures <- load_admixture_information(directed_edges)
            for (l in seq(1, length(admixtures))) {
              if (admixtures[[l]][2] == admixtures[[l]][3]) {
                flow_result$problem <- TRUE
              }
            }
            if (flow_result$problem == FALSE) {
              graph <- root_graph(leaves, inner_nodes, flow_result$edges, flow_result$directed_edges, admixtures, outgroup)
              graph_list[[length(graph_list) + 1]] <- graph
            }
          }
        }
      }
    }
  }
  if (length(original_directed_edges) > 0) {
    for (i in seq(1, length(original_directed_edges))) {
      # Admixture event added in a directed edge:
      # Detach the selected edge from its old parent and merge some edges together.
      we_can_continue <- TRUE
      unnecessary_node <- original_directed_edges[[i]][1]
      temporary_directed_edges <- original_directed_edges
      temporary_edges <- original_edges
      starting_directed_edge <- temporary_directed_edges[[i]]
      starting_directed_edge[1] <- admix_name
      temporary_directed_edges[[i]] <- NULL
      delete_indices <- numeric(0)
      endpoints <- character(0)
      for (l in seq(1, length(temporary_edges))) {
        if (temporary_edges[[l]][1] == unnecessary_node) {
          delete_indices <- c(delete_indices, l)
          endpoints <- c(endpoints, temporary_edges[[l]][2])
        }
        if (temporary_edges[[l]][2] == unnecessary_node) {
          delete_indices <- c(delete_indices, l)
          endpoints <- c(endpoints, temporary_edges[[l]][1])
        }
      }
      if (length(delete_indices) == 0) {
        we_can_continue <- FALSE
        temporary_inner_nodes <- original_inner_nodes
        temporary_inner_nodes <- c(temporary_inner_nodes, admix_name, second_name)
        temporary_directed_edges[[length(temporary_directed_edges) + 1]] <- temporary_directed_edges[[1]]
        temporary_directed_edges[[1]] <- starting_directed_edge
        starting_directed_edge <- starting_directed_edge[1:2]
        for (k in seq(1, length(temporary_edges))) {
          # SPECIAL_D-N
          edges <- temporary_edges
          directed_edges <- temporary_directed_edges
          inner_nodes <- temporary_inner_nodes
          edges[[length(edges) + 1]] <- c(temporary_edges[[k]][1], second_name)
          edges[[k]][1] <- second_name
          directed_edges[[length(directed_edges) + 1]] <- c(unnecessary_node, admix_name, admixture_variable_name, "")
          directed_edges[[length(directed_edges) + 1]] <- c(second_name, admix_name, "", admixture_variable_name)
          flow_result <- flow(leaves, edges, directed_edges, starting_directed_edge)
          admixtures <- load_admixture_information(directed_edges)
          for (l in seq(1, length(admixtures))) {
            if (admixtures[[l]][2] == admixtures[[l]][3]) {
              flow_result$problem <- TRUE
            }
          }
          if (flow_result$problem == FALSE) {
            graph <- root_graph(leaves, inner_nodes, flow_result$edges, flow_result$directed_edges, admixtures, outgroup)
            graph_list[[length(graph_list) + 1]] <- graph
          }
        }
        for (k in seq(2, length(temporary_directed_edges))) {
          # SPECIAL_D-D
          edges <- temporary_edges
          directed_edges <- temporary_directed_edges
          inner_nodes <- temporary_inner_nodes
          edges[[length(edges) + 1]] <- c(temporary_directed_edges[[k]][1], second_name)
          directed_edges[[k]][1] <- second_name
          directed_edges[[length(directed_edges) + 1]] <- c(unnecessary_node, admix_name, admixture_variable_name, "")
          directed_edges[[length(directed_edges) + 1]] <- c(second_name, admix_name, "", admixture_variable_name)
          flow_result <- flow(leaves, edges, directed_edges, starting_directed_edge)
          admixtures <- load_admixture_information(directed_edges)
          for (l in seq(1, length(admixtures))) {
            if (admixtures[[l]][2] == admixtures[[l]][3]) {
              flow_result$problem <- TRUE
            }
          }
          if (flow_result$problem == FALSE) {
            graph <- root_graph(leaves, inner_nodes, flow_result$edges, flow_result$directed_edges, admixtures, outgroup)
            graph_list[[length(graph_list) + 1]] <- graph
          }
        }
      } else if (length(delete_indices) == 1) {
        for (l in seq(1, length(temporary_directed_edges))) {  
          if (temporary_directed_edges[[l]][1] == unnecessary_node) {
            temporary_edges[[delete_indices[1]]] <- NULL
            temporary_directed_edges[[l]][1] <- endpoints[1]
          }
        }
      } else if (length(delete_indices) == 2) {
        temporary_edges[[delete_indices[2]]] <- NULL
        temporary_edges[[delete_indices[1]]] <- NULL
        temporary_edges[[length(temporary_edges) + 1]] <- c(endpoints[1], endpoints[2])
      }
      temporary_inner_nodes <- original_inner_nodes[original_inner_nodes != unnecessary_node]
      temporary_inner_nodes <- c(temporary_inner_nodes, admix_name, first_name, second_name)
      temporary_directed_edges[[length(temporary_directed_edges) + 1]] <- temporary_directed_edges[[1]]
      temporary_directed_edges[[1]] <- starting_directed_edge
      starting_directed_edge <- starting_directed_edge[1:2]
      # Attach the selected egde to two new parents.
      if (we_can_continue == TRUE) {
        for (j in seq(1, length(temporary_edges))) {
          for (k in seq(j, length(temporary_edges))) {
            # D-NN
            edges <- temporary_edges
            directed_edges <- temporary_directed_edges
            inner_nodes <- temporary_inner_nodes
            edges[[j]][2] <- first_name
            edges[[length(edges) + 1]] <- c(first_name, temporary_edges[[j]][2])
            edges[[length(edges) + 1]] <- c(temporary_edges[[k]][1], second_name)
            edges[[k]][1] <- second_name
            # Looks inconsistent but this way it works even when j = k.
            directed_edges[[length(directed_edges) + 1]] <- c(first_name, admix_name, admixture_variable_name, "")
            directed_edges[[length(directed_edges) + 1]] <- c(second_name, admix_name, "", admixture_variable_name)
            flow_result <- flow(leaves, edges, directed_edges, starting_directed_edge)
            admixtures <- load_admixture_information(directed_edges)
            for (l in seq(1, length(admixtures))) {
              if (admixtures[[l]][2] == admixtures[[l]][3]) {
                flow_result$problem <- TRUE
              }
            }
            if (flow_result$problem == FALSE) {
              graph <- root_graph(leaves, inner_nodes, flow_result$edges, flow_result$directed_edges, admixtures, outgroup)
              graph_list[[length(graph_list) + 1]] <- graph
            }
          }
          for (k in seq(2, length(temporary_directed_edges))) {
            # D-ND
            edges <- temporary_edges
            directed_edges <- temporary_directed_edges
            inner_nodes <- temporary_inner_nodes
            edges[[j]] <- c(temporary_edges[[j]][1], first_name)
            edges[[length(edges) + 1]] <- c(first_name, temporary_edges[[j]][2])
            edges[[length(edges) + 1]] <- c(temporary_directed_edges[[k]][1], second_name)
            directed_edges[[k]][1] <- second_name
            directed_edges[[length(directed_edges) + 1]] <- c(first_name, admix_name, admixture_variable_name, "")
            directed_edges[[length(directed_edges) + 1]] <- c(second_name, admix_name, "", admixture_variable_name)
            flow_result <- flow(leaves, edges, directed_edges, starting_directed_edge)
            admixtures <- load_admixture_information(directed_edges)
            for (l in seq(1, length(admixtures))) {
              if (admixtures[[l]][2] == admixtures[[l]][3]) {
                flow_result$problem <- TRUE
              }
            }
            if (flow_result$problem == FALSE) {
              graph <- root_graph(leaves, inner_nodes, flow_result$edges, flow_result$directed_edges, admixtures, outgroup)
              graph_list[[length(graph_list) + 1]] <- graph
            }
          }
        }
        for (j in seq(2, length(temporary_directed_edges))) {
          for (k in seq(j, length(temporary_directed_edges))) {
            # D-DD
            edges <- temporary_edges
            directed_edges <- temporary_directed_edges
            inner_nodes <- temporary_inner_nodes
            edges[[length(edges) + 1]] <- c(directed_edges[[j]][1], first_name)
            directed_edges[[j]][1] <- first_name
            edges[[length(edges) + 1]] <- c(directed_edges[[k]][1], second_name)
            directed_edges[[k]][1] <- second_name
            # Looks inconsistent but this way it works even when j = k.
            directed_edges[[length(directed_edges) + 1]] <- c(first_name, admix_name, admixture_variable_name, "")
            directed_edges[[length(directed_edges) + 1]] <- c(second_name, admix_name, "", admixture_variable_name)
            flow_result <- flow(leaves, edges, directed_edges, starting_directed_edge)
            admixtures <- load_admixture_information(directed_edges)
            for (l in seq(1, length(admixtures))) {
              if (admixtures[[l]][2] == admixtures[[l]][3]) {
                flow_result$problem <- TRUE
              }
            }
            if (flow_result$problem == FALSE) {
              graph <- root_graph(leaves, inner_nodes, flow_result$edges, flow_result$directed_edges, admixtures, outgroup)
              graph_list[[length(graph_list) + 1]] <- graph
            }
          }
        }
      }
    }
  }
  return(graph_list)
}

#' Make an outgroup.
#' 
#' Given a graph and a leaf, tries to put the root of the graph on the edge leading to the leaf.
#' If not possible (\emph{i. e.} if the leaf has admixture in its ancestry), puts the root
#' somewhere else.
#' 
#' @param graph        An admixture graph.
#' @param outgroup     A leaf we want to be the outgroup.
#' @param all_neutral  For when other functions need to root graphs in a neutral way.
#' 
#' @return An admixture graph with the given leaf as an outgroup, if possible.
#'
#' @seealso \code{\link{make_permutations}}
#' @seealso \code{\link{four_leaves_graphs}}
#' @seealso \code{\link{five_leaves_graphs}}
#' @seealso \code{\link{six_leaves_graphs}}
#' @seealso \code{\link{seven_leaves_graphs}}
#' @seealso \code{\link{eight_leaves_trees}}
#' @seealso \code{\link{fit_permutations_and_graphs}}
#' @seealso \code{\link{fit_graph_list}}
#' @seealso \code{\link{add_a_leaf}}
#' @seealso \code{\link{add_an_admixture}}
#' @seealso \code{\link{add_an_admixture2}}
#' 
#' @examples
#' \donttest{
#' # Here is a little family tree of some dinosaur-like animals.
#' 
#' species <- c("triceratops", "crocodile", "diplodocus", "tyrannosaurus", "chicken")
#' graph <- five_leaves_graphs[[1]](species)
#' plot(graph)
#' 
#' # Of course we know that while this is correct as an undirected graph, "crocodile"
#' # should really be the outgroup.
#' 
#' graph <- make_an_outgroup(graph, "crocodile")
#' plot(graph)
#' 
#' # Strictly speaking the graph is still a little misleading because unfortunately
#' # the (non-bird) dinosaurs are extinct :-(
#' }
#' 
#' @export
make_an_outgroup <- function(graph, outgroup = "", all_neutral = FALSE) {
  broken_graph <- break_graph(graph)
  leaves <- broken_graph$leaves
  inner_nodes <- broken_graph$inner_nodes
  edges <- broken_graph$edges
  admixtures <- broken_graph$admixtures
  root <- broken_graph$root
  for (i in seq(1, length(inner_nodes))) {
    node <- inner_nodes[i]
    if (node == root) {
      inner_nodes <- inner_nodes[-i]
      break
    }
  }
  memory <- ""
  to_be_deleted <- nchar(0)
  for (i in seq(1, length(edges))) {
    edge <- edges[[i]]
    if (edge[2] == root) {
      if (nchar(memory) == 0) {
        memory <- edge[1]
        to_be_deleted <- i
      } else {
        edges[[i]] <- c(memory, edge[1])
      }
    }
  }
  edges[[to_be_deleted]] <- NULL
  directed_edges <- list()
  for (admixture in admixtures) {
    directed_edges[[length(directed_edges) + 1]] <- c(admixture[2], admixture[1])
    directed_edges[[length(directed_edges) + 1]] <- c(admixture[3], admixture[1])
  }
  graph <- root_graph(leaves, inner_nodes, edges, directed_edges, admixtures, outgroup, all_neutral)
  return(graph)
}

# Breaking up a graph back to its components.
break_graph <- function(graph) {
  nodes <- graph$nodes
  edges <- list()
  admixtures <- list()
  parents <- graph$parents
  probs <- graph$probs
  for (i in seq(1, NROW(parents))) {
    match <- which(parents[i, ] == TRUE)
    if (length(match) == 0) {
      root <- nodes[i]
    } else if (length(match) == 1) {
      edges[[length(edges) + 1]] <- c(nodes[i], nodes[match[1]])
    } else if (length(match) == 2) {
      if (nchar(probs[i, match[1]]) > nchar(probs[i, match[2]])) {
        admixtures[[length(admixtures) + 1]] <- c(nodes[i], nodes[match[2]], nodes[match[1]], probs[i, match[2]])
      } else {
        admixtures[[length(admixtures) + 1]] <- c(nodes[i], nodes[match[1]], nodes[match[2]], probs[i, match[1]])
      }
    }
  }
  return(list(leaves = graph$leaves, inner_nodes = graph$inner_nodes, edges = edges, admixtures = admixtures, root = root))
}

# Given material for a graph known to been all right, choose a suitable root and build the graph.
root_graph <- function(leaves, inner_nodes, edges, directed_edges, admixtures, outgroup = "", all_neutral = FALSE) {
  for (admixture in admixtures) {
    flow_result <- flow(leaves, edges, directed_edges, c(admixture[2], admixture[1]))
    edges <- flow_result$edges
    directed_edges <- flow_result$directed_edges
  }
  inner_nodes <- c(inner_nodes, "R")
  roots <- list()
  graph_list <- list()
  if (nchar(outgroup) == 0) {
    if (all_neutral == FALSE) {
      root <- list(index = 1, edge = edges[[1]])
      roots[[1]] <- root
    } else {
      # This is a bit ugly but we need a way to choose the root in a neutral way.
      for (j in seq(1, length(edges))) {
        roots[[j]] <- list(index = j, edge = edges[[j]])
      }
    }
  } else {
    root_index <- 1
    for (j in seq(1, length(edges))) {
      edge <- edges[[j]]
      if (edge[1] == outgroup || edge[2] == outgroup) {
        root_index <- j
      }
    }
    roots[[1]] <- list(index = root_index, edge = edges[[root_index]])
  }
  original_edges <- edges
  original_directed_edges <- directed_edges
  original_admixtures <- admixtures
  for(j in seq(1, length(roots))) {
    new_edges <- original_edges
    new_directed_edges <- original_directed_edges
    new_admixtures <- original_admixtures
    edge <- roots[[j]]$edge
    new_edges[[roots[[j]]$index]] <- NULL
    new_directed_edges[[length(new_directed_edges) + 1]] <- c("R", edge[1])
    new_directed_edges[[length(new_directed_edges) + 1]] <- c("R", edge[2])
    flow_result <- flow(leaves, new_edges, new_directed_edges, c("R", edge[1]))
    new_edges <- flow_result$edges
    new_directed_edges <- flow_result$directed_edges
    flow_result <- flow(leaves, new_edges, new_directed_edges, c("R", edge[2]))
    # There should be no undirected edges anymore.
    new_directed_edges <- flow_result$directed_edges
    # This is a little embarassing but the edge directions are in fact wrong at the moment.
    edge_argument <- character(0)
    admix_argument <- character(0)
    for (j in seq(1, length(new_directed_edges))) {
      edge_argument <- c(edge_argument, edge(new_directed_edges[[j]][2], new_directed_edges[[j]][1]))
    }
    if (length(new_admixtures) > 0) {
      for (k in seq(1, length(new_admixtures))) {
        edge_argument <- c(edge_argument, admixture_edge(new_admixtures[[k]][1], new_admixtures[[k]][2],
                                                         new_admixtures[[k]][3]))
        admix_argument <- c(admix_argument, admix_props(new_admixtures[[k]][1], new_admixtures[[k]][2],
                                                        new_admixtures[[k]][3], new_admixtures[[k]][4]))
      }
    }
    new_edges <- parent_edges(edge_argument)
    if (length(new_admixtures) > 0) {
      new_admixtures <- admixture_proportions(admix_argument)
    } else {
      new_admixtures <- NULL
    }
    new_graph <- agraph(leaves, inner_nodes, new_edges, new_admixtures)
    graph_list[[length(graph_list) + 1]] <- new_graph
  }
  if (all_neutral == FALSE) {graph_list <- graph_list[[1]]}
  return(graph_list)
}

# Detecting directed loops and collisions, gives a direction to undirected edges.
# The variable default_problem allows the output problem to be forced to TRUE.
flow <- function(leaves, edges, directed_edges, starting_directed_edge, default_problem = FALSE) {
  forgive <- TRUE
  problem <- FALSE
  if (default_problem == TRUE) {
    problem <- TRUE
  }
  active_directed_edges <- list(starting_directed_edge)
  # We might need to remove the starting directed edge from the normal edge list and add it to directed edge list:
  erase <- 0
  if (length(edges) > 0) {
    for (i in seq(1, length(edges))) {
      edge <- edges[[i]]
      if (edge[1] == starting_directed_edge[1] && edge[2] == starting_directed_edge[2]) {
        erase <- i
        directed_edges[[length(directed_edges) + 1]] <- starting_directed_edge
      }
      if (edge[1] == starting_directed_edge[2] && edge[2] == starting_directed_edge[1]) {
        erase <- i
        directed_edges[[length(directed_edges) + 1]] <- starting_directed_edge
      }
    }
  }
  if (erase > 0) {
    edges[[erase]] <- NULL
  }
  # Now start the iteration:
  while (problem == FALSE && length(active_directed_edges) > 0) {
    active_directed_edge <- active_directed_edges[[1]]
    # Detect directed loops:
    if (active_directed_edge[1] == starting_directed_edge[1] && active_directed_edge[2] == starting_directed_edge[2]) {
      if (forgive == TRUE) {
        forgive <- FALSE
      } else {
        problem <- TRUE
      }
    }
    continues <- FALSE
    # Flow to directed edges:
    if (length(directed_edges) > 0) {
      for (i in seq(1, length(directed_edges))) {
        selected_directed_edge <- directed_edges[[i]]
        if (active_directed_edge[2] == selected_directed_edge[1]) {
          active_directed_edges[[length(active_directed_edges) + 1]] <- selected_directed_edge
          continues <- TRUE
        }
      }
    }
    # Flow to undirected edges and direct them:
    if (length(edges) > 0) {
      to_be_erased <- integer(0)
      for (i in seq(1, length(edges))) {
        selected_edge <- edges[[i]]
        if (active_directed_edge[2] == selected_edge[1]) {
          active_directed_edges[[length(active_directed_edges) + 1]] <- selected_edge
          directed_edges[[length(directed_edges) + 1]] <- selected_edge
          to_be_erased <- c(to_be_erased, i)
          continues <- TRUE
        }
        if (active_directed_edge[2] == selected_edge[2]) {
          active_directed_edges[[length(active_directed_edges) + 1]] <- c(selected_edge[2], selected_edge[1])
          directed_edges[[length(directed_edges) + 1]] <- c(selected_edge[2], selected_edge[1])
          to_be_erased <- c(to_be_erased, i)
          continues <- TRUE
        }
      }
      while (length(to_be_erased) > 0) {
        edges[[to_be_erased[length(to_be_erased)]]] <- NULL
        to_be_erased <- to_be_erased[-length(to_be_erased)]
      }
    }
    # If the flow stops we must be at a leaf:
    if (continues == FALSE) {
      if (active_directed_edge[2] %in% leaves) {
      } else {
        problem <- TRUE  
      }
    }
    active_directed_edges[[1]] <- NULL
  }
  return(list(problem = problem, edges = edges, directed_edges = directed_edges))
}

# It is vastly more convenient to forget the admixture events for now and just let the directed edges carry
# all the information.
save_admixture_information <- function(admixtures, directed_edges) {
  if (length(admixtures) > 0) {
    for (i in seq(1, length(admixtures))) {
      for (j in seq(1, length(directed_edges))) {
        if (directed_edges[[j]][1] == admixtures[[i]][2] && directed_edges[[j]][2] == admixtures[[i]][1]) {
          directed_edges[[j]][3] <- admixtures[[i]][4]
          directed_edges[[j]][4] <- ""
        }
        if (directed_edges[[j]][1] == admixtures[[i]][3] && directed_edges[[j]][2] == admixtures[[i]][1]) {
          directed_edges[[j]][3] <- ""
          directed_edges[[j]][4] <- admixtures[[i]][4]
        }
      }
    }
  }
  return(directed_edges)
}

# Then we can rebuild the admixture events from the directed edges.
load_admixture_information <- function(directed_edges) {
  admixtures <- list()
  for (i in seq(1, length(directed_edges))) {
    found_a_match <- FALSE
    if (length(admixtures) > 0) {
      for (j in seq(1, length(admixtures))) {
        if (directed_edges[[i]][2] == admixtures[[j]][1]) {
          found_a_match <- TRUE
          if (directed_edges[[i]][3] == "") {
            admixtures[[j]][3] <- directed_edges[[i]][1]
          } else {
            admixtures[[j]][2] <- directed_edges[[i]][1]
          }
        }
      }
    }
    if (found_a_match == FALSE) {
      if (directed_edges[[i]][3] == "") {
        admixtures[[length(admixtures) + 1]] <- c(directed_edges[[i]][2], "", directed_edges[[i]][1], directed_edges[[i]][4])
      } else {
        admixtures[[length(admixtures) + 1]] <- c(directed_edges[[i]][2], directed_edges[[i]][1], "", directed_edges[[i]][3])
      }
    }
  }
  return(admixtures)
}

# Recursive renaming function that unfortunately produces also duplicate names.
invent_name <- function(graph, node) {
  if (node %in% graph$leaves) {
    newname <- paste(match(node, sort(graph$leaves)))
  } else {
    children <- character(0)
    for (j in seq(1, NCOL(graph$children))) {
      if (graph$children[node, j] == TRUE) {
        children <- c(children, invent_name(graph, colnames(graph$children)[j]))
      }
    }
    children <- sort(children)
    newname <- "("
    for (j in seq(1, length(children))) {
      if (j != 1) {newname <- paste(newname, ",", sep = "")}
      newname <- paste(newname, children[j], sep = "")
    }
    newname <- paste(newname, ")", sep = "")
  }
  return(newname)
}

#' Rename nodes.
#'
#' Changes the names of the nodes of a graph.
#' Capable of giving new standard names to the inner nodes in a way that only depends on the graph topology
#' (without the root) and the leaf names. This is necessary when detecting when graphs are identical up to inner
#' node and admixture proportion names, see \code{\link{canonise_graph}} and \code{\link{remove_duplicates}}.
#' 
#' @param graph     The graph to be renamed.
#' @param newnames  A list of new names, given in the form \code{list(old = "new")}. Nodes not listed will keep
#'                  their old name, unless no list is given at all, in which case the leaves keep their old
#'                  names while the inner nodes get new standardised names.
#'                  
#' @return A graph with new node names.
#' 
#' @export
rename_nodes <- function(graph, newnames = list()) {
  if (length(newnames) == 0) {
    # Here we make a default naming system of the inner nodes (as a function of leaf names and graph topology).
    # That is for the graph comparison and maybe isomorphism functions.
    # First, put the outgroup in some neutral place.
    graph_list <- make_an_outgroup(graph, all_neutral = TRUE)
    root_names <- character(0)
    for (j in seq(1, length(graph_list))) {
      root_names <- c(root_names, invent_name(graph_list[[j]], "R"))
    }
    graph <- graph_list[[match(sort(root_names)[1], root_names)]]
    # Now that the graph shape is neutral, rename the inner nodes and rebuild it.
    for (node in graph$inner_nodes) {
      newnames[[node]] <- invent_name(graph, node)
    }
    # But there is a problem. With some graphs the inner node names are not unique now. So we need to name them
    # in a way that depends not only on their descendants, and does not depend on the original inner node order.
    to_do_list <- newnames
    while (length(to_do_list) > 0) {
      inspect <- names(to_do_list)[1] # The original name of the node we study. Replicants must be dealt with.
      multi <- list()
      for (j in seq(1, length(to_do_list))) {
        if (newnames[[names(to_do_list)[j]]] == newnames[[inspect]]) {
          # The current name of the node is the same as the current name of "inspect".
          # Divide into bins according to "family name" depending (current names of) parents.
          row <- graph$parents[names(to_do_list)[j], ]
          family_name <- ""
          if (length(row[row == TRUE]) == 1) {
            family_name <- paste("*", newnames[[colnames(graph$parents)[which(row == TRUE)[1]]]], sep = "")
          }
          if (length(row[row == TRUE]) == 2) {
            parent_names <- sort(c(newnames[[colnames(graph$parents)[which(row == TRUE)[1]]]],
                                   newnames[[colnames(graph$parents)[which(row == TRUE)[2]]]]))
            family_name <- paste(parent_names[1], parent_names[2], sep = "")
          }
          if (family_name %in% names(multi)) {
            multi[[family_name]] <- c(multi[[family_name]], j)
          } else {
            multi[[family_name]] <- c(j)
          }
        }
      }
      # The situation is not actually that complicated: multi has either:
      # 1) A single index (remove).
      # 2) Two indices from a single parent (rename how ever and continue).
      # 3) Two indices from similarly named different parents (wait for the parents to get new names).
      # 4) Two indices from differently named parents (rename alphabetically and remove).
      away <- numeric(0) # Keep track of indices that we can remove.
      end <- numeric(0) # Keep track of indices that we must investigate again later.
      if (length(multi) == 1) {
        if (length(multi[[1]]) == 1) {
          # Case 1).
          away <- c(away, multi[[1]])
        } else {
          if (substr(names(multi)[1], 1, 1) == "*") {
            # Case 2).
            newnames[[names(to_do_list)[multi[[1]][2]]]] <- paste(newnames[[names(to_do_list)[multi[[1]][2]]]],
                                                                  "+", sep = "")
            away <- c(away, multi[[1]])
          } else {
            # Case 3).
            end <- c(end, multi[[1]])
          }
        }
      } else {
        # Case 4).
        if (names(multi)[1] > names(multi)[2]) {
          newnames[[names(to_do_list)[multi[[1]][1]]]] <- paste(newnames[[names(to_do_list)[multi[[1]][1]]]],
                                                                "+", sep = "")
        } else {
          newnames[[names(to_do_list)[multi[[2]][1]]]] <- paste(newnames[[names(to_do_list)[multi[[2]][1]]]],
                                                                "+", sep = "")
        }
        away <- c(away, multi[[1]], multi[[2]])
      }
      # Then remove the indices corresponding to dealt-with nodes.
      if (length(away) > 0) {
        for (j in seq(1, length(away))) {
          to_do_list[[away[length(away) + 1 - j]]] <- NULL
        }
      }
      # And finally move the unfinished nodes to the end of the line.
      if (length(end) > 0) {
        for (j in seq(1, length(end))) {
          to_do_list[[length(to_do_list) + 1]] <- to_do_list[[end[length(end) + 1 - j]]]
          names(to_do_list)[length(to_do_list)] <- names(to_do_list)[end[length(end) + 1 - j]]
          to_do_list[[end[length(end) + 1 - j]]] <- NULL
        }
      }
    }
    # Now that every name is unique in a way that depends only on the leaf names and the graph topology, let's
    # simplify them for clarity (50 letter words don't look nice as column names of a matrix).
    # Suppose no one names their leaves as numbers inside parenthesis, if they do the plotting function at least 
    # will crash.
    word_order <- sort(unlist(newnames))
    for (j in seq(1, length(newnames))) {
      newnames[[j]] <- paste("(", match(newnames[[j]], word_order), ")", sep = "")
    }
  }
  leaf_amount <- length(graph$leaves)
  nodes <- list()
  for (node in graph$nodes) {
    nodes[[node]] <- node
  }
  nodes <- utils::modifyList(nodes, newnames)
  nodes <- unlist(nodes)
  leaves <- sort(nodes[1:leaf_amount])
  inner_nodes <- sort(nodes[leaf_amount + 1:(length(nodes) - leaf_amount)])
  nodes <- c(leaves, inner_nodes)
  edge_vector <- character(0)
  for (j in seq(1, NROW(graph$parents))) {
    row <- graph$parents[j, ]
    if (length(row[row == TRUE]) == 1) {
      edge_vector <- c(edge_vector, edge(nodes[rownames(graph$parents)[j]],
                                         nodes[rownames(graph$parents)[which(row == TRUE)[1]]]))
    }
    if (length(row[row == TRUE]) == 2) {
      if (nchar(graph$probs[j, which(row == TRUE)[1]]) < nchar(graph$probs[j, which(row == TRUE)[2]])) {
        k <- which(row == TRUE)[1]
        l <- which(row == TRUE)[2]
        p <- graph$probs[j, which(row == TRUE)[1]]
      } else {
        k <- which(row == TRUE)[2]
        l <- which(row == TRUE)[1]
        p <- graph$probs[j, which(row == TRUE)[2]]
      }
      edge_vector <- c(edge_vector, admixture_edge(nodes[rownames(graph$parents)[j]],
                                                   nodes[rownames(graph$parents)[k]],
                                                   nodes[rownames(graph$parents)[l]],
                                                   p))
    }
  }
  return(agraph(leaves, inner_nodes, parent_edges(edge_vector)))
}

#' Canonise graph.
#' 
#' Given a graph builds a unique logical vector depending only on the leaf names and the graph
#' topology (not the inner node names, root position or the input order of edges or inner nodes).
#' Can be used to detect graph isomorphism, that is, to weed out duplicates from a graph list.
#' Only for comparison of graphs with the same leaf set!
#' 
#' @param graph  An admixture graph.
#' 
#' @return A logical vector coding the parental matrix of a canonised version of the input graph.
#' 
#' @export
canonise_graph <- function(graph) {
  V <- logical(0)
  matrix <- rename_nodes(graph)$parents
  for (i in seq(1, NCOL(matrix))) {
    for (j in seq(1, NROW(matrix))) {
      V <- c(V, matrix[i, j])
    }
  }
  return(V)
}

#' Remove duplicate graphs from a list.
#' 
#' Using \code{\link{canonise_graph}} to calculate unique characteristic logical vector for each
#' graph in a given list of graphs, then sorts the list according to this attribute and remove
#' repeated graphs.
#' Leaf names count so that graphs with permuted leaves are considered different.
#' Also organises similar graphs next to each other if \code{organise} is \code{TRUE}, but this
#' is extremely slow.
#' 
#' @param graph_list    A list of graphs with the same leaf set.
#' @param organise      If \code{TRUE} also organises isomorphic graphs (now disregarding leaf names)
#'                      next to each other.
#' @param return_piles  If \code{TRUE} and \code{organise} is also \code{TRUE}, the output will be a
#'                      list of lists of isomorphic graphs instead of one big list.
#' 
#' @return A list of graphs that are all different (or a list of lists of graphs if both \code{organise}
#'         and \code{return_piles} are \code{TRUE}).
#' 
#' @export
remove_duplicates <- function(graph_list, organise = FALSE, return_piles = FALSE) {
  # Sorting:
  M <- as.numeric(canonise_graph(graph_list[[1]]))
  M <- rbind(M, as.numeric(canonise_graph(graph_list[[1]]))) # R is stupid.
  for (graph in graph_list) {
    word <- as.numeric(canonise_graph(graph))
    while (length(word) > NCOL(M)) {
      M <- cbind(M, rep(3, NROW(M)))
    }
    if (length(word) < NCOL(M)) {
      word <- c(word, rep(3, NCOL(M) - length(word)))
    }
    M <- rbind(M, word)
  }
  M <- M[-1, ]
  M <- M[-1, ]
  rownames(M) <- 1:NROW(M)
  N <- as.data.frame(M)
  N <- N[do.call(order, as.list(N)), ] 
  # Marking the indices of the duplicates in a vector:
  remove <- numeric(0)
  if (NROW(M) > 1) {
    for (j in seq(2, NROW(M))) {
      if (all(N[j, ] == N[j - 1, ])) {
        remove <- c(remove, as.numeric(rownames(N)[j]))
      }
    }
  }
  remove <- -sort(-remove)
  # Removing the duplicates:
  if (length(remove) > 0) {
    for (j in remove) {
      graph_list[[j]] <- NULL
    }
  }
  if (organise == TRUE) {
    leaves <- graph_list[[1]]$leaves
    permutations <- make_permutations(leaves)
    P <- 1
    piles <- list()
    piles[[1]] <- list()
    piles[[1]][[1]] <- graph_list[[1]]
    if (length(graph_list) > 1) {
      for (j in seq(2, length(graph_list))) {
        graph <- graph_list[[j]]
        home_found <- FALSE
        for (permutation in permutations) {
          newnames <- list()
          for (k in seq(1, length(leaves))) {
            newnames[[leaves[k]]] <- permutation[k]
          }
          permuted_graph <- rename_nodes(graph, newnames)
          for (p in seq(1, P)) {
            if (length(canonise_graph(piles[[p]][[1]])) == length(canonise_graph(permuted_graph))
                && home_found == FALSE) {
              if (all(canonise_graph(piles[[p]][[1]]) == canonise_graph(permuted_graph))) {
                piles[[p]][[length(piles[[p]]) + 1]] <- graph
                home_found <- TRUE
              }
            }
          }
        }
        if (home_found == FALSE) {
          P <- P + 1
          piles[[P]] <- list()
          piles[[P]][[1]] <- graph
        }
      }
    }
    new_list <- list()
    for (p in seq(1, P)) {
      for (j in seq(1, length(piles[[p]]))) {
        new_list[[length(new_list) + 1]] <- piles[[p]][[j]]
      }
    }
    if (return_piles == TRUE) {
      graph_list <- piles
    } else {
      graph_list <- new_list
    }
  }
  return(graph_list)
}

#' All graphs.
#' 
#' Gives a list of all the graphs with at most a given number of admixture events.
#' No duplicates.
#' 
#' @param populations       A vector of populations (leaf names).
#' @param admixture_events  The maximum number of admixture events allowed.
#' 
#' @return A list of admixture graphs.
#' 
#' @seealso \code{\link{fit_graph_list}}
#' 
#' @export
all_graphs <- function(populations, admixture_events) {
  # Base structure.
  leaves <- populations[1:2]
  inner_nodes <- c("R")
  edges <- parent_edges(c(edge(leaves[1], "R"),
                          edge(leaves[2], "R")))
  admixtures <- NULL
  base_structure <- agraph(leaves, inner_nodes, edges, admixtures)
  # Build a list of trees (no duplicates).
  tree_list <- list(base_structure)
  if (length(populations > 2)) {
    for (j in seq(3, length(populations))) {
      temp_list <- list()
      for (tree in tree_list) {
        temp_list <- c(temp_list, add_a_leaf(tree, populations[j]))
      }
      tree_list <- temp_list
    }
  }
  # Add admixture events while weeding out the duplicates as soon as they appear.
  graph_list <- tree_list
  if (admixture_events > 0) {
    old_list <- tree_list
    for (j in seq(1, admixture_events)) {
      new_list <- list()
      canon_list <- list()
      name <- paste("p", j, sep = "")
      for (graph in old_list) {
        temp_list <- add_an_admixture(graph, name)
        # add_an_admixture2() sometimes makes eyes but this doesn't.
        # We're sorting them here in order to not create one giant list to drain all the memory.
        for (candidate in temp_list) {
          canon <- paste(as.character(as.numeric(canonise_graph(candidate))), sep = "", collapse = "")
          if (length(new_list) > 0) {
            old_length <- length(canon_list)
            canon_list <- try_to_add(canon_list, canon)
            new_length <- length(canon_list)
            if (new_length > old_length) {
              new_list[[length(new_list) + 1]] <- candidate
            }
          } else {
            canon_list <- list(canon)
            new_list <- list(candidate)
          }
        }
      }
      old_list <- new_list
      graph_list <- c(graph_list, new_list)
    }
  }
  return(graph_list)
}

try_to_add <- function(canon_list, canon) {
  B <- length(canon_list) # Current block size.
  N <- 1 + floor(B/2) # Middle element of the current block.
  while (B > 0) {
    if (canon == canon_list[[N]]) {
      B <- 0
    } else if (canon < canon_list[[N]]) {
      if (B < 2) {
        canon_list <- append(canon_list, canon, N - 1)
      }
      B <- ceiling((B - 1)/2)
      N <- N - ceiling(B/2)
    } else if (canon > canon_list[[N]]) {
      if (B < 3) {
        canon_list <- append(canon_list, canon, N)
      }
      B <- floor((B - 1)/2)
      N <- N + 1 + floor(B/2)
    }
  }
  return(canon_list)
}

#' Graph to vector.
#' 
#' Encodes an \code{\link{agraph}} object into a logical vector for saving memory. The admixture
#' proportion names will be lost.
#' 
#' @param graph  A graph.
#' 
#' @return The logical vector representing the graph.
#' 
#' @export
graph_to_vector <- function(graph) {
  matrix <- graph$parents
  vector <- as.vector(matrix)
  names(vector) <- colnames(matrix)
  return(vector)
}

#' Vector to graph.
#' 
#' Interprets a logical vector back to an \code{\link{agraph}} object. The admixture
#' proportion names are now lost.
#' 
#' @param vector  A logical vector.
#' 
#' @return The graph corresponding to the vector.
#' 
#' @export
vector_to_graph <- function(vector) {
  column_amount <- 1
  while(is.na(names(vector)[column_amount]) == FALSE) {
    column_amount <- column_amount + 1
  }
  column_amount <- column_amount - 1
  parents <- matrix(vector, ncol = column_amount)
  colnames(parents) <- names(vector)[1:column_amount]
  rownames(parents) <- names(vector)[1:column_amount]
  leaves <- character(0)
  inner_nodes <- character(0)
  edges <- rep("", 3)
  admixtures <- rbind(rep("", 3), rep("", 3))
  admix_name_count <- 1
  for (j in seq(1, column_amount)) {
    C <- parents[, j]
    R <- parents[j, ]
    if (length(C[C == TRUE]) == 0) {
      leaves <- c(leaves, rownames(parents)[j])
    } else {
      inner_nodes <- c(inner_nodes, rownames(parents)[j])
    }
    if (length(R[R == TRUE]) == 1) {
      edges <- rbind(edges, parent_edges(edge(rownames(parents)[j], rownames(parents)[which(R == T)])))
    } else if (length(R[R == TRUE]) == 2) {
      admix_name <- paste("p", admix_name_count, sep = "")
      admix_name_count <- admix_name_count + 1
      edges <- rbind(edges, parent_edges(admixture_edge(rownames(parents)[j],
                                                        rownames(parents)[which(R == T)[1]],
                                                        rownames(parents)[which(R == T)[2]])))
      admixtures <- rbind(admixtures, admixture_proportions(c(admix_props(rownames(parents)[j],
                                                                          rownames(parents)[which(R == T)[1]], 
                                                                          rownames(parents)[which(R == T)[2]], 
                                                                          admix_name))))
    }
  }
  edges <- edges[-1, ]
  if (NROW(admixtures) == 2) {
    admixtures <- NULL
  } else {
    admixtures <- admixtures[-1, ]
    admixtures <- admixtures[-1, ]
  }
  graph <- agraph(leaves, inner_nodes, edges, admixtures)
  return(graph)
}

#' Seven leaves trees.
#' 
#' The function \code{\link{seven_leaves_graphs}} is better than this as it also contains
#' graphs with one admixture. (This function is only kept for legacy reasons.)
#' Even more obsolete since the introduction of \code{\link{all_graphs}}.
#' 
#' @format A list of functions on seven leaves.
#'         The outputs of these functions are \code{\link{agraph}} objects.
#'   
#' @export
seven_leaves_trees <- list(
  tree_1 = function(leaves) {
    inner_nodes <- c("R", "x", "y", "z", "u", "v")
    edges <- parent_edges(c(
      edge("x", "R"),
      edge("y", "x"),
      edge("z", "y"),
      edge("u", "z"),
      edge("v", "u"),
      edge(leaves[1], "R"),
      edge(leaves[2], "x"),
      edge(leaves[3], "y"),
      edge(leaves[4], "z"),
      edge(leaves[5], "u"),
      edge(leaves[6], "v"),
      edge(leaves[7], "v")
    ))
    admixtures <- NULL
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  tree_2 = function(leaves) {
    inner_nodes <- c("R", "x", "y", "z", "u", "v")
    edges <- parent_edges(c(
      edge("x", "y"),
      edge("y", "R"),
      edge("z", "u"),
      edge("u", "R"),
      edge("v", "u"),
      edge(leaves[1], "x"),
      edge(leaves[2], "x"),
      edge(leaves[3], "y"),
      edge(leaves[4], "z"),
      edge(leaves[5], "z"),
      edge(leaves[6], "v"),
      edge(leaves[7], "v")
    ))
    admixtures <- NULL
    agraph(leaves, inner_nodes, edges, admixtures)
  }
)