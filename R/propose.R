
#' Make a list of all permutations of populations.
#' 
#' Right now it is hardwired to create permutations for four or five populations because
#' I just copied code from the two fit functions, but this should be generalized if we
#' need that at some point.
#'
#' @export
make_permutations <- function(populations) {
  P <- list()
  if (length(populations) == 4) {
    for (i in seq(1, 4)) {
      permutation <- rep("", 4)
      permutation[1] <- populations[i]
      for (j in seq(1, 3)) {
        permutation[2] <- populations[-i][j]
        for (k in seq(1, 2)) {
          permutation[3] <- populations[-i][-j][k]
          permutation[4] <- populations[-i][-j][-k][1]
          P[[length(P) + 1]] <- permutation
        }
      }
    }  
  } else if (length(populations) == 5) {
    for (i in seq(1, 5)) {
      permutation <- rep("", 5)
      permutation[1] <- populations[i]
      for (j in seq(1, 4)) {
        permutation[2] <- populations[-i][j]
        for (k in seq(1, 3)) {
          permutation[3] <- populations[-i][-j][k]
          for (l in seq(1, 2)) {
            permutation[4] <- populations[-i][-j][-k][l]
            permutation[5] <- populations[-i][-j][-k][-l][1]
            P[[length(P) + 1]] <- permutation
          }
        }
      }
    }    
  } else if (length(populations) == 6) {
    for (i in seq(1, 6)) {
      permutation <- rep("", 6)
      permutation[1] <- populations[i]
      for (j in seq(1, 5)) {
        permutation[2] <- populations[-i][j]
        for (k in seq(1, 4)) {
          permutation[3] <- populations[-i][-j][k]
          for (l in seq(1, 3)) {
            permutation[4] <- populations[-i][-j][-k][l]
            for (m in seq(1, 2)) {
              permutation[5] <- populations[-i][-j][-k][-l][m]  
              permutation[6] <- populations[-i][-j][-k][-l][-m][1]  
              P[[length(P) + 1]] <- permutation
            }
          }
        }
      }
    }
  } else if (length(populations) == 7) {
    for (i in seq(1, 7)) {
      permutation <- rep("", 7)
      permutation[1] <- populations[i]
      for (j in seq(1, 6)) {
        permutation[2] <- populations[-i][j]
        for (k in seq(1, 5)) {
          permutation[3] <- populations[-i][-j][k]
          for (l in seq(1, 4)) {
            permutation[4] <- populations[-i][-j][-k][l]
            for (m in seq(1, 3)) {
              permutation[5] <- populations[-i][-j][-k][-l][m]
              for (n in seq(1, 2)) {
                permutation[6] <- populations[-i][-j][-k][-l][-m][n]
                permutation[7] <- populations[-i][-j][-k][-l][-m][-n][1]  
                P[[length(P) + 1]] <- permutation
              }
            }
          }
        }
      }
    }
  } else if (length(populations) == 8) {
    for (i in seq(1, 8)) {
      permutation <- rep("", 8)
      permutation[1] <- populations[i]
      for (j in seq(1, 7)) {
        permutation[2] <- populations[-i][j]
        for (k in seq(1, 6)) {
          permutation[3] <- populations[-i][-j][k]
          for (l in seq(1, 5)) {
            permutation[4] <- populations[-i][-j][-k][l]
            for (m in seq(1, 4)) {
              permutation[5] <- populations[-i][-j][-k][-l][m]
              for (n in seq(1, 3)) {
                permutation[6] <- populations[-i][-j][-k][-l][-m][n]
                for (o in seq(1, 2)) {
                  permutation[7] <- populations[-i][-j][-k][-l][-m][-n][o]
                  permutation[8] <- populations[-i][-j][-k][-l][-m][-n][-o][1]  
                  P[[length(P) + 1]] <- permutation
                }
              }
            }
          }
        }
      }
    }
  }
  return(P)
}

#' @export
four_leaves_graphs <- list(
  # single tree
  tree = function(leaves) {
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
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  # one admixture event
  one_admixture_1 = function(leaves) {
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
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  one_admixture_2 = function(leaves) {
    inner_nodes <- c("R", "y", "z", "w", "M")
    edges <- parent_edges(c(
      edge("y", "R"),
      edge("z", "R"),
      edge("w", "z"),
      edge(leaves[1], "R"), 
      edge(leaves[2], "M"), 
      edge(leaves[3], "w"),
      edge(leaves[4], "w"),
      admixture_edge("M", "y", "z")
    ))
    admixtures <- admixture_proportions(c(
      admix_props("M", "y", "z", "a")
    ))
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  one_admixture_3 = function(leaves) {
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
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  # two admixture events
  two_admixtures_1 = function(leaves) {
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
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  two_admixtures_2 = function(leaves) {
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
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  two_admixtures_3 = function(leaves) {
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
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  two_admixtures_4 = function(leaves) {
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
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  two_admixtures_5 = function(leaves) {
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
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  two_admixtures_6 = function(leaves) {
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
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  two_admixtures_7 = function(leaves) {
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
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  two_admixtures_8 = function(leaves) {
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
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  two_admixtures_9 = function(leaves) {
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
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  two_admixtures_10 = function(leaves) {
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
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  two_admixtures_11 = function(leaves) {
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
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  two_admixtures_12 = function(leaves) {
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
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  two_admixtures_13 = function(leaves) {
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
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  two_admixtures_14 = function(leaves) {
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
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  two_admixtures_15 = function(leaves) {
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
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  two_admixtures_16 = function(leaves) {
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
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  two_admixtures_17 = function(leaves) {
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
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  two_admixtures_18 = function(leaves) {
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
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  two_admixtures_19 = function(leaves) {
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
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  two_admixtures_20 = function(leaves) {
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
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  two_admixtures_21 = function(leaves) {
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
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  two_admixtures_22 = function(leaves) {
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
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  two_admixtures_23 = function(leaves) {
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
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  two_admixtures_24 = function(leaves) {
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
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  two_admixtures_25 = function(leaves) {
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
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  two_admixtures_26 = function(leaves) {
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
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  two_admixtures_27 = function(leaves) {
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
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  two_admixtures_28 = function(leaves) {
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
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  two_admixtures_29 = function(leaves) {
    inner_nodes <- c("R", "x", "y", "z", "w", "M", "N")
    edges <- parent_edges(c(
      edge("x", "R"),
      edge("y", "R"),
      edge("z", "y"),
      edge("w", "z"),
      edge(leaves[1], "x"), 
      edge(leaves[2], "N"), 
      edge(leaves[3], "w"),
      edge(leaves[4], "w"),
      admixture_edge("M", "x", "y"),
      admixture_edge("N", "M", "z")
    ))
    admixtures <- admixture_proportions(c(
      admix_props("M", "x", "y", "a"),
      admix_props("N", "M", "z", "b")
    ))
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  two_admixtures_30 = function(leaves) {
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
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  two_admixtures_31 = function(leaves) {
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
    agraph(leaves, inner_nodes, edges, admixtures)
  }
)

#' @export
five_leaves_graphs <- list(
  # tree
  tree = function(leaves) {
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
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  # one admixture
  one_admixture_1 = function(leaves) {
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
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  one_admixture_2 = function(leaves) {
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
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  one_admixture_3 = function(leaves) {
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
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  one_admixture_4 = function(leaves) {
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
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  one_admixture_5 = function(leaves) {
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
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  one_admixture_6 = function(leaves) {
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
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  one_admixture_7 = function(leaves) {
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
    agraph(leaves, inner_nodes, edges, admixtures)
  }
)

#' @export
six_leaves_graphs <- list(
  # tree
  tree_1 = function(leaves) {
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
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  tree_2 = function(leaves) {
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
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  # one admixture
  one_admixture_1 = function(leaves) {
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
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  one_admixture_2 = function(leaves) {
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
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  one_admixture_3 = function(leaves) {
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
    agraph(leaves, inner_nodes, edges, admixtures)
  },

  one_admixture_4 = function(leaves) {
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
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  one_admixture_5 = function(leaves) {
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
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  one_admixture_6 = function(leaves) {
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
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  one_admixture_7 = function(leaves) {
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
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  one_admixture_8 = function(leaves) {
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
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  one_admixture_9 = function(leaves) {
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
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  one_admixture_10 = function(leaves) {
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
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  one_admixture_11 = function(leaves) {
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
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  one_admixture_12 = function(leaves) {
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
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  one_admixture_13 = function(leaves) {
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
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  one_admixture_14 = function(leaves) {
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
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  one_admixture_15 = function(leaves) {
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
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  one_admixture_16 = function(leaves) {
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
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  one_admixture_17 = function(leaves) {
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
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  one_admixture_18 = function(leaves) {
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
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  one_admixture_19 = function(leaves) {
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
    agraph(leaves, inner_nodes, edges, admixtures)
  }
)

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

#' @export
eight_leaves_trees <- list(
  tree_1 = function(leaves) {
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
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  tree_2 = function(leaves) {
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
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  tree_3 = function(leaves) {
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
    agraph(leaves, inner_nodes, edges, admixtures)
  }
)

#' Combine a list of permutations with a list of (parameterized) graphs to fit them
#' 
#' @param permutations List of population permutations
#' @param graphs       List of functions for producing graphs
#' 
#' @return a list of fitted graphs.
#'
#' @export
fit_permutations_and_graphs <- function(data, permutations, graphs) {
  library(foreach)
  library(doParallel)
  cl <- makeCluster(2)
  registerDoParallel(cl)
  foreach(i = seq(1, length(permutations)), .combine = c, .packages = "admixturegraph") %:%
  foreach(j = seq(1, length(graphs)), .packages = "admixturegraph") %dopar% { 
    permutation <- permutations[[i]]
    graph_function <- graphs[[j]]
    graph <- graph_function(permutation)
    result <- fast_fit(filter_on_leaves(data, graph), graph)
    return(result)
  }
}

#' @export
add_a_leaf <- function(graph, leaf_name) {
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
    for (k in seq(1, length(new_admixtures))) {
      edge_argument <- c(edge_argument, admixture_edge(new_admixtures[[k]][1], new_admixtures[[k]][2],
                                                       new_admixtures[[k]][3]))
      admix_argument <- c(admix_argument, admix_props(new_admixtures[[k]][1], new_admixtures[[k]][2],
                                                      new_admixtures[[k]][3], new_admixtures[[k]][4]))
    }
    edges <- parent_edges(edge_argument)
    admixtures <- admixture_proportions(admix_argument)
    skip <- FALSE
    if (chosen[2] == root) {
      if (forgive == TRUE) {
        forgive <- FALSE
      } else {
        skip <- TRUE
      }
    }
    if (skip == FALSE) {
      graph_list[[length(graph_list) + 1]] <- agraph(leaves, inner_nodes, edges, admixtures)
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
      graph_list[[length(graph_list) + 1]] <- agraph(leaves, inner_nodes, edges, admixtures)
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
      graph_list[[length(graph_list) + 1]] <- agraph(leaves, inner_nodes, edges, admixtures)
    }  
  }
  return(graph_list)
}

#' @export
add_an_admixture <- function(graph, admixture_variable_name, labels_matter = FALSE) {
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
        # (I'm planning to write a complete description of how this works as a pdf.)
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
          graph <- root_graph(leaves, inner_nodes, flow_result$edges, flow_result$directed_edges, admixtures)
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
          graph <- root_graph(leaves, inner_nodes, flow_result$edges, flow_result$directed_edges, admixtures)
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
          graph <- root_graph(leaves, inner_nodes, flow_result$edges, flow_result$directed_edges, admixtures)
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
              default_problem <- TRUE
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
          graph <- root_graph(leaves, inner_nodes, flow_result$edges, flow_result$directed_edges, admixtures)
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
          graph <- root_graph(leaves, inner_nodes, flow_result$edges, flow_result$directed_edges, admixtures)
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
                default_problem <- TRUE
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
            graph <- root_graph(leaves, inner_nodes, flow_result$edges, flow_result$directed_edges, admixtures)
            graph_list[[length(graph_list) + 1]] <- graph
          }
        }
      }
    }
  }
  return(graph_list)
}

#' @export
add_an_admixture2 <- function(graph, admixture_variable_name) {
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
          graph <- root_graph(leaves, inner_nodes, flow_result$edges, flow_result$directed_edges, admixtures)
          graph_list[[length(graph_list) + 1]] <- graph
        }
      }
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
          graph <- root_graph(leaves, inner_nodes, flow_result$edges, flow_result$directed_edges, admixtures)
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
            graph <- root_graph(leaves, inner_nodes, flow_result$edges, flow_result$directed_edges, admixtures)
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
              graph <- root_graph(leaves, inner_nodes, flow_result$edges, flow_result$directed_edges, admixtures)
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
              graph <- root_graph(leaves, inner_nodes, flow_result$edges, flow_result$directed_edges, admixtures)
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
          graph <- root_graph(leaves, inner_nodes, flow_result$edges, flow_result$directed_edges, admixtures)
          graph_list[[length(graph_list) + 1]] <- graph
        }
      }
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
          graph <- root_graph(leaves, inner_nodes, flow_result$edges, flow_result$directed_edges, admixtures)
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
            graph <- root_graph(leaves, inner_nodes, flow_result$edges, flow_result$directed_edges, admixtures)
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
              graph <- root_graph(leaves, inner_nodes, flow_result$edges, flow_result$directed_edges, admixtures)
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
              graph <- root_graph(leaves, inner_nodes, flow_result$edges, flow_result$directed_edges, admixtures)
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
            graph <- root_graph(leaves, inner_nodes, flow_result$edges, flow_result$directed_edges, admixtures)
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
            graph <- root_graph(leaves, inner_nodes, flow_result$edges, flow_result$directed_edges, admixtures)
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
              graph <- root_graph(leaves, inner_nodes, flow_result$edges, flow_result$directed_edges, admixtures)
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
              graph <- root_graph(leaves, inner_nodes, flow_result$edges, flow_result$directed_edges, admixtures)
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
              graph <- root_graph(leaves, inner_nodes, flow_result$edges, flow_result$directed_edges, admixtures)
              graph_list[[length(graph_list) + 1]] <- graph
            }
          }
        }
      }
    }
  }
  return(graph_list)
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
root_graph <- function(leaves, inner_nodes, edges, directed_edges, admixtures) {
  for (admixture in admixtures) {
    flow_result <- flow(leaves, edges, directed_edges, c(admixture[2], admixture[1]))
    edges <- flow_result$edges
    directed_edges <- flow_result$directed_edges
  }
  inner_nodes <- c(inner_nodes, "R")
  edge <- edges[[1]]
  edges[[1]] <- NULL
  directed_edges[[length(directed_edges) + 1]] <- c("R", edge[1])
  directed_edges[[length(directed_edges) + 1]] <- c("R", edge[2])
  flow_result <- flow(leaves, edges, directed_edges, c("R", edge[1]))
  edges <- flow_result$edges
  directed_edges <- flow_result$directed_edges
  flow_result <- flow(leaves, edges, directed_edges, c("R", edge[2]))
  # There should be no undirected edges anymore.
  directed_edges <- flow_result$directed_edges
  # This is a little embarassing but the edge directions are in fact wrong at the moment.
  edge_argument <- character(0)
  admix_argument <- character(0)
  for (j in seq(1, length(directed_edges))) {
    edge_argument <- c(edge_argument, edge(directed_edges[[j]][2], directed_edges[[j]][1]))
  }
  for (k in seq(1, length(admixtures))) {
    edge_argument <- c(edge_argument, admixture_edge(admixtures[[k]][1], admixtures[[k]][2],
                                                     admixtures[[k]][3]))
    admix_argument <- c(admix_argument, admix_props(admixtures[[k]][1], admixtures[[k]][2],
                                                    admixtures[[k]][3], admixtures[[k]][4]))
  }
  edges <- parent_edges(edge_argument)
  admixtures <- admixture_proportions(admix_argument)
  return(agraph(leaves, inner_nodes, edges, admixtures))
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
