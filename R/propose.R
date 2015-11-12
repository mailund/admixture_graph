
#' Make a list of all permutations of populations.
#' 
#' Right now it is hardwired to create permutations for four or five populations because
#' I just copied code from the two fit functions, but this should be generalized if we
#' need that at some point.
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
  }
  return(P)
}



four_leaves_graphs <- list(
  # single tree
  tree = function(A,B,C,D) {
    leaves <- c(A, B, C, D)
    inner_nodes <- c("R", "x", "y")
    edges <- parent_edges(c(
      edge("x", "R"),
      edge("y", "x"),
      edge(A, "R"),
      edge(B, "x"),
      edge(C, "y"),
      edge(D, "y")
    ))
    admixtures <- NULL
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  # one admixture event
  one_admixture_1 = function(A,B,C,D) {
    leaves <- c(A, B, C, D)
    inner_nodes <- c("R", "x", "y", "z", "M")
    edges <- parent_edges(c(
      edge("x", "R"),
      edge("y", "x"),
      edge("z", "x"),
      edge(A, "R"),
      edge(B, "y"), 
      edge(C, "M"),
      edge(D, "z"),
      admixture_edge("M", "y", "z")
    ))
    admixtures <- admixture_proportions(c(
      admix_props("M", "y", "z", "a")
    ))
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  one_admixture_2 = function(A, B, C, D) {
    leaves <- c(A, B, C, D)
    inner_nodes <- c("R", "y", "z", "w", "M")
    edges <- parent_edges(c(
      edge("y", "R"),
      edge("z", "R"),
      edge("w", "z"),
      edge(A, "R"), 
      edge(B, "M"), 
      edge(C, "w"),
      edge(D, "w"),
      admixture_edge("M", "y", "z")
    ))
    admixtures <- admixture_proportions(c(
      admix_props("M", "y", "z", "a")
    ))
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  one_admixture_3 = function(A, B, C, D) {
    leaves <- c(A, B, C, D)
    inner_nodes <- c("R", "y", "z", "w", "M")
    edges <- parent_edges(c(
      edge("y", "R"),
      edge("z", "R"),
      edge("w", "M"),
      edge(A, "y"), 
      edge(B, "w"), 
      edge(C, "w"),
      edge(D, "z"),
      admixture_edge("M", "y", "z")
    ))
    admixtures <- admixture_proportions(c(
      admix_props("M", "y", "z", "a")
    ))
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  # two admixture events
  two_admixtures_1 = function(A, B, C, D) {
    leaves <- c(A, B, C, D)
    inner_nodes <- c("R", "x", "y", "z", "w", "M", "N")
    edges <- parent_edges(c(
      edge("x", "R"),
      edge("y", "R"),
      edge("z", "y"),
      edge("w", "z"),
      edge(A, "x"), 
      edge(B, "M"), 
      edge(C, "N"),
      edge(D, "w"),
      admixture_edge("M", "x", "y"),
      admixture_edge("N", "z", "w")
    ))
    admixtures <- admixture_proportions(c(
      admix_props("M", "x", "y", "a"),
      admix_props("N", "z", "w", "b")
    ))
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  two_admixtures_2 = function(A, B, C, D) {
    leaves <- c(A, B, C, D)
    inner_nodes <- c("R", "x", "y", "z", "w", "M", "N")
    edges <- parent_edges(c(
      edge("x", "R"),
      edge("y", "R"),
      edge("z", "M"),
      edge("w", "z"),
      edge(A, "x"), 
      edge(B, "w"), 
      edge(C, "N"),
      edge(D, "y"),
      admixture_edge("M", "x", "y"),
      admixture_edge("N", "z", "w")
    ))
    admixtures <- admixture_proportions(c(
      admix_props("M", "x", "y", "a"),
      admix_props("N", "z", "w", "b")
    ))
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  two_admixtures_3 = function(A, B, C, D) {
    leaves <- c(A, B, C, D)
    inner_nodes <- c("R", "x", "y", "z", "w", "M", "N")
    edges <- parent_edges(c(
      edge("x", "R"),
      edge("y", "x"),
      edge("z", "x"),
      edge("w", "M"),
      edge(A, "R"), 
      edge(B, "N"), 
      edge(C, "w"),
      edge(D, "z"),
      admixture_edge("M", "y", "z"),
      admixture_edge("N", "y", "w")
    ))
    admixtures <- admixture_proportions(c(
      admix_props("M", "y", "z", "a"),
      admix_props("N", "y", "w", "b")
    ))
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  two_admixtures_4 = function(A, B, C, D) {
    leaves <- c(A, B, C, D)
    inner_nodes <- c("R", "x", "y", "z", "w", "M", "N")
    edges <- parent_edges(c(
      edge("x", "R"),
      edge("y", "R"),
      edge("z", "y"),
      edge("w", "z"),
      edge(A, "x"), 
      edge(B, "M"), 
      edge(C, "w"),
      edge(D, "N"),
      admixture_edge("M", "x", "z"),
      admixture_edge("N", "w", "y")
    ))
    admixtures <- admixture_proportions(c(
      admix_props("M", "x", "z", "a"),
      admix_props("N", "w", "y", "b")
    ))
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  two_admixtures_5 = function(A, B, C, D) {
    leaves <- c(A, B, C, D)
    inner_nodes <- c("R", "x", "y", "z", "w", "M", "N")
    edges <- parent_edges(c(
      edge("x", "R"),
      edge("y", "x"),
      edge("z", "x"),
      edge("w", "z"),
      edge(A, "R"), 
      edge(B, "y"), 
      edge(C, "N"),
      edge(D, "w"),
      admixture_edge("M", "y", "z"),
      admixture_edge("N", "M", "w")
    ))
    admixtures <- admixture_proportions(c(
      admix_props("M", "y", "z", "a"),
      admix_props("N", "M", "w", "b")
    ))
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  two_admixtures_6 = function(A, B, C, D) {
    leaves <- c(A, B, C, D)
    inner_nodes <- c("R", "x", "y", "z", "w", "M", "N")
    edges <- parent_edges(c(
      edge("x", "R"),
      edge("y", "R"),
      edge("z", "y"),
      edge("w", "y"),
      edge(A, "x"), 
      edge(B, "M"), 
      edge(C, "N"),
      edge(D, "w"),
      admixture_edge("M", "x", "z"),
      admixture_edge("N", "z", "w")
    ))
    admixtures <- admixture_proportions(c(
      admix_props("M", "x", "z", "a"),
      admix_props("N", "z", "w", "b")
    ))
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  two_admixtures_7 = function(A, B, C, D) {
    leaves <- c(A, B, C, D)
    inner_nodes <- c("R", "x", "y", "z", "w", "M", "N")
    edges <- parent_edges(c(
      edge("x", "R"),
      edge("y", "R"),
      edge("z", "x"),
      edge("w", "z"),
      edge(A, "x"), 
      edge(B, "M"), 
      edge(C, "N"),
      edge(D, "y"),
      admixture_edge("M", "z", "w"),
      admixture_edge("N", "w", "y")
    ))
    admixtures <- admixture_proportions(c(
      admix_props("M", "z", "w", "a"),
      admix_props("N", "w", "y", "b")
    ))
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  two_admixtures_8 = function(A, B, C, D) {
    leaves <- c(A, B, C, D)
    inner_nodes <- c("R", "x", "y", "z", "w", "M", "N")
    edges <- parent_edges(c(
      edge("x", "R"),
      edge("y", "R"),
      edge("z", "y"),
      edge("w", "z"),
      edge(A, "x"), 
      edge(B, "M"), 
      edge(C, "w"),
      edge(D, "N"),
      admixture_edge("M", "x", "w"),
      admixture_edge("N", "z", "y")
    ))
    admixtures <- admixture_proportions(c(
      admix_props("M", "x", "w", "a"),
      admix_props("N", "z", "y", "b")
    ))
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  two_admixtures_9 = function(A, B, C, D) {
    leaves <- c(A, B, C, D)
    inner_nodes <- c("R", "x", "y", "z", "w", "M", "N")
    edges <- parent_edges(c(
      edge("x", "R"),
      edge("y", "x"),
      edge("z", "x"),
      edge("w", "z"),
      edge(A, "R"), 
      edge(B, "y"), 
      edge(C, "N"),
      edge(D, "z"),
      admixture_edge("M", "y", "w"),
      admixture_edge("N", "M", "w")
    ))
    admixtures <- admixture_proportions(c(
      admix_props("M", "y", "w", "a"),
      admix_props("N", "M", "w", "b")
    ))
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  two_admixtures_10 = function(A, B, C, D) {
    leaves <- c(A, B, C, D)
    inner_nodes <- c("R", "x", "y", "z", "w", "M", "N")
    edges <- parent_edges(c(
      edge("x", "R"),
      edge("y", "R"),
      edge("z", "y"),
      edge("w", "z"),
      edge(A, "x"), 
      edge(B, "N"), 
      edge(C, "w"),
      edge(D, "z"),
      admixture_edge("M", "x", "y"),
      admixture_edge("N", "M", "w")
    ))
    admixtures <- admixture_proportions(c(
      admix_props("M", "x", "y", "a"),
      admix_props("N", "M", "w", "b")
    ))
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  two_admixtures_11 = function(A, B, C, D) {
    leaves <- c(A, B, C, D)
    inner_nodes <- c("R", "x", "y", "z", "w", "M", "N")
    edges <- parent_edges(c(
      edge("x", "R"),
      edge("y", "R"),
      edge("z", "M"),
      edge("w", "y"),
      edge(A, "x"), 
      edge(B, "z"), 
      edge(C, "N"),
      edge(D, "w"),
      admixture_edge("M", "x", "y"),
      admixture_edge("N", "z", "w")
    ))
    admixtures <- admixture_proportions(c(
      admix_props("M", "x", "y", "a"),
      admix_props("N", "z", "w", "b")
    ))
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  two_admixtures_12 = function(A, B, C, D) {
    leaves <- c(A, B, C, D)
    inner_nodes <- c("R", "x", "y", "z", "w", "M", "N")
    edges <- parent_edges(c(
      edge("x", "R"),
      edge("y", "x"),
      edge("z", "M"),
      edge("w", "z"),
      edge(A, "R"), 
      edge(B, "z"), 
      edge(C, "w"),
      edge(D, "N"),
      admixture_edge("M", "x", "y"),
      admixture_edge("N", "w", "y")
    ))
    admixtures <- admixture_proportions(c(
      admix_props("M", "x", "y", "a"),
      admix_props("N", "w", "y", "b")
    ))
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  two_admixtures_13 = function(A, B, C, D) {
    leaves <- c(A, B, C, D)
    inner_nodes <- c("R", "x", "y", "z", "w", "M", "N")
    edges <- parent_edges(c(
      edge("x", "R"),
      edge("y", "R"),
      edge("z", "x"),
      edge("w", "y"),
      edge(A, "x"), 
      edge(B, "M"), 
      edge(C, "N"),
      edge(D, "y"),
      admixture_edge("M", "z", "w"),
      admixture_edge("N", "z", "w")
    ))
    admixtures <- admixture_proportions(c(
      admix_props("M", "z", "w", "a"),
      admix_props("N", "z", "w", "b")
    ))
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  two_admixtures_14 = function(A, B, C, D) {
    leaves <- c(A, B, C, D)
    inner_nodes <- c("R", "x", "y", "z", "w", "M", "N")
    edges <- parent_edges(c(
      edge("x", "R"),
      edge("y", "x"),
      edge("z", "x"),
      edge("w", "z"),
      edge(A, "R"), 
      edge(B, "M"), 
      edge(C, "N"),
      edge(D, "w"),
      admixture_edge("M", "y", "z"),
      admixture_edge("N", "y", "w")
    ))
    admixtures <- admixture_proportions(c(
      admix_props("M", "y", "z", "a"),
      admix_props("N", "y", "w", "b")
    ))
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  two_admixtures_15 = function(A, B, C, D) {
    leaves <- c(A, B, C, D)
    inner_nodes <- c("R", "x", "y", "z", "w", "M", "N")
    edges <- parent_edges(c(
      edge("x", "R"),
      edge("y", "x"),
      edge("z", "x"),
      edge("w", "y"),
      edge(A, "R"), 
      edge(B, "w"), 
      edge(C, "M"),
      edge(D, "z"),
      admixture_edge("M", "y", "N"),
      admixture_edge("N", "w", "z")
    ))
    admixtures <- admixture_proportions(c(
      admix_props("M", "y", "N", "a"),
      admix_props("N", "w", "z", "b")
    ))
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  two_admixtures_16 = function(A, B, C, D) {
    leaves <- c(A, B, C, D)
    inner_nodes <- c("R", "x", "y", "z", "w", "M", "N")
    edges <- parent_edges(c(
      edge("x", "R"),
      edge("y", "x"),
      edge("z", "y"),
      edge("w", "y"),
      edge(A, "R"), 
      edge(B, "z"), 
      edge(C, "w"),
      edge(D, "N"),
      admixture_edge("M", "z", "w"),
      admixture_edge("N", "M", "x")
    ))
    admixtures <- admixture_proportions(c(
      admix_props("M", "z", "w", "a"),
      admix_props("N", "M", "x", "b")
    ))
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  two_admixtures_17 = function(A, B, C, D) {
    leaves <- c(A, B, C, D)
    inner_nodes <- c("R", "x", "y", "z", "w", "M", "N")
    edges <- parent_edges(c(
      edge("x", "R"),
      edge("y", "R"),
      edge("z", "y"),
      edge("w", "M"),
      edge(A, "x"), 
      edge(B, "w"), 
      edge(C, "N"),
      edge(D, "z"),
      admixture_edge("M", "x", "z"),
      admixture_edge("N", "w", "y")
    ))
    admixtures <- admixture_proportions(c(
      admix_props("M", "x", "z", "a"),
      admix_props("N", "w", "y", "b")
    ))
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  two_admixtures_18 = function(A, B, C, D) {
    leaves <- c(A, B, C, D)
    inner_nodes <- c("R", "x", "y", "z", "w", "M", "N")
    edges <- parent_edges(c(
      edge("x", "R"),
      edge("y", "x"),
      edge("z", "y"),
      edge("w", "y"),
      edge(A, "R"), 
      edge(B, "M"), 
      edge(C, "N"),
      edge(D, "x"),
      admixture_edge("M", "z", "w"),
      admixture_edge("N", "z", "w")
    ))
    admixtures <- admixture_proportions(c(
      admix_props("M", "z", "w", "a"),
      admix_props("N", "z", "w", "b")
    ))
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  two_admixtures_19 = function(A, B, C, D) {
    leaves <- c(A, B, C, D)
    inner_nodes <- c("R", "x", "y", "z", "w", "M", "N")
    edges <- parent_edges(c(
      edge("x", "R"),
      edge("y", "x"),
      edge("z", "y"),
      edge("w", "z"),
      edge(A, "R"), 
      edge(B, "w"), 
      edge(C, "N"),
      edge(D, "x"),
      admixture_edge("M", "w", "y"),
      admixture_edge("N", "z", "M")
    ))
    admixtures <- admixture_proportions(c(
      admix_props("M", "w", "y", "a"),
      admix_props("N", "z", "M", "b")
    ))
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  two_admixtures_20 = function(A, B, C, D) {
    leaves <- c(A, B, C, D)
    inner_nodes <- c("R", "x", "y", "z", "w", "M", "N")
    edges <- parent_edges(c(
      edge("x", "R"),
      edge("y", "x"),
      edge("z", "x"),
      edge("w", "N"),
      edge(A, "R"), 
      edge(B, "M"), 
      edge(C, "w"),
      edge(D, "w"),
      admixture_edge("M", "y", "z"),
      admixture_edge("N", "y", "z")
    ))
    admixtures <- admixture_proportions(c(
      admix_props("M", "y", "z", "a"),
      admix_props("N", "y", "z", "b")
    ))
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  two_admixtures_21 = function(A, B, C, D) {
    leaves <- c(A, B, C, D)
    inner_nodes <- c("R", "x", "y", "z", "w", "M", "N")
    edges <- parent_edges(c(
      edge("x", "R"),
      edge("y", "x"),
      edge("z", "y"),
      edge("w", "N"),
      edge(A, "R"), 
      edge(B, "z"), 
      edge(C, "w"),
      edge(D, "w"),
      admixture_edge("M", "z", "x"),
      admixture_edge("N", "y", "M")
    ))
    admixtures <- admixture_proportions(c(
      admix_props("M", "z", "x", "a"),
      admix_props("N", "y", "M", "b")
    ))
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  two_admixtures_22 = function(A, B, C, D) {
    leaves <- c(A, B, C, D)
    inner_nodes <- c("R", "x", "y", "z", "w", "M", "N")
    edges <- parent_edges(c(
      edge("x", "R"),
      edge("y", "R"),
      edge("z", "y"),
      edge("w", "N"),
      edge(A, "x"), 
      edge(B, "w"), 
      edge(C, "w"),
      edge(D, "y"),
      admixture_edge("M", "x", "z"),
      admixture_edge("N", "M", "z")
    ))
    admixtures <- admixture_proportions(c(
      admix_props("M", "x", "z", "a"),
      admix_props("N", "M", "z", "b")
    ))
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  two_admixtures_23 = function(A, B, C, D) {
    leaves <- c(A, B, C, D)
    inner_nodes <- c("R", "x", "y", "z", "w", "M", "N")
    edges <- parent_edges(c(
      edge("x", "R"),
      edge("y", "R"),
      edge("z", "y"),
      edge("w", "N"),
      edge(A, "x"), 
      edge(B, "M"), 
      edge(C, "w"),
      edge(D, "w"),
      admixture_edge("M", "x", "z"),
      admixture_edge("N", "y", "z")
    ))
    admixtures <- admixture_proportions(c(
      admix_props("M", "x", "z", "a"),
      admix_props("N", "y", "z", "b")
    ))
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  two_admixtures_24 = function(A, B, C, D) {
    leaves <- c(A, B, C, D)
    inner_nodes <- c("R", "x", "y", "z", "w", "M", "N")
    edges <- parent_edges(c(
      edge("x", "R"),
      edge("y", "x"),
      edge("z", "y"),
      edge("w", "z"),
      edge(A, "R"), 
      edge(B, "N"), 
      edge(C, "w"),
      edge(D, "x"),
      admixture_edge("M", "y", "z"),
      admixture_edge("N", "M", "w")
    ))
    admixtures <- admixture_proportions(c(
      admix_props("M", "y", "z", "a"),
      admix_props("N", "M", "w", "b")
    ))
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  two_admixtures_25 = function(A, B, C, D) {
    leaves <- c(A, B, C, D)
    inner_nodes <- c("R", "x", "y", "z", "w", "M", "N")
    edges <- parent_edges(c(
      edge("x", "R"),
      edge("y", "x"),
      edge("z", "y"),
      edge("w", "y"),
      edge(A, "R"), 
      edge(B, "z"), 
      edge(C, "N"),
      edge(D, "x"),
      admixture_edge("M", "z", "w"),
      admixture_edge("N", "M", "w")
    ))
    admixtures <- admixture_proportions(c(
      admix_props("M", "z", "w", "a"),
      admix_props("N", "M", "w", "b")
    ))
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  two_admixtures_26 = function(A, B, C, D) {
    leaves <- c(A, B, C, D)
    inner_nodes <- c("R", "x", "y", "z", "w", "M", "N")
    edges <- parent_edges(c(
      edge("x", "R"),
      edge("y", "R"),
      edge("z", "y"),
      edge("w", "M"),
      edge(A, "x"), 
      edge(B, "w"), 
      edge(C, "w"),
      edge(D, "N"),
      admixture_edge("M", "x", "z"),
      admixture_edge("N", "z", "y")
    ))
    admixtures <- admixture_proportions(c(
      admix_props("M", "x", "z", "a"),
      admix_props("N", "z", "y", "b")
    ))
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  two_admixtures_27 = function(A, B, C, D) {
    leaves <- c(A, B, C, D)
    inner_nodes <- c("R", "x", "y", "z", "w", "M", "N")
    edges <- parent_edges(c(
      edge("x", "R"),
      edge("y", "x"),
      edge("z", "y"),
      edge("w", "z"),
      edge(A, "R"), 
      edge(B, "M"), 
      edge(C, "N"),
      edge(D, "x"),
      admixture_edge("M", "y", "w"),
      admixture_edge("N", "w", "z")
    ))
    admixtures <- admixture_proportions(c(
      admix_props("M", "y", "w", "a"),
      admix_props("N", "w", "z", "b")
    ))
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  two_admixtures_28 = function(A, B, C, D) {
    leaves <- c(A, B, C, D)
    inner_nodes <- c("R", "x", "y", "z", "w", "M", "N")
    edges <- parent_edges(c(
      edge("x", "R"),
      edge("y", "R"),
      edge("z", "y"),
      edge("w", "N"),
      edge(A, "x"), 
      edge(B, "w"), 
      edge(C, "w"),
      edge(D, "z"),
      admixture_edge("M", "x", "y"),
      admixture_edge("N", "M", "z")
    ))
    admixtures <- admixture_proportions(c(
      admix_props("M", "x", "y", "a"),
      admix_props("N", "M", "z", "b")
    ))
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  two_admixtures_29 = function(A, B, C, D) {
    leaves <- c(A, B, C, D)
    inner_nodes <- c("R", "x", "y", "z", "w", "M", "N")
    edges <- parent_edges(c(
      edge("x", "R"),
      edge("y", "R"),
      edge("z", "y"),
      edge("w", "z"),
      edge(A, "x"), 
      edge(B, "N"), 
      edge(C, "w"),
      edge(D, "w"),
      admixture_edge("M", "x", "y"),
      admixture_edge("N", "M", "z")
    ))
    admixtures <- admixture_proportions(c(
      admix_props("M", "x", "y", "a"),
      admix_props("N", "M", "z", "b")
    ))
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  two_admixtures_30 = function(A, B, C, D) {
    leaves <- c(A, B, C, D)
    inner_nodes <- c("R", "x", "y", "z", "w", "M", "N")
    edges <- parent_edges(c(
      edge("x", "R"),
      edge("y", "R"),
      edge("z", "M"),
      edge("w", "N"),
      edge(A, "x"), 
      edge(B, "z"), 
      edge(C, "w"),
      edge(D, "w"),
      admixture_edge("M", "x", "y"),
      admixture_edge("N", "z", "y")
    ))
    admixtures <- admixture_proportions(c(
      admix_props("M", "x", "y", "a"),
      admix_props("N", "z", "y", "b")
    ))
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  two_admixtures_31 = function(A, B, C, D) {
    leaves <- c(A, B, C, D)
    inner_nodes <- c("R", "x", "y", "z", "w", "M", "N")
    edges <- parent_edges(c(
      edge("x", "R"),
      edge("y", "x"),
      edge("z", "y"),
      edge("w", "z"),
      edge(A, "R"), 
      edge(B, "w"), 
      edge(C, "N"),
      edge(D, "x"),
      admixture_edge("M", "w", "z"),
      admixture_edge("N", "M", "y")
    ))
    admixtures <- admixture_proportions(c(
      admix_props("M", "w", "z", "a"),
      admix_props("N", "M", "y", "b")
    ))
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  two_admixtures_32 = function(A, B, C, D) {
    leaves <- c(A, B, C, D)
    inner_nodes <- c("R", "x", "y", "z", "w", "M", "N")
    edges <- parent_edges(c(
      edge("x", "R"),
      edge("y", "R"),
      edge("z", "M"),
      edge("w", "N"),
      edge(A, "z"), 
      edge(B, "z"), 
      edge(C, "w"),
      edge(D, "w"),
      admixture_edge("M", "x", "y"),
      admixture_edge("N", "x", "y")
    ))
    admixtures <- admixture_proportions(c(
      admix_props("M", "x", "y", "a"),
      admix_props("N", "x", "y", "b")
    ))
    agraph(leaves, inner_nodes, edges, admixtures)
  }
)

five_leaves_graphs <- list(
  # tree
  tree = function(A, B, C, D, E) {
    leaves <- c(A, B, C, D, E)
    inner_nodes <- c("R", "x", "y", "z")
    edges <- parent_edges(c(
      edge("x", "R"),
      edge("y", "x"),
      edge("z", "y"),
      edge(A, "R"),
      edge(B, "x"),
      edge(C, "y"),
      edge(D, "z"),
      edge(E, "z")
    ))
    admixtures <- NULL
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  # one admixture
  one_admixture_1 = function(A, B, C, D, E) {
    leaves <- c(A, B, C, D, E)
    inner_nodes <- c("R", "x", "y", "z", "w", "M")
    edges <- parent_edges(c(
      edge("x", "R"),
      edge("y", "x"),
      edge("z", "x"),
      edge("w", "z"),
      edge(A, "R"),
      edge(B, "y"), 
      edge(C, "M"),
      edge(D, "w"),
      edge(E, "z"),
      admixture_edge("M", "y", "w")
    ))
    admixtures <- admixture_proportions(c(
      admix_props("M", "y", "w", "a")
    ))
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  one_admixture_2 = function(A, B, C, D, E) {
    leaves <- c(A, B, C, D, E)
    inner_nodes <- c("R", "x", "y", "z", "w", "M")
    edges <- parent_edges(c(
      edge("x", "R"),
      edge("y", "x"),
      edge("z", "x"),
      edge("w", "M"),
      edge(A, "R"),
      edge(B, "y"), 
      edge(C, "w"),
      edge(D, "w"),
      edge(E, "z"),
      admixture_edge("M", "y", "z")
    ))
    admixtures <- admixture_proportions(c(
      admix_props("M", "y", "z", "a")
    ))
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  one_admixture_3 = function(A, B, C, D, E) {
    leaves <- c(A, B, C, D, E)
    inner_nodes <- c("R", "x", "y", "z", "w", "M")
    edges <- parent_edges(c(
      edge("x", "R"),
      edge("y", "x"),
      edge("z", "x"),
      edge("w", "z"),
      edge(A, "R"),
      edge(B, "y"), 
      edge(C, "M"),
      edge(D, "w"),
      edge(E, "w"),
      admixture_edge("M", "y", "z")
    ))
    admixtures <- admixture_proportions(c(
      admix_props("M", "y", "z", "a")
    ))
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  one_admixture_4 = function(A, B, C, D, E) {
    leaves <- c(A, B, C, D, E)
    inner_nodes <- c("R", "x", "y", "z", "w", "M")
    edges <- parent_edges(c(
      edge("x", "R"),
      edge("y", "R"),
      edge("z", "x"),
      edge("w", "y"),
      edge(A, "z"),
      edge(B, "z"), 
      edge(C, "M"),
      edge(D, "w"),
      edge(E, "w"),
      admixture_edge("M", "x", "y")
    ))
    admixtures <- admixture_proportions(c(
      admix_props("M", "x", "y", "a")
    ))
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  one_admixture_5 = function(A, B, C, D, E) {
    leaves <- c(A, B, C, D, E)
    inner_nodes <- c("R", "x", "y", "z", "w", "M")
    edges <- parent_edges(c(
      edge("x", "R"),
      edge("y", "R"),
      edge("z", "M"),
      edge("w", "y"),
      edge(A, "x"),
      edge(B, "z"), 
      edge(C, "z"),
      edge(D, "w"),
      edge(E, "w"),
      admixture_edge("M", "x", "y")
    ))
    admixtures <- admixture_proportions(c(
      admix_props("M", "x", "y", "a")
    ))
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  one_admixture_6 = function(A, B, C, D, E) {
    leaves <- c(A, B, C, D, E)
    inner_nodes <- c("R", "x", "y", "z", "w", "M")
    edges <- parent_edges(c(
      edge("x", "R"),
      edge("y", "R"),
      edge("z", "M"),
      edge("w", "z"),
      edge(A, "x"),
      edge(B, "z"), 
      edge(C, "w"),
      edge(D, "w"),
      edge(E, "y"),
      admixture_edge("M", "x", "y")
    ))
    admixtures <- admixture_proportions(c(
      admix_props("M", "x", "y", "a")
    ))
    agraph(leaves, inner_nodes, edges, admixtures)
  },
  
  one_admixture_7 = function(A, B, C, D, E) {
    leaves <- c(A, B, C, D, E)
    inner_nodes <- c("R", "x", "y", "z", "w", "M")
    edges <- parent_edges(c(
      edge("x", "R"),
      edge("y", "R"),
      edge("z", "y"),
      edge("w", "z"),
      edge(A, "x"),
      edge(B, "M"), 
      edge(C, "z"),
      edge(D, "w"),
      edge(E, "w"),
      admixture_edge("M", "x", "y")
    ))
    admixtures <- admixture_proportions(c(
      admix_props("M", "x", "y", "a")
    ))
    agraph(leaves, inner_nodes, edges, admixtures)
  }
)



#' Combine a list of permutations with a list of (parameterized) graphs to fit them
#' 
#' @param permutations List of population permutations
#' @param graphs       List of functions for producing graphs
#' 
#' @return a list of fitted graphs.
fit_permutations_and_graphs <- function(data, permutations, graphs) {
  result <- vector("list", length(permutations) * length(graphs))
  idx <- 1
  for (permutation in permutations) {
    for (graph_constructure in graphs) {
      graph <- do.call(graph_constructure, as.list(permutation))
      fit <- fit_graph(filter_on_leaves(data, graph), graph)
      result[[idx]] <- fit
      idx <- idx + 1
    }
  }
  result
}

#' Helping to untangle the relationships between four populations based on given data.
#' 
#' Note that we do not say anything about the exact position of the root: the graphs
#' drawn are just chosen in a way they look nice. In all but one cases it is possible
#' to set the root on an outgroup branch, but sometimes the plot is so weird that we
#' nonetheless prefer to set it somewhere else. The first element of populations is
#' the leaf often drawn as the outgroup.
#'
#' @param data         The data set.
#' @param populations  A four or element vector of population names.
#'
#' @return Printing stuff on console and "plots" tab. Making something fancy later.
#'
#' @export
fit_all_graphs <- function(data, populations) {
  if (!(length(populations) %in% c(4,5)) ) {
    stop("We can currently only explore graph spaces for four or five leaves")
  }
    
  P <- make_permutations(populations)
  graphs <- if (length(populations) == 4) four_leaves_graphs else five_leaves_graphs
  fits <- fit_permutations_and_graphs(data, P, graphs)
  structure(fits, class = c("agraph_fit_list", "list"))
}



