#' Helping to untangle the relationships between four populations based on given data.
#' 
#' Note that we do not say anything about the exact position of the root: the graphs
#' drawn are just chosen in a way they look nice. In all but one cases it is possible
#' to set the root on an outgroup branch, but often the plot is so weird that we
#' nonetheless prefer to set it somewhere else.
#'
#' @param data         The data set.
#' @param populations  A four element vector of population names.
#'
#' @return Printing stuff on console and "plots" tab. Making something fancy later.
#'
#' @export
fit_four <- function(data, populations, display = 10) {
  # We start by building a set of all 24 permutations of the 4 populations.
  P <- list()
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
  tree_list <- list()
  one_admixture_list <- list()
  two_admixture_list <- list()
  for (permutation in P) {
    A <- permutation[1]
    B <- permutation[2]
    C <- permutation[3]
    D <- permutation[4]
    # First the only tree.
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
    graph <- agraph(leaves, inner_nodes, edges, admixtures)
    description <- "Tree 1, each combination listed 8 times"
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    tree_list[[length(tree_list) + 1]] <- list(graph = graph, description = description,
                                               error = fit$best_error)
    # Then the three graphs with one admixture event.
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
    graph <- agraph(leaves, inner_nodes, edges, admixtures)
    description <- "Graph with one admixture 1, each combination listed 2 times"
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    one_admixture_list[[length(one_admixture_list) + 1]] <- list(graph = graph, description = description,
                                                                 error = fit$best_error)
    leaves <- c(A, B, C, D)
    inner_nodes <- c("R", "y", "z", "w", "M")
    edges <- parent_edges(c(
      edge("y", "R"),
      edge("z", "R"),
      edge("w", "z"),
      edge(A, "y"), 
      edge(B, "M"), 
      edge(C, "w"),
      edge(D, "w"),
      admixture_edge("M", "y", "z")
    ))
    admixtures <- admixture_proportions(c(
      admix_props("M", "y", "z", "a")
    ))
    graph <- agraph(leaves, inner_nodes, edges, admixtures)
    description <- "Graph with one admixture 2, each combination listed 2 times"
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    one_admixture_list[[length(one_admixture_list) + 1]] <- list(graph = graph, description = description,
                                                                 error = fit$best_error)
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
    graph <- agraph(leaves, inner_nodes, edges, admixtures)
    description <- "Graph with one admixture 3, each combination listed 4 times"
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    one_admixture_list[[length(one_admixture_list) + 1]] <- list(graph = graph, description = description,
                                                                 error = fit$best_error)
    # Finally the 32 graphs with two admixture events, in no particular order.
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
    graph <- agraph(leaves, inner_nodes, edges, admixtures)
    description <- "Graph with two admixtures 1, each combination listed 2 times"
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    two_admixture_list[[length(two_admixture_list) + 1]] <- list(graph = graph, description = description,
                                                                 error = fit$best_error)
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
    graph <- agraph(leaves, inner_nodes, edges, admixtures)
    description <- "Graph with two admixtures 2, each combination listed once"
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    two_admixture_list[[length(two_admixture_list) + 1]] <- list(graph = graph, description = description,
                                                                 error = fit$best_error)
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
    graph <- agraph(leaves, inner_nodes, edges, admixtures)
    description <- "Graph with two admixtures 3, each combination listed once"
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    two_admixture_list[[length(two_admixture_list) + 1]] <- list(graph = graph, description = description,
                                                                 error = fit$best_error)
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
    graph <- agraph(leaves, inner_nodes, edges, admixtures)
    description <- "Graph with two admixtures 4, each combination listed once"
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    two_admixture_list[[length(two_admixture_list) + 1]] <- list(graph = graph, description = description,
                                                                 error = fit$best_error)
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
    graph <- agraph(leaves, inner_nodes, edges, admixtures)
    description <- "Graph with two admixtures 5, each combination listed once"
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    two_admixture_list[[length(two_admixture_list) + 1]] <- list(graph = graph, description = description,
                                                                 error = fit$best_error)
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
    graph <- agraph(leaves, inner_nodes, edges, admixtures)
    description <- "Graph with two admixtures 6, each combination listed 2 times"
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    two_admixture_list[[length(two_admixture_list) + 1]] <- list(graph = graph, description = description,
                                                                 error = fit$best_error)
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
    graph <- agraph(leaves, inner_nodes, edges, admixtures)
    description <- "Graph with two admixtures 7, each combination listed once"
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    two_admixture_list[[length(two_admixture_list) + 1]] <- list(graph = graph, description = description,
                                                                 error = fit$best_error)
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
    graph <- agraph(leaves, inner_nodes, edges, admixtures)
    description <- "Graph with two admixtures 8, each combination listed once"
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    two_admixture_list[[length(two_admixture_list) + 1]] <- list(graph = graph, description = description,
                                                                 error = fit$best_error)
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
    graph <- agraph(leaves, inner_nodes, edges, admixtures)
    description <- "Graph with two admixtures 9, each combination listed once"
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    two_admixture_list[[length(two_admixture_list) + 1]] <- list(graph = graph, description = description,
                                                                 error = fit$best_error)
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
    graph <- agraph(leaves, inner_nodes, edges, admixtures)
    description <- "Graph with two admixtures 10, each combination listed once"
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    two_admixture_list[[length(two_admixture_list) + 1]] <- list(graph = graph, description = description,
                                                                 error = fit$best_error)
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
    graph <- agraph(leaves, inner_nodes, edges, admixtures)
    description <- "Graph with two admixtures 11, each combination listed once"
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    two_admixture_list[[length(two_admixture_list) + 1]] <- list(graph = graph, description = description,
                                                                 error = fit$best_error)
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
    graph <- agraph(leaves, inner_nodes, edges, admixtures)
    description <- "Graph with two admixtures 12, each combination listed once"
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    two_admixture_list[[length(two_admixture_list) + 1]] <- list(graph = graph, description = description,
                                                                 error = fit$best_error)
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
    graph <- agraph(leaves, inner_nodes, edges, admixtures)
    description <- "Graph with two admixtures 13, each combination listed 4 times"
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    two_admixture_list[[length(two_admixture_list) + 1]] <- list(graph = graph, description = description,
                                                                 error = fit$best_error)
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
    graph <- agraph(leaves, inner_nodes, edges, admixtures)
    description <- "Graph with two admixtures 14, each combination listed once"
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    two_admixture_list[[length(two_admixture_list) + 1]] <- list(graph = graph, description = description,
                                                                 error = fit$best_error)
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
    graph <- agraph(leaves, inner_nodes, edges, admixtures)
    description <- "Graph with two admixtures 15, each combination listed once"
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    two_admixture_list[[length(two_admixture_list) + 1]] <- list(graph = graph, description = description,
                                                                 error = fit$best_error)
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
    graph <- agraph(leaves, inner_nodes, edges, admixtures)
    description <- "Graph with two admixtures 16, each combination listed once"
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    two_admixture_list[[length(two_admixture_list) + 1]] <- list(graph = graph, description = description,
                                                                 error = fit$best_error)
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
    graph <- agraph(leaves, inner_nodes, edges, admixtures)
    description <- "Graph with two admixtures 17, each combination listed once"
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    two_admixture_list[[length(two_admixture_list) + 1]] <- list(graph = graph, description = description,
                                                                 error = fit$best_error)
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
    graph <- agraph(leaves, inner_nodes, edges, admixtures)
    description <- "Graph with two admixtures 18, each combination listed 4 times"
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    two_admixture_list[[length(two_admixture_list) + 1]] <- list(graph = graph, description = description,
                                                                 error = fit$best_error)
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
    graph <- agraph(leaves, inner_nodes, edges, admixtures)
    description <- "Graph with two admixtures 19, each combination listed 2 times"
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    two_admixture_list[[length(two_admixture_list) + 1]] <- list(graph = graph, description = description,
                                                                 error = fit$best_error)
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
    graph <- agraph(leaves, inner_nodes, edges, admixtures)
    description <- "Graph with two admixtures 20, each combination listed 2 times"
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    two_admixture_list[[length(two_admixture_list) + 1]] <- list(graph = graph, description = description,
                                                                 error = fit$best_error)
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
    graph <- agraph(leaves, inner_nodes, edges, admixtures)
    description <- "Graph with two admixtures 21, each combination listed 2 times"
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    two_admixture_list[[length(two_admixture_list) + 1]] <- list(graph = graph, description = description,
                                                                 error = fit$best_error)
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
    graph <- agraph(leaves, inner_nodes, edges, admixtures)
    description <- "Graph with two admixtures 22, each combination listed 2 times"
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    two_admixture_list[[length(two_admixture_list) + 1]] <- list(graph = graph, description = description,
                                                                 error = fit$best_error)
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
    graph <- agraph(leaves, inner_nodes, edges, admixtures)
    description <- "Graph with two admixtures 23, each combination listed 2 times"
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    two_admixture_list[[length(two_admixture_list) + 1]] <- list(graph = graph, description = description,
                                                                 error = fit$best_error)
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
    graph <- agraph(leaves, inner_nodes, edges, admixtures)
    description <- "Graph with two admixtures 24, each combination listed 2 times"
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    two_admixture_list[[length(two_admixture_list) + 1]] <- list(graph = graph, description = description,
                                                                 error = fit$best_error)
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
    graph <- agraph(leaves, inner_nodes, edges, admixtures)
    description <- "Graph with two admixtures 25, each combination listed 2 times"
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    two_admixture_list[[length(two_admixture_list) + 1]] <- list(graph = graph, description = description,
                                                                 error = fit$best_error)
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
    graph <- agraph(leaves, inner_nodes, edges, admixtures)
    description <- "Graph with two admixtures 26, each combination listed 2 times"
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    two_admixture_list[[length(two_admixture_list) + 1]] <- list(graph = graph, description = description,
                                                                 error = fit$best_error)
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
    graph <- agraph(leaves, inner_nodes, edges, admixtures)
    description <- "Graph with two admixtures 27, each combination listed 2 times"
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    two_admixture_list[[length(two_admixture_list) + 1]] <- list(graph = graph, description = description,
                                                                 error = fit$best_error)
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
    graph <- agraph(leaves, inner_nodes, edges, admixtures)
    description <- "Graph with two admixtures 28, each combination listed 2 times"
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    two_admixture_list[[length(two_admixture_list) + 1]] <- list(graph = graph, description = description,
                                                                 error = fit$best_error)
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
    graph <- agraph(leaves, inner_nodes, edges, admixtures)
    description <- "Graph with two admixtures 29, each combination listed 2 times"
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    two_admixture_list[[length(two_admixture_list) + 1]] <- list(graph = graph, description = description,
                                                                 error = fit$best_error)
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
    graph <- agraph(leaves, inner_nodes, edges, admixtures)
    description <- "Graph with two admixtures 30, each combination listed 2 times"
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    two_admixture_list[[length(two_admixture_list) + 1]] <- list(graph = graph, description = description,
                                                                 error = fit$best_error)
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
    graph <- agraph(leaves, inner_nodes, edges, admixtures)
    description <- "Graph with two admixtures 31, each combination listed 2 times"
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    two_admixture_list[[length(two_admixture_list) + 1]] <- list(graph = graph, description = description,
                                                                 error = fit$best_error)
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
    graph <- agraph(leaves, inner_nodes, edges, admixtures)
    description <- "Graph with two admixtures 32, each combination listed 8 times"
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    two_admixture_list[[length(two_admixture_list) + 1]] <- list(graph = graph, description = description,
                                                                 error = fit$best_error)
  }
  tree_list <- tree_list[order(sapply(tree_list, "[[", 3))]
  one_admixture_list <- one_admixture_list[order(sapply(one_admixture_list, "[[", 3))]
  two_admixture_list <- two_admixture_list[order(sapply(two_admixture_list, "[[", 3))]
  if (display > 0) {
    for (j in seq(1, display)) {
      stuff <- paste(tree_list[[j]]$description, ", error: ", tree_list[[j]]$error, sep = "")
      plot(tree_list[[j]]$graph, main = stuff)
    }
    for (j in seq(1, display)) {
      stuff <- paste(one_admixture_list[[j]]$description, ", error: ", one_admixture_list[[j]]$error, sep = "")
      plot(one_admixture_list[[j]]$graph, main = stuff)
    }
    for (j in seq(1, display)) {
      stuff <- paste(two_admixture_list[[j]]$description, ", error: ", two_admixture_list[[j]]$error, sep = "")
      plot(two_admixture_list[[j]]$graph, main = stuff)
    }
  }
  list(tree = tree_list, one_admixture = one_admixture_list, two_admixture = two_admixture_list)
}

#' Helping to untangle the relationships between five populations based on given data.
#' 
#' Note that we do not say anything about the exact position of the root: the graphs
#' drawn are just chosen in a way they look nice.
#'
#' @param data         The data set.
#' @param populations  A five element vector of population names.
#'
#' @return Printing stuff on console and "plots" tab. Making something fancy later.
#'
#' @export
fit_five <- function(data, populations, display = 10) {
  # We start by building a set of all 120 permutations of the 5 populations.
  P <- list()
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
  tree_list <- list()
  one_admixture_list <- list()
  for (permutation in P) {
    A <- permutation[1]
    B <- permutation[2]
    C <- permutation[3]
    D <- permutation[4]
    E <- permutation[5]
    # First the only tree.
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
    graph <- agraph(leaves, inner_nodes, edges, admixtures)
    description <- "Tree 1, each combination listed 8 times"
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    tree_list[[length(tree_list) + 1]] <- list(graph = graph, description = description,
                                               error = fit$best_error)
    # Then the seven graphs with one admixture event.
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
    graph <- agraph(leaves, inner_nodes, edges, admixtures)
    description <- "Graph with one admixture 1, each combination listed 2 times"
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    one_admixture_list[[length(one_admixture_list) + 1]] <- list(graph = graph, description = description,
                                                                 error = fit$best_error)
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
    graph <- agraph(leaves, inner_nodes, edges, admixtures)
    description <- "Graph with one admixture 2, each combination listed 4 times"
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    one_admixture_list[[length(one_admixture_list) + 1]] <- list(graph = graph, description = description,
                                                                 error = fit$best_error)
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
    graph <- agraph(leaves, inner_nodes, edges, admixtures)
    description <- "Graph with one admixture 3, each combination listed 2 times"
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    one_admixture_list[[length(one_admixture_list) + 1]] <- list(graph = graph, description = description,
                                                                 error = fit$best_error)
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
    graph <- agraph(leaves, inner_nodes, edges, admixtures)
    description <- "Graph with one admixture 4, each combination listed 8 times"
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    one_admixture_list[[length(one_admixture_list) + 1]] <- list(graph = graph, description = description,
                                                                 error = fit$best_error)
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
    graph <- agraph(leaves, inner_nodes, edges, admixtures)
    description <- "Graph with one admixture 5, each combination listed 4 times"
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    one_admixture_list[[length(one_admixture_list) + 1]] <- list(graph = graph, description = description,
                                                                 error = fit$best_error)
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
    graph <- agraph(leaves, inner_nodes, edges, admixtures)
    description <- "Graph with one admixture 6, each combination listed 4 times"
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    one_admixture_list[[length(one_admixture_list) + 1]] <- list(graph = graph, description = description,
                                                                 error = fit$best_error)
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
    graph <- agraph(leaves, inner_nodes, edges, admixtures)
    description <- "Graph with one admixture 7, each combination listed 2 times"
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    one_admixture_list[[length(one_admixture_list) + 1]] <- list(graph = graph, description = description,
                                                                 error = fit$best_error)
  }
  tree_list <- tree_list[order(sapply(tree_list, "[[", 3))]
  one_admixture_list <- one_admixture_list[order(sapply(one_admixture_list, "[[", 3))]
  if (display > 0) {
    for (j in seq(1, display)) {
      stuff <- paste(tree_list[[j]]$description, ", error: ", tree_list[[j]]$error, sep = "")
      plot(tree_list[[j]]$graph, main = stuff)
    }
    for (j in seq(1, display)) {
      stuff <- paste(one_admixture_list[[j]]$description, ", error: ", one_admixture_list[[j]]$error, sep = "")
      plot(one_admixture_list[[j]]$graph, main = stuff)
    }
  }
  list(tree = tree_list, one_admixture = one_admixture_list)
}