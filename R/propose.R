#' Helping to untangle the relationships between four populations based on given data.
#' 
#' Note that we do not say anything about the exact position of the root: the graphs
#' drawn are just chosen in a way they look nice. In all but one cases it is possible
#' to set the root on an outgroup branch, but sometimes the plot is so weird that we
#' nonetheless prefer to set it somewhere else. The first element of populations is
#' the leaf often drawn as the outgroup.
#'
#' @param data         The data set.
#' @param populations  A four element vector of population names.
#'
#' @return Printing stuff on console and "plots" tab. Making something fancy later.
#'
#' @export
fit_four <- function(data, populations) {
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
  best_tree_error <- Inf
  tree_errors <- numeric()
  best_one_admixture_error <- Inf
  one_admixture_errors <- numeric()
  best_two_admixture_error <- Inf
  two_admixture_errors <- numeric()
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
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    plot(graph, main = fit$best_error)
    tree_errors[length(tree_errors) + 1] <- fit$best_error
    if (fit$best_error < best_tree_error) {
      best_tree <- graph
      best_tree_error <- fit$best_error
    }
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
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    plot(graph, main = fit$best_error)
    one_admixture_errors[length(one_admixture_errors) + 1] <- fit$best_error
    if (fit$best_error <= best_one_admixture_error) {
      best_one_admixture <- graph
      best_one_admixture_error <- fit$best_error
    }
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
    graph <- agraph(leaves, inner_nodes, edges, admixtures)
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    plot(graph, main = fit$best_error)
    one_admixture_errors[length(one_admixture_errors) + 1] <- fit$best_error
    if (fit$best_error < best_one_admixture_error) {
      best_one_admixture <- graph
      best_one_admixture_error <- fit$best_error
    }
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
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    plot(graph, main = fit$best_error)
    one_admixture_errors[length(one_admixture_errors) + 1] <- fit$best_error
    if (fit$best_error < best_one_admixture_error) {
      best_one_admixture <- graph
      best_one_admixture_error <- fit$best_error
    }
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
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    plot(graph, main = fit$best_error)
    two_admixture_errors[length(two_admixture_errors) + 1] <- fit$best_error
    if (fit$best_error < best_two_admixture_error) {
      best_two_admixture <- graph
      best_two_admixture_error <- fit$best_error
    }
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
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    plot(graph, main = fit$best_error)
    two_admixture_errors[length(two_admixture_errors) + 1] <- fit$best_error
    if (fit$best_error < best_two_admixture_error) {
      best_two_admixture <- graph
      best_two_admixture_error <- fit$best_error
    }
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
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    plot(graph, main = fit$best_error)
    two_admixture_errors[length(two_admixture_errors) + 1] <- fit$best_error
    if (fit$best_error < best_two_admixture_error) {
      best_two_admixture <- graph
      best_two_admixture_error <- fit$best_error
    }
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
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    plot(graph, main = fit$best_error)
    two_admixture_errors[length(two_admixture_errors) + 1] <- fit$best_error
    if (fit$best_error < best_two_admixture_error) {
      best_two_admixture <- graph
      best_two_admixture_error <- fit$best_error
    }
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
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    plot(graph, main = fit$best_error)
    two_admixture_errors[length(two_admixture_errors) + 1] <- fit$best_error
    if (fit$best_error < best_two_admixture_error) {
      best_two_admixture <- graph
      best_two_admixture_error <- fit$best_error
    }
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
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    plot(graph, main = fit$best_error)
    two_admixture_errors[length(two_admixture_errors) + 1] <- fit$best_error
    if (fit$best_error < best_two_admixture_error) {
      best_two_admixture <- graph
      best_two_admixture_error <- fit$best_error
    }
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
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    plot(graph, main = fit$best_error)
    two_admixture_errors[length(two_admixture_errors) + 1] <- fit$best_error
    if (fit$best_error < best_two_admixture_error) {
      best_two_admixture <- graph
      best_two_admixture_error <- fit$best_error
    }
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
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    plot(graph, main = fit$best_error)
    two_admixture_errors[length(two_admixture_errors) + 1] <- fit$best_error
    if (fit$best_error < best_two_admixture_error) {
      best_two_admixture <- graph
      best_two_admixture_error <- fit$best_error
    }
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
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    plot(graph, main = fit$best_error)
    two_admixture_errors[length(two_admixture_errors) + 1] <- fit$best_error
    if (fit$best_error < best_two_admixture_error) {
      best_two_admixture <- graph
      best_two_admixture_error <- fit$best_error
    }
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
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    plot(graph, main = fit$best_error)
    two_admixture_errors[length(two_admixture_errors) + 1] <- fit$best_error
    if (fit$best_error < best_two_admixture_error) {
      best_two_admixture <- graph
      best_two_admixture_error <- fit$best_error
    }
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
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    plot(graph, main = fit$best_error)
    two_admixture_errors[length(two_admixture_errors) + 1] <- fit$best_error
    if (fit$best_error < best_two_admixture_error) {
      best_two_admixture <- graph
      best_two_admixture_error <- fit$best_error
    }
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
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    plot(graph, main = fit$best_error)
    two_admixture_errors[length(two_admixture_errors) + 1] <- fit$best_error
    if (fit$best_error < best_two_admixture_error) {
      best_two_admixture <- graph
      best_two_admixture_error <- fit$best_error
    }
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
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    plot(graph, main = fit$best_error)
    two_admixture_errors[length(two_admixture_errors) + 1] <- fit$best_error
    if (fit$best_error < best_two_admixture_error) {
      best_two_admixture <- graph
      best_two_admixture_error <- fit$best_error
    }
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
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    plot(graph, main = fit$best_error)
    two_admixture_errors[length(two_admixture_errors) + 1] <- fit$best_error
    if (fit$best_error < best_two_admixture_error) {
      best_two_admixture <- graph
      best_two_admixture_error <- fit$best_error
    }
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
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    plot(graph, main = fit$best_error)
    two_admixture_errors[length(two_admixture_errors) + 1] <- fit$best_error
    if (fit$best_error < best_two_admixture_error) {
      best_two_admixture <- graph
      best_two_admixture_error <- fit$best_error
    }
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
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    plot(graph, main = fit$best_error)
    two_admixture_errors[length(two_admixture_errors) + 1] <- fit$best_error
    if (fit$best_error < best_two_admixture_error) {
      best_two_admixture <- graph
      best_two_admixture_error <- fit$best_error
    }
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
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    plot(graph, main = fit$best_error)
    two_admixture_errors[length(two_admixture_errors) + 1] <- fit$best_error
    if (fit$best_error < best_two_admixture_error) {
      best_two_admixture <- graph
      best_two_admixture_error <- fit$best_error
    }
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
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    plot(graph, main = fit$best_error)
    two_admixture_errors[length(two_admixture_errors) + 1] <- fit$best_error
    if (fit$best_error < best_two_admixture_error) {
      best_two_admixture <- graph
      best_two_admixture_error <- fit$best_error
    }
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
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    plot(graph, main = fit$best_error)
    two_admixture_errors[length(two_admixture_errors) + 1] <- fit$best_error
    if (fit$best_error < best_two_admixture_error) {
      best_two_admixture <- graph
      best_two_admixture_error <- fit$best_error
    }
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
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    plot(graph, main = fit$best_error)
    two_admixture_errors[length(two_admixture_errors) + 1] <- fit$best_error
    if (fit$best_error < best_two_admixture_error) {
      best_two_admixture <- graph
      best_two_admixture_error <- fit$best_error
    }
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
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    plot(graph, main = fit$best_error)
    two_admixture_errors[length(two_admixture_errors) + 1] <- fit$best_error
    if (fit$best_error < best_two_admixture_error) {
      best_two_admixture <- graph
      best_two_admixture_error <- fit$best_error
    }
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
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    plot(graph, main = fit$best_error)
    two_admixture_errors[length(two_admixture_errors) + 1] <- fit$best_error
    if (fit$best_error < best_two_admixture_error) {
      best_two_admixture <- graph
      best_two_admixture_error <- fit$best_error
    }
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
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    plot(graph, main = fit$best_error)
    two_admixture_errors[length(two_admixture_errors) + 1] <- fit$best_error
    if (fit$best_error < best_two_admixture_error) {
      best_two_admixture <- graph
      best_two_admixture_error <- fit$best_error
    }
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
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    plot(graph, main = fit$best_error)
    two_admixture_errors[length(two_admixture_errors) + 1] <- fit$best_error
    if (fit$best_error < best_two_admixture_error) {
      best_two_admixture <- graph
      best_two_admixture_error <- fit$best_error
    }
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
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    plot(graph, main = fit$best_error)
    two_admixture_errors[length(two_admixture_errors) + 1] <- fit$best_error
    if (fit$best_error < best_two_admixture_error) {
      best_two_admixture <- graph
      best_two_admixture_error <- fit$best_error
    }
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
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    plot(graph, main = fit$best_error)
    two_admixture_errors[length(two_admixture_errors) + 1] <- fit$best_error
    if (fit$best_error < best_two_admixture_error) {
      best_two_admixture <- graph
      best_two_admixture_error <- fit$best_error
    }
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
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    plot(graph, main = fit$best_error)
    two_admixture_errors[length(two_admixture_errors) + 1] <- fit$best_error
    if (fit$best_error < best_two_admixture_error) {
      best_two_admixture <- graph
      best_two_admixture_error <- fit$best_error
    }
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
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    plot(graph, main = fit$best_error)
    two_admixture_errors[length(two_admixture_errors) + 1] <- fit$best_error
    if (fit$best_error < best_two_admixture_error) {
      best_two_admixture <- graph
      best_two_admixture_error <- fit$best_error
    }
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
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    plot(graph, main = fit$best_error)
    two_admixture_errors[length(two_admixture_errors) + 1] <- fit$best_error
    if (fit$best_error < best_two_admixture_error) {
      best_two_admixture <- graph
      best_two_admixture_error <- fit$best_error
    }
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
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    plot(graph, main = fit$best_error)
    two_admixture_errors[length(two_admixture_errors) + 1] <- fit$best_error
    if (fit$best_error < best_two_admixture_error) {
      best_two_admixture <- graph
      best_two_admixture_error <- fit$best_error
    }
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
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    plot(graph, main = fit$best_error)
    two_admixture_errors[length(two_admixture_errors) + 1] <- fit$best_error
    if (fit$best_error < best_two_admixture_error) {
      best_two_admixture <- graph
      best_two_admixture_error <- fit$best_error
    }
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
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    plot(graph, main = fit$best_error)
    two_admixture_errors[length(two_admixture_errors) + 1] <- fit$best_error
    if (fit$best_error < best_two_admixture_error) {
      best_two_admixture <- graph
      best_two_admixture_error <- fit$best_error
    }
  }
  plot(best_tree, main = "The best tree")
  print("Best error among trees (the graph is plotted):")
  print(best_tree_error)
  print("The other errors for comparison:")
  print(tree_errors)
  plot(best_one_admixture, main = "The best graph with one admixture event")
  print("Best error allowing one admixture event (the graph is plotted):")
  print(best_one_admixture_error)
  print("The other errors for comparison:")
  print(one_admixture_errors)
  plot(best_two_admixture, main = "The best graph with two admixture events")
  print("Best error allowing two admixture events (the graph is plotted):")
  print(best_two_admixture_error)
  print("The other errors for comparison:")
  print(two_admixture_errors)
}

#' Helping to untangle the relationships between five populations based on given data.
#' 
#' Note that we do not say anything about the exact position of the root: the graphs
#' drawn are just chosen in a way they look nice. The first element of populations is
#' the leaf often drawn as the outgroup.
#'
#' @param data         The data set.
#' @param populations  A five element vector of population names.
#'
#' @return Printing stuff on console and "plots" tab. Making something fancy later.
#'
#' @export
fit_five <- function(data, populations) {
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
  best_tree_error <- Inf
  tree_errors <- numeric()
  best_one_admixture_error <- Inf
  one_admixture_errors <- numeric()
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
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    plot(graph, main = fit$best_error)
    tree_errors[length(tree_errors) + 1] <- fit$best_error
    if (fit$best_error < best_tree_error) {
      best_tree <- graph
      best_tree_error <- fit$best_error
    }
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
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    plot(graph, main = fit$best_error)
    one_admixture_errors[length(one_admixture_errors) + 1] <- fit$best_error
    if (fit$best_error <= best_one_admixture_error) {
      best_one_admixture <- graph
      best_one_admixture_error <- fit$best_error
    }
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
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    plot(graph, main = fit$best_error)
    one_admixture_errors[length(one_admixture_errors) + 1] <- fit$best_error
    if (fit$best_error <= best_one_admixture_error) {
      best_one_admixture <- graph
      best_one_admixture_error <- fit$best_error
    }
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
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    plot(graph, main = fit$best_error)
    one_admixture_errors[length(one_admixture_errors) + 1] <- fit$best_error
    if (fit$best_error <= best_one_admixture_error) {
      best_one_admixture <- graph
      best_one_admixture_error <- fit$best_error
    }
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
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    plot(graph, main = fit$best_error)
    one_admixture_errors[length(one_admixture_errors) + 1] <- fit$best_error
    if (fit$best_error <= best_one_admixture_error) {
      best_one_admixture <- graph
      best_one_admixture_error <- fit$best_error
    }
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
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    plot(graph, main = fit$best_error)
    one_admixture_errors[length(one_admixture_errors) + 1] <- fit$best_error
    if (fit$best_error <= best_one_admixture_error) {
      best_one_admixture <- graph
      best_one_admixture_error <- fit$best_error
    }
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
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    plot(graph, main = fit$best_error)
    one_admixture_errors[length(one_admixture_errors) + 1] <- fit$best_error
    if (fit$best_error <= best_one_admixture_error) {
      best_one_admixture <- graph
      best_one_admixture_error <- fit$best_error
    }
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
    fit <- fit_graph(filter_on_leaves(data, graph), graph)
    plot(graph, main = fit$best_error)
    one_admixture_errors[length(one_admixture_errors) + 1] <- fit$best_error
    if (fit$best_error <= best_one_admixture_error) {
      best_one_admixture <- graph
      best_one_admixture_error <- fit$best_error
    }
  }
  plot(best_tree, main = "The best tree")
  print("Best error among trees (the graph is plotted):")
  print(best_tree_error)
  print("The other errors for comparison:")
  print(tree_errors)
  plot(best_one_admixture, main = "The best graph with one admixture event")
  print("Best error allowing one admixture event (the graph is plotted):")
  print(best_one_admixture_error)
  print("The other errors for comparison:")
  print(one_admixture_errors)
}
