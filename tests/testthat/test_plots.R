context("Graph plotting")

test_that("trees plot fine.", {
  leaves <- c("A")
  inner_nodes <- c("Root")
  edges <- parent_edges(c(edge("A", "Root")))
  admix <- NULL
  graph <- agraph(leaves, inner_nodes, edges, admix)
  plot(graph)
  
  leaves <- c("A", "B", "C")
  inner_nodes <- c("Root", "z", "y")
  edges <- parent_edges(c(edge("z", "Root"),
                          edge("y", "z"),
                          edge("A", "y"),
                          edge("B", "y"),
                          edge("C", "y")))
  admix <- NULL
  graph <- agraph(leaves, inner_nodes, edges, admix)
  plot(graph)
})

test_that("complex non-planar grahs plot fine.", {
  leaves <- c("A", "B")
  inner_nodes <- c("Root", "L", "R")
  edges <- parent_edges(c(edge("L", "Root"),
                          edge("R", "Root"),
                          admixture_edge("A", "L", "R"),
                          admixture_edge("B", "L", "R")))
  admix <- admixture_proportions(c(admix_props("A", "L", "R", "a"),
                                   admix_props("B", "L", "R", "b")))
  graph <- agraph(leaves, inner_nodes, edges, admix)
  plot(graph)
  
  leaves <- c("A")
  inner_nodes <- c("Root", "L1", "R1", "L2", "R2")
  edges <- parent_edges(c(edge("L1", "Root"),
                          edge("R1", "Root"),
                          admixture_edge("L2", "L1", "R1"),
                          admixture_edge("R2", "L1", "R1"),
                          admixture_edge("A", "L2", "R2")))
  admix <- admixture_proportions(c(admix_props("L2", "L1", "R1", "l"),
                                   admix_props("R2", "L1", "R1", "r"),
                                   admix_props("A", "L2", "R2", "a")))
  graph <- agraph(leaves, inner_nodes, edges, admix)
  plot(graph, platform = 0)
  
  leaves <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J")
  inner_nodes <- c("Root", "L1", "R1", "L2", "R2", "Z")
  edges <- parent_edges(c(edge("L1", "Root"),
                          edge("R1", "Root"),
                          admixture_edge("L2", "L1", "R1"),
                          admixture_edge("R2", "L1", "R1"),
                          admixture_edge("Z", "L2", "R2"),
                          edge("A", "Root"),
                          edge("B", "Root"),
                          edge("C", "Root"),
                          edge("D", "L1"),
                          edge("E", "R1"),
                          edge("F", "L2"),
                          edge("G", "R2"),
                          edge("H", "Z"),
                          edge("I", "Z"),
                          edge("J", "Z")))
  admix <- admixture_proportions(c(admix_props("L2", "L1", "R1", "l"),
                                   admix_props("R2", "L1", "R1", "r"),
                                   admix_props("Z", "L2", "R2", "z")))
  graph <- agraph(leaves, inner_nodes, edges, admix)
  plot(graph)
})