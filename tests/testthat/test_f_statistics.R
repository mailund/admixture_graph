context("f statistics")

test_that("that the f statistics computes the correct paths overlaps.", {
  leaves <- c("A", "B", "C", "D")
  inner_nodes <- c("AC", "BC", "ABC")
  edges <- parent_edges(c(edge("A", "AC"),
                          edge("B", "BC"),
                          admixture_edge("C", "AC", "BC"),
                          edge("AC", "ABC"),
                          edge("BC", "ABC"),
                          edge("D", "ABC")))
  admix <- admixture_proportions(c(admix_props("C", "AC", "BC", "a")))
  graph <- agraph(leaves, inner_nodes, edges, admix)
  
  AB_paths <- all_paths(graph, "A", "B")
  CB_paths <- all_paths(graph, "C", "B")
  CD_paths <- all_paths(graph, "C", "D")
    
  expect_equal(all_path_overlaps(AB_paths, AB_paths), f2(graph, "A", "B"))
  expect_equal(all_path_overlaps(CB_paths, CD_paths), f3(graph, "C", "B", "D"))
  expect_equal(all_path_overlaps(AB_paths, CD_paths), f4(graph, "A", "B", "C", "D"))
})
  
test_that("we get the right symbolic representation of f statistics.", {
  leaves <- c("A", "B", "C", "D")
  inner_nodes <- c("AC", "BC", "ABC")
  edges <- parent_edges(c(edge("A", "AC"),
                          edge("B", "BC"),
                          admixture_edge("C", "AC", "BC"),
                          edge("AC", "ABC"),
                          edge("BC", "ABC"),
                          edge("D", "ABC")))
  admix <- admixture_proportions(c(admix_props("C", "AC", "BC", "a")))
  graph <- agraph(leaves, inner_nodes, edges, admix)

  expect_equal(toString(sf2(graph, "A", "B")), "edge_AC_A + edge_ABC_BC + edge_ABC_AC + edge_BC_B")
  expect_equal(toString(sf3(graph, "A", "B", "D")), "edge_AC_A + edge_ABC_AC")
  expect_equal(toString(sf3(graph, "C", "B", "D")),
  "a * a * (edge_ABC_AC + edge_AC_C) + a * (1 - a) * (-edge_ABC_BC) + (1 - a) * (1 - a) * (edge_BC_C)")
  expect_equal(toString(sf4(graph, "A", "B", "C", "D")),
               "a * (edge_ABC_AC) + (1 - a) * (-edge_ABC_BC)")
  expect_equal(toString(sf4(graph, "B", "D", "A", "C")), "(1 - a) * (-edge_ABC_BC)")
  expect_equal(toString(sf4(graph, "B", "D", "C", "A")), "(1 - a) * (edge_ABC_BC)")
})