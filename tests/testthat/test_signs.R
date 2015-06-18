context("Signs")

test_that("we extract the right sign for f statistics.", {
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
  
  expect_equal(overlaps_sign(f2(graph, "A", "B")), +1)
  expect_equal(overlaps_sign(f3(graph, "A", "B", "D")), +1)
  expect_that(overlaps_sign(f4(graph, "A", "B", "C", "D")), equals(NA))
  expect_equal(overlaps_sign(f4(graph, "B", "D", "A", "C")), -1)
  expect_equal(overlaps_sign(f4(graph, "B", "D", "C", "A")), +1)
})
