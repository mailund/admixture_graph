context("Graph construction")

test_that("we can build a simple tree.", {
  leaves <- c("A", "B", "C")
  inner_nodes <- c("AB", "ABC")
  edges <- parent_edges(c(edge("A", "AB"),
                          edge("B", "AB"),
                          edge("AB", "ABC"),
                          edge("C", "ABC")))
  graph <- agraph(leaves, inner_nodes, edges, NULL)
  
  expect_equal(graph$leaves, leaves)
  expect_equal(graph$inner_nodes, inner_nodes)
  expect_equal(graph$nodes, c(leaves, inner_nodes))
  
  expect_equal(names(which(graph$children["A",])), character(0))
  expect_equal(names(which(graph$children["B",])), character(0))
  expect_equal(names(which(graph$children["C",])), character(0))
  expect_equal(names(which(graph$children["AB",])), c("A", "B"))
  expect_equal(names(which(graph$children["ABC",])), c("C", "AB"))
  
  expect_equal(names(which(graph$parents["A",])), c("AB"))
  expect_equal(names(which(graph$parents["B",])), c("AB"))
  expect_equal(names(which(graph$parents["C",])), c("ABC"))
  expect_equal(names(which(graph$parents["AB",])), c("ABC"))
  expect_equal(names(which(graph$parents["ABC",])), character(0))
})

test_that("we can build an unresolved simple tree.", {
  leaves <- c("A", "B", "C")
  inner_nodes <- c("ABC")
  edges <- parent_edges(c(edge("A", "ABC"),
                          edge("B", "ABC"),
                          edge("C", "ABC")))
  graph <- agraph(leaves, inner_nodes, edges, NULL)
  
  expect_equal(graph$leaves, leaves)
  expect_equal(graph$inner_nodes, inner_nodes)
  expect_equal(graph$nodes, c(leaves, inner_nodes))
  
  expect_equal(names(which(graph$children["A",])), character(0))
  expect_equal(names(which(graph$children["B",])), character(0))
  expect_equal(names(which(graph$children["C",])), character(0))
  expect_equal(names(which(graph$children["ABC",])), c("A", "B", "C"))
  
  expect_equal(names(which(graph$parents["A",])), c("ABC"))
  expect_equal(names(which(graph$parents["B",])), c("ABC"))
  expect_equal(names(which(graph$parents["C",])), c("ABC"))
  expect_equal(names(which(graph$parents["ABC",])), character(0))
})

test_that("we can build a simple graph.", {
  leaves <- c("A", "B", "C")
  inner_nodes <- c("AC", "BC", "ABC")
  edges <- parent_edges(c(edge("A", "AC"),
                          edge("B", "BC"),
                          admixture_edge("C", "AC", "BC"),
                          edge("AC", "ABC"),
                          edge("BC", "ABC")))
  admixtures <- admixture_proportions(c(admix_props("C", "AC", "BC", "a")))
  
  graph <- agraph(leaves, inner_nodes, edges, admixtures)
  
  expect_equal(graph$leaves, leaves)
  expect_equal(graph$inner_nodes, inner_nodes)
  expect_equal(graph$nodes, c(leaves, inner_nodes))
  
  expect_equal(names(which(graph$children["A",])), character(0))
  expect_equal(names(which(graph$children["B",])), character(0))
  expect_equal(names(which(graph$children["C",])), character(0))
  expect_equal(names(which(graph$children["AC",])), c("A", "C"))
  expect_equal(names(which(graph$children["BC",])), c("B", "C"))
  expect_equal(names(which(graph$children["ABC",])), c("AC", "BC"))
  
  expect_equal(names(which(graph$parents["A",])), c("AC"))
  expect_equal(names(which(graph$parents["B",])), c("BC"))
  expect_equal(names(which(graph$parents["C",])), c("AC", "BC"))
  expect_equal(names(which(graph$parents["AC",])), c("ABC"))
  expect_equal(names(which(graph$parents["BC",])), c("ABC"))
  expect_equal(names(which(graph$parents["ABC",])), character(0))
  
  expect_equal(graph$probs["C","AC"], c("a"))
  expect_equal(graph$probs["C","BC"], c("(1 - a)"))
  expect_equal(graph$probs["AC","C"], c("a"))
  expect_equal(graph$probs["BC","C"], c("(1 - a)"))
  
})

test_that("we get errors when specifying an edge between non-existing nodes.", {
  leaves <- c("A", "B", "C")
  inner_nodes <- c("AC", "BC", "ABC")
  edges <- parent_edges(c(edge("A", "AC"),
                          edge("B", "BC"),
                          admixture_edge("C", "AC", "BC"),
                          edge("AC", "x"),
                          edge("BC", "ABC")))
  admixtures <- admixture_proportions(c(admix_props("C", "AC", "BC", "a")))
  
  expect_that(agraph(leaves, inner_nodes, edges, admixtures), throws_error())
})

test_that("we get errors when specifying admixture proportions for non-existing edges.", {
  leaves <- c("A", "B", "C")
  inner_nodes <- c("AC", "BC", "ABC")
  edges <- parent_edges(c(edge("A", "AC"),
                          edge("B", "BC"),
                          admixture_edge("C", "AC", "BC"),
                          edge("AC", "ABC"),
                          edge("BC", "ABC")))
  admixtures <- admixture_proportions(c(admix_props("C", "ABC", "BC", "a")))
  
  expect_that(agraph(leaves, inner_nodes, edges, admixtures), throws_error())
})


