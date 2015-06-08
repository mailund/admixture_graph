context("Graph construction")

test_that("we can build a simple tree.", {
  leaves <- c("A", "B", "C")
  inner_nodes <- c("AB", "ABC")
  edges <- matrix(ncol = 2, byrow=TRUE,
                  data = c("A", "AB",
                           "B", "AB",
                           "AB", "ABC",
                           "C", "ABC"))
  graph <- agraph(leaves, inner_nodes, edges, NULL)
  
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
  edges <- matrix(ncol = 2, byrow=TRUE,
                  data = c("A", "ABC",
                           "B", "ABC",
                           "C", "ABC"))
  graph <- agraph(leaves, inner_nodes, edges, NULL)
  
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
  edges <- matrix(ncol = 2, byrow=TRUE,
                  data = c("A", "AC",
                           "B", "BC",
                           "C", "AC", "C", "BC",
                           "AC", "ABC",
                           "BC", "ABC"))
  admixture_proportions <- matrix(ncol = 3, byrow=TRUE,
                                  data = c("C", "AC", "a", "C", "BC", "(1-a)"))
  
  graph <- agraph(leaves, inner_nodes, edges, admixture_proportions)
  
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
  expect_equal(graph$probs["C","BC"], c("(1-a)"))
  expect_equal(graph$probs["AC","C"], c("a"))
  expect_equal(graph$probs["BC","C"], c("(1-a)"))
  
})
