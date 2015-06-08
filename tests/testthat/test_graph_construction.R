context("Graph construction")

test_that("we can build a simple tree.", {
  nodes <- c("A", "B", "C", "AB", "ABC")
  edges <- matrix(ncol = 2, byrow=TRUE,
                  data = c("A", "AB",
                           "B", "AB",
                           "AB", "ABC",
                           "C", "ABC"))
  admixture_proportions <- NULL
  
  graph <- agraph(nodes, edges, admixture_proportions)
  
  expect_equal(graph$nodes, nodes)
  
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