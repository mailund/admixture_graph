context("Paths")

test_that("we can extract paths from a graph.", {
  
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
  
  AB_paths <- all_paths(graph, "A", "B")
  CB_paths <- all_paths(graph, "C", "B")
  
  expect_equal(length(AB_paths), 1)
  expect_equal(AB_paths[[1]]$from, c("A", "AC", "ABC", "BC"))
  expect_equal(AB_paths[[1]]$to, c("AC", "ABC", "BC", "B"))
  
  expect_equal(length(CB_paths), 2)
  expect_equal(CB_paths[[1]]$from, c("C", "AC", "ABC", "BC"))
  expect_equal(CB_paths[[1]]$to, c("AC", "ABC", "BC", "B"))
  expect_equal(CB_paths[[1]]$prob, c("a", "", "", ""))
  
  expect_equal(CB_paths[[2]]$from, c("C", "BC"))
  expect_equal(CB_paths[[2]]$to, c("BC", "B"))
  expect_equal(CB_paths[[2]]$prob, c("(1-a)", ""))
  
})
