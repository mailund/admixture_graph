context("Path overlaps")

test_that("we can compute the overlap between paths.", {
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
  expect_equal(length(CB_paths), 2)
  
  overlaps1 <- path_overlap(AB_paths[[1]], CB_paths[[1]])
  overlaps2 <- path_overlap(AB_paths[[1]], CB_paths[[2]])
  
  expect_equal(overlaps1$prob, c("a"))
  expect_equal(nrow(overlaps1$positive), 3)
  expect_equal(overlaps1$positive$from, c("ABC", "AC",  "BC"))
  expect_equal(overlaps1$positive$to,   c("BC",  "ABC", "B"))
  expect_equal(nrow(overlaps1$negative), 0)
  
  expect_equal(overlaps2$prob, c("(1-a)"))
  expect_equal(nrow(overlaps2$positive), 1)
  expect_equal(overlaps2$positive$from, c("BC"))
  expect_equal(overlaps2$positive$to,   c("B"))
  expect_equal(nrow(overlaps2$negative), 0)
  
  all_overlaps <- all_path_overlaps(AB_paths, CB_paths)
  
  expect_equal(all_overlaps[[1]], overlaps1)
  expect_equal(all_overlaps[[2]], overlaps2)
})
  
  