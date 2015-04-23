


nodes <- c("BLK", "APB", "PB", "pb", "BC", "A", "YB", "BB", "bb", "EBB",
           "R", "a", "b", "c", "d", "e", "f", "g", "h")

edges <- matrix(ncol = 2, byrow=TRUE,
                data = c("BLK", "R",
                         "APB", "b",
                         "PB", "pb", "pb", "b", "pb", "f",
                         "BC", "f",
                         "f", "e",
                         "A", "e",
                         "YB", "g",
                         "BB", "bb", "bb", "g", "bb", "h",
                         "EBB", "h",
                         "b", "a",
                         "e", "d",
                         "g", "d",
                         "d", "c",
                         "h", "c",
                         "c", "a",
                         "a", "R"
                         ))

admixture_proportions <- matrix(ncol = 3, byrow=TRUE,
                                data = c("pb", "b", "alpha", "pb", "f", "(1-alpha)",
                                         "bb", "g", "beta",  "bb", "h", "(1-beta)"))


graph <- agraph(nodes, edges, admixture_proportions)
f4(graph, "BLK", "APB", "BB", "EBB")
