


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

pdf('bear-graph.pdf', width = 6)
plot(graph, c("BLK", "APB", "PB", "BC", "A", "YB", "BB", "EBB"))
dev.off()

f4(graph, "BLK", "PB", "BC", "EBB")
f4(graph, "BLK", "PB", "A", "EBB")
f4(graph, "BLK", "PB", "YB", "EBB")
f4(graph, "BLK", "PB", "BB", "EBB")

f4(graph, "BLK", "PB", "A", "BB")
f4(graph, "BLK", "PB", "BC", "BB")
f4(graph, "BLK", "PB", "YB", "BB")

f4(graph, "BLK", "PB", "A", "YB")
f4(graph, "BLK", "PB", "BC", "YB")

f4(graph, "BLK", "PB", "BC", "A")


