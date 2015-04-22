


nodes <- c("A", "B", "C", "AB", "BC", "ABC", "R", "O")
edges <- matrix(ncol = 2, byrow=TRUE,
                data = c("A", "AB",
                         "C", "BC", 
                         "B", "AB", 
                         "B", "BC", 
                         "AB", "ABC",
                         "BC", "ABC",
                         "ABC", "R",
                         "O", "R"))
admixture_proportions <- matrix(ncol = 3, byrow=TRUE,
                                data = c("B", "AB", "a", "B", "BC", "(1-a)"))


graph <- agraph(nodes, edges, admixture_proportions)
f4(graph, "O", "A", "B", "C")
f3(graph, "O", "A", "C")
f3(graph, "O", "A", "B")
