### flow polar -> brown
nodes <- c("BLK", "PB", "AK",
           "ABC_BC", "ABC_A", "YB", "BB", "EBB",
           "R", "a", "b", "c", "d", "e", "f", "g", "h", "i", "j",
           "abc_bc", "I", "E", "bb")

edges <- matrix(ncol = 2, byrow=TRUE,
                data = c("BLK", "R",
                         "PB", "j", "AK", "j",
                         "ABC_BC", "abc_bc", "abc_bc", "h", "abc_bc", "i",
                         "ABC_A", "i", 
                         "YB", "f",
                         "BB", "bb", "bb", "e", "bb", "g",
                         "EBB", "g",
                         "j", "h",
                         "h", "d",
                         "i", "I", "I", "d", "I", "f",
                         "d", "b",
                         "f", "e",
                         "e", "E", "E", "b", "E", "c",
                         "b", "a",
                         "g", "c",
                         "c", "a",
                         "a", "R"
                ))

admixture_proportions <- matrix(ncol = 3, byrow=TRUE,
                                data = c(
                                  "abc_bc", "h", "a", "abc_bc", "i", "(1-a)",
                                  "I", "d", "b", "I", "f", "(1-b)",
                                  "bb", "e", "c", "bb", "g", "(1-c)",
                                  "E", "b", "d", "E", "c", "(1-d)"
                                ))


graph <- agraph(nodes, edges, admixture_proportions)
plot(graph, ordered_leaves = c("BLK", "PB", "AK", "ABC_BC", "ABC_A", "YB", "BB", "EBB"), 
     show_admixture_labels = TRUE, show_inner_node_labels = TRUE)

(xx <- sf4(graph, "BLK", "PB", "AK", "BB"))
(yy <- sf4(graph, "PB", "BLK", "BB", "AK"))

env <- graph_environment(extract_graph_parameters(graph))
eval(xx, env)
eval(yy, env)

data <- read.table('testdata.txt', header=TRUE)
data %>% add_graph_f4_sign(graph) %>% add_graph_f4(graph, env)
