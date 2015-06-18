
data(bears)

leaves <- c("BLK", "PB", "AK", "ABC_A", "ABC_BC", "YB", "EBB") 
inner_nodes <- c("R", "a", "b", "c", "d", "e", "f", "g", "h",
                 "abc_a", "G", "E")

edges <- parent_edges(c(edge("BLK", "R"),
                        edge("PB", "h"),
                        edge("AK", "h"),
                        edge("ABC_A", "abc_a"),
                        admixture_edge("abc_a", "f", "g"),
                        edge("ABC_BC", "g"), 
                        edge("YB", "e"),
                        edge("EBB", "c"),
                        edge("h", "f"),
                        edge("f", "d"),
                        edge("g", "G"),
                        admixture_edge("G", "d", "e"),
                        edge("d", "b"),
                        edge("e", "E"),
                        admixture_edge("E", "b", "c"),
                        edge("b", "a"),
                        edge("c", "a"),
                        edge("a", "R")))


admixtures <- admixture_proportions(c(admix_props("abc_a", "f", "g", "alpha"),
                                      admix_props("G", "d", "e", "beta"),
                                      admix_props("E", "b", "c", "gamma")))

bears_graph <- agraph(leaves, inner_nodes, edges, admixtures)
plot(bears_graph, show_inner_node_labels = TRUE, show_admixture_labels = TRUE)

(x <- fit_graph(bears, bears_graph))
