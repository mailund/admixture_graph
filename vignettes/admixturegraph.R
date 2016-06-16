## ----preamble, echo = FALSE, warning=FALSE, message=FALSE----------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

library(admixturegraph)
library(neldermead)

## ---- echo=FALSE---------------------------------------------------------
leaves <- c("A", "B", "C")
inner_nodes <- c("AB", "ABC")
edges <- parent_edges(c(edge("A", "AB"),
                        edge("B", "AB"),
                        edge("AB", "ABC"),
                        edge("C", "ABC")))
graph <- agraph(leaves, inner_nodes, edges, NULL)
plot(graph)

## ------------------------------------------------------------------------
leaves <- c("A", "B", "C")
inner_nodes <- c("AB", "ABC")
edges <- parent_edges(c(edge("A", "AB"),
                        edge("B", "AB"),
                        edge("AB", "ABC"),
                        edge("C", "ABC")))
graph <- agraph(leaves, inner_nodes, edges, NULL)

## ------------------------------------------------------------------------
plot(graph)  

## ------------------------------------------------------------------------
plot(graph, show_inner_node_labels = TRUE)  

## ------------------------------------------------------------------------
leaves <- c("A", "B", "C")
inner_nodes <- c("a", "c", "ABC")
edges <- parent_edges(c(edge("A", "a"), edge("a", "ABC"),
                        edge("C", "c"), edge("c", "ABC"),
                        admixture_edge("B", "a", "c")))
graph <- agraph(leaves, inner_nodes, edges, NULL)

plot(graph, show_inner_node_labels = TRUE)

## ------------------------------------------------------------------------
leaves <- c("A", "B", "C")
inner_nodes <- c("a", "c", "b", "ABC")
edges <- parent_edges(c(edge("A", "a"), edge("a", "ABC"),
                        edge("C", "c"), edge("c", "ABC"),
                        edge("B", "b"),
                        admixture_edge("b", "a", "c")))
graph <- agraph(leaves, inner_nodes, edges, NULL)

plot(graph, show_inner_node_labels = TRUE)

## ------------------------------------------------------------------------
leaves <- c("A", "B", "C")
inner_nodes <- c("a", "c", "b", "ABC")
edges <- parent_edges(c(edge("A", "a"), edge("a", "ABC"),
                        edge("C", "c"), edge("c", "ABC"),
                        edge("B", "b"),
                        admixture_edge("b", "a", "c")))
admixtures <- admixture_proportions(c(
    admix_props("b", "a", "c", "x")
    ))
graph <- agraph(leaves, inner_nodes, edges, admixtures)

plot(graph, show_inner_node_labels = TRUE, show_admixture_labels = TRUE)

## ----bears_graph, fig.width=6, fig.height=5, cache=TRUE------------------
leaves <- c("BLK", "PB", "Adm1", "Adm2", "Bar", "Chi1", "Chi2", "Denali", "Kenai", "Sweden") 
inner_nodes <- c("R", "r", "s", "t", "u", "v", "w", "x", "y", "z", "M")
edges <- parent_edges(c(edge("BLK", "R"),
                        edge("PB", "s"),
                        edge("Adm1", "x"),
                        edge("Adm2", "x"),
                        edge("Bar", "y"),
                        edge("Chi1", "z"),
                        edge("Chi2", "z"),
                        edge("Denali", "v"),
                        edge("Kenai", "u"),
                        edge("Sweden", "t"),
                        edge("r", "R"),
                        edge("s", "r"),
                        edge("t", "r"),
                        edge("u", "t"),
                        edge("v", "u"),
                        edge("w", "M"),
                        edge("x", "w"),
                        edge("y", "w"),
                        edge("z", "y"),
                        admixture_edge("M", "s", "v")))
admixtures <- admixture_proportions(c(admix_props("M", "s", "v", "a")))
bears_graph <- agraph(leaves, inner_nodes, edges, admixtures)

plot(bears_graph, platform = 1.4, show_admixture_labels = TRUE)

## ------------------------------------------------------------------------
sf2(bears_graph, "Bar", "Chi1")
sf3(bears_graph, "Bar", "Chi1", "Chi2")
sf4(bears_graph, "BLK", "Chi1", "Bar", "Chi2")

## ----bears_data----------------------------------------------------------
data(bears)
bears

## ---- fig.height=5, fig.width=6------------------------------------------
plot(f4stats(bears))

## ----fitting_data--------------------------------------------------------
bears_fit <- fit_graph(bears, bears_graph)

## ------------------------------------------------------------------------
summary(bears_fit)

## ---- fig.height=5, fig.width=6------------------------------------------
plot(bears_fit)

## ------------------------------------------------------------------------
length(four_leaves_graphs)

length(five_leaves_graphs)

length(six_leaves_graphs)

length(seven_leaves_trees)

length(eight_leaves_trees)

example_graph <- four_leaves_graphs[[15]](c("A", "B", "C", "D"))
plot(example_graph, draw_inner_nodes = FALSE)

## ----new_leaves, fig.height=15, fig.width=6------------------------------
example_list_1 <- add_a_leaf(example_graph, "E")

original <- graphics::par()$mfrow
graphics::par(mfrow = c(6, 2))
for (graph in example_list_1) {
  plot(graph, draw_inner_nodes = FALSE)
}

graphics::par(mfrow = original)

## ----new_admixtures------------------------------------------------------
example_list_2 <- add_an_admixture(example_graph, "b")
length(example_list_2)
example_list_3 <- add_an_admixture2(example_graph, "b")
length(example_list_3)

plot(example_list_2[[5]], draw_inner_nodes = FALSE)

## ----new_root------------------------------------------------------------
graph_with_a_new_root <- make_an_outgroup(example_list_2[[5]], "D")
plot(graph_with_a_new_root, draw_inner_nodes = FALSE)

