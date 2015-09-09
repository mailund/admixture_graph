## ----preamble, echo = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

library(admixturegraph)

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
    admix_props("b", "a", "c", "alpha")
    ))
graph <- agraph(leaves, inner_nodes, edges, admixtures)

plot(graph, show_inner_node_labels = TRUE, show_admixture_labels = TRUE)

