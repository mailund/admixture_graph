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
plot(graph)  

