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
    admix_props("b", "a", "c", "alpha")
    ))
graph <- agraph(leaves, inner_nodes, edges, admixtures)

plot(graph, show_inner_node_labels = TRUE, show_admixture_labels = TRUE)

## ----bears_graph, fig.width=6, fig.height=5, cache=TRUE------------------
leaves <- c("BLK", "PB",
            "Bar", "Chi1", "Chi2", "Adm1", "Adm2",
            "Denali", "Kenai", "Sweden") 
inner_nodes <- c("R", "PBBB",
                 "Adm", "Chi", "BC", "ABC",
                 "x", "y", "z",
                 "pb_a1", "pb_a2", "pb_a3", "pb_a4",
                 "bc_a1", "abc_a2", "x_a3", "y_a4")

edges <- parent_edges(c(edge("BLK", "R"),
                        edge("PB", "pb_a1"),
                        edge("pb_a1", "pb_a2"),
                        edge("pb_a2", "pb_a3"),
                        edge("pb_a3", "pb_a4"),
                        edge("pb_a4", "PBBB"),
                        
                        edge("Chi1", "Chi"),
                        edge("Chi2", "Chi"),
                        edge("Chi", "BC"),
                        edge("Bar", "BC"),
                        edge("BC", "bc_a1"),
                        
                        edge("Adm1", "Adm"),
                        edge("Adm2", "Adm"),
                        
                        admixture_edge("bc_a1", "pb_a1", "ABC"),
                        edge("Adm", "ABC"),
                        
                        edge("ABC", "abc_a2"),
                        admixture_edge("abc_a2", "pb_a2", "x"),
                        
                        edge("Denali", "x"),
                        edge("x", "x_a3"),
                        admixture_edge("x_a3", "pb_a3", "y"),
                      
                        edge("Kenai", "y"),
                        edge("y", "y_a4"),                        
                        admixture_edge("y_a4", "pb_a4", "z"),
                        
                        edge("Sweden", "z"),
                        
                        edge("z", "PBBB"),
                        edge("PBBB", "R")))
 

admixtures <- admixture_proportions(c(admix_props("bc_a1", "pb_a1", "ABC", "a"),
                                      admix_props("abc_a2", "pb_a2", "x", "b"),
                                      admix_props("x_a3", "pb_a3", "y", "c"),
                                      admix_props("y_a4", "pb_a4", "z", "d")))
                                
bears_graph <- agraph(leaves, inner_nodes, edges, admixtures)
plot(bears_graph, show_admixture_labels = TRUE)

## ------------------------------------------------------------------------
sf4(bears_graph, "BLK", "PB", "Adm1", "Adm2")
sf2(bears_graph, "Bar", "Chi1")
sf3(bears_graph, "Bar", "Chi1", "Chi2")
sf4(bears_graph, "BLK", "PB", "Bar", "Adm2")

## ------------------------------------------------------------------------
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
coef(bears_fit)
fitted(bears_fit)

