source("preamble.R")
f4 <- f4 %>% filter(D >= 0) %>% data.frame

suppressPackageStartupMessages(library(neldermead, quietly = TRUE))

#library(devtools)
#install_github("mailund/admixture_graph")
suppressPackageStartupMessages(library(admixturegraph, quietly = TRUE))

map_to_pop <- . %>%
  mutate(W_ind = W, X_ind = X, Y_ind = Y, Z_ind = Z) %>%
  mutate(W = population(W_ind),
         X = population(X_ind),
         Y = population(Y_ind),
         Z = population(Z_ind))

map_to_indv <- . %>%
  mutate(W = W_ind, X = X_ind, Y = Y_ind, Z = Z_ind) %>%
  select(-W_ind, -X_ind, -Y_ind, -Z_ind)

filter_on_leaves <- function(x, graph)
  filter(x, W %in% graph$leaves,
         X %in% graph$leaves,
         Y %in% graph$leaves,
         Z %in% graph$leaves)
filter_unique_pops <- . %>%
  filter(W != X, W != Y, W != Z, X != Y, X != Z, Y != Z)

#remap_fit_graph <- function(x, graph, ...) {
#  x%>% map_to_pop %>% filter_unique_pops %>% filter_on_leaves(graph) %>% fit_graph(graph, ...) -> fit
#  fit$fit_data <- fit$fit_data %>% map_to_indv
#  fit
#}

remap_fit_graph <- function(x, graph, ...) {
  x %>% project_to_population(population) %>%
    filter_on_leaves(graph) %>%
    fit_graph(graph, ...) %>%
    split_population
}

project_to_population_fit <- function(fit, f) {
  fit$data <- project_to_population(fit$data, f)
  fit
}


strict_populations <- c("AK", "PB", "EBB", "YB", "BB", "ABC_A", "ABC_BC")
populations <- c("APB", strict_populations)

f4_pop <- f4 %>%
  filter(X %in% populations, Y %in% populations, Z %in% populations) %>%
  distinct

f4_indv <- f4 %>%
  filter(!(X %in% strict_populations),
         !(Y %in% strict_populations),
         !(Z %in% strict_populations)) %>%
  distinct




leaves <- c("PANDA", "BLK", "APB", "PB", "AK", "ABC_BC", "ABC_A", "YB", "BB",  "EBB")
inner_nodes <- c("RR", "ghost", "R",
                 "blk", "x", "y", "z", "w", "v", "vv",
                 "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L",
                 "abc_a", "xx", "AA", "blk2")

edges <- parent_edges(c(edge("PANDA", "RR"),

                        edge("BLK", "blk"), edge("blk", "blk2"), edge("blk2", "R"),
                        edge("R", "RR"),

                        edge("APB", "x"),
                        admixture_edge("x", "y", "C"),
                        admixture_edge("y", "blk", "z"),
                        edge("z", "A"),
                        edge("A", "AA"),
                        admixture_edge("AA", "blk2", "ghost"),
                        edge("ghost", "R"),

                        edge("PB", "v"),
                        edge("AK", "vv"), edge("vv", "v"),
                        edge("v", "w"),
                        admixture_edge("w", "z", "G"),

                        edge("ABC_BC", "F"),
                        edge("F", "G"),
                        edge("G", "D"),
                        edge("D", "C"),
                        edge("C", "xx"),   edge("xx", "B"),
                        edge("B", "A"),

                        edge("ABC_A", "abc_a"),
                        admixture_edge("abc_a", "F", "xx"),

                        edge("YB", "E"),
                        edge("E", "D"),

                        edge("BB", "L"),
                        admixture_edge("L", "vv", "K"),
                        admixture_edge("K", "I", "H"),
                        admixture_edge("I", "E", "J"),
                        edge("H", "ghost"),
                        edge("J", "B"),

                        edge("EBB", "J")))

admix_weights <- admixture_proportions(c(admix_props("y", "blk", "z", "a"),
                                         admix_props("x", "y", "C", "b"),
                                         admix_props("w", "z", "G", "c"),
                                         admix_props("I", "E", "J", "d"),
                                         admix_props("K", "I", "H", "e"),
                                         admix_props("L", "vv", "K", "f"),
                                         admix_props("abc_a", "F", "xx", "g"),
                                         admix_props("AA", "blk2", "ghost", "h")))

alt_graph_26 <- agraph(leaves, inner_nodes, edges, admix_weights)
#pdf("alg-graph-26.pdf")
plot(alt_graph_26, show_admixture_labels = TRUE, show_inner_node_labels = TRUE, main = 'Alternative Graph 26')
#dev.off()

#f4_indv %>% remap_fit_graph(alt_graph_26) -> alt_graph_26_indv_fit
#f4_pop %>% remap_fit_graph(alt_graph_26) -> alt_graph_26_pop_fit


leaves <- c("BLK", "APB", "PB", "AK", "ABC_BC", "ABC_A", "YB", "BB", "EBB")
inner_nodes <- c("R", "a", "b", "c", "d", "e", "f", "g", "h", "i", "j",
                 "abc_bc", "I", "E", "bb", "x")

edges <- parent_edges(c(edge("BLK", "R"),

                        edge("APB", "x"),
                        edge("PB", "j"),
                        edge("AK", "j"),

                        edge("ABC_BC", "abc_bc"),
                        admixture_edge("abc_bc", "h", "i"),
                        edge("ABC_A", "i"),
                        edge("YB", "f"),
                        edge("BB", "bb"),
                        admixture_edge("bb", "e", "g"),
                        edge("EBB", "g"),
                        edge("j", "h"),
                        edge("h", "d"),
                        edge("i", "I"),
                        admixture_edge("I", "d", "f"),
                        edge("d", "b"),
                        edge("f", "e"),
                        edge("e", "E"),
                        admixture_edge("E", "b", "c"),
                        edge("b", "x"),
                        edge("x", "a"),
                        edge("g", "c"),
                        edge("c", "a"),
                        edge("a", "R")))

admix_weights <- admixture_proportions(c(
  admix_props("abc_bc", "h", "i", "a"),
  admix_props("I", "d", "f", "b"),
  admix_props("bb", "e", "g", "c"),
  admix_props("E", "b", "c", "d")
))


graph_PB_to_BB_APB_3 <- agraph(leaves, inner_nodes, edges, admix_weights)
pdf("Fig3-B1.pdf")
plot(graph_PB_to_BB_APB_3, show_admixture_labels = TRUE)
dev.off()

leaves <- c("BLK", "APB", "PB", "AK", "ABC_BC", "ABC_A", "YB", "BB", "EBB")
inner_nodes <- c("pbr", "pba", "bb", "R", "a", "b", "c", "d", "e", "f", "g", "h")

edges <- parent_edges(c(edge("BLK", "R"),
                        edge("APB", "b"),
                        edge("PB", "pbr"),
                        edge("AK", "pbr"),
                        edge("pbr", "pba"),
                        admixture_edge("pba", "b", "f"),
                        edge("ABC_BC", "f"),
                        edge("f", "e"),
                        edge("ABC_A", "e"),
                        edge("YB", "g"),
                        edge("BB", "bb"),
                        admixture_edge("bb", "g", "h"),
                        edge("EBB", "h"),
                        edge("b", "a"),
                        edge("e", "d"),
                        edge("g", "d"),
                        edge("d", "c"),
                        edge("h", "c"),
                        edge("c", "a"),
                        edge("a", "R")))

admix_weights <- admixture_proportions(c(admix_props("pba", "b", "f", "a"),
                                         admix_props("bb", "g", "h", "b")))

graph_BB_to_PB <- agraph(leaves, inner_nodes, edges, admix_weights)
pdf("Fig3-B2.pdf")
plot(graph_BB_to_PB, show_admixture_labels = TRUE)
dev.off()


## Fitting
fit_BB_to_PB_pop <- f4_pop %>% remap_fit_graph(graph_BB_to_PB)
fit_BB_to_PB_indv <- f4_indv %>% remap_fit_graph(graph_BB_to_PB)

fit_PB_to_BB_APB_3_pop <- f4_pop %>% remap_fit_graph(graph_PB_to_BB_APB_3)
fit_PB_to_BB_APB_3_indv <- f4_indv %>% remap_fit_graph(graph_PB_to_BB_APB_3)

