
data(bears)

leaves <- c("BLK", "PB",
            "Bar", "Chi1", "Chi2", "Adm1", "Adm2",
            "Denali", "Kenai", "Sweden") 
inner_nodes <- c("R", "PBBB",
                 "Adm", "Chi", "BC", "ABC",
                 "x", "y", "z",
                 "bc_a1", "pb_a1", "abc_a2", "pb_a2")

edges <- parent_edges(c(edge("BLK", "R"),
                        edge("PB", "pb_a1"),
                        edge("pb_a1", "pb_a2"),
                        edge("pb_a2", "PBBB"),
                        
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
                        
                        edge("x", "y"),
                        edge("Kenai", "y"),
                        
                        edge("y", "z"),
                        edge("Sweden", "z"),
                        
                        edge("z", "PBBB"),
                        edge("PBBB", "R")))


admixtures <- admixture_proportions(c(admix_props("bc_a1", "pb_a1", "ABC", "a"),
                                      admix_props("abc_a2", "pb_a2", "x", "b")))

bears_graph <- agraph(leaves, inner_nodes, edges, admixtures)
plot(bears_graph, show_inner_node_labels = TRUE, show_admixture_labels = TRUE)



unpack_environment <- function(parameters, x) {
  n_edges <- length(parameters$edges)
  n_admix_prop <- length(parameters$admix_prop)
  edges <- x[1:n_edges]
  admix_prop <- x[ (n_edges + 1) : (n_edges + n_admix_prop)]
  graph_environment(parameters, edges, admix_prop)
}


make_cost_function <- function(data, graph,
                               parameters = extract_graph_parameters(graph)) {
  force(data)
  force(graph)
  force(parameters)
  
  goal <- data$D
  expressions <- Map(function(W,X,Y,Z) sf4(graph, W, X, Y, Z), data$W, data$X, data$Y, data$Z)
  
  function(x) {
    env <- unpack_environment(parameters, x)
    predictions <- unlist(Map(function(expression) eval(expression, env), expressions), use.names = FALSE)
    sum( (goal - predictions) ** 2)
  }
}



(params <- extract_graph_parameters(bears_graph))
(init_env <- graph_environment(params))

make_fun <- function(fit, variable) {
  
  force(variable)
  
  data <- fitted(fit)
  goal <- data$D
  expressions <- Map(function(W,X,Y,Z) sf4(fit$graph, W, X, Y, Z), data$W, data$X, data$Y, data$Z)
  
  function(value) {
    
    cfun <- function(x) {
      env <- unpack_environment(parameters, x)
      assign(variable, value, env)
      predictions <- unlist(Map(function(expression) eval(expression, env), expressions), use.names = FALSE)
      sum( (goal - predictions) ** 2)
    }  
    
    params <- extract_graph_parameters(fit$graph)
    init_env <- graph_environment(params)
    x0 <- pack_environment(params, init_env)
    cfunc <- make_cost_function(data, fit$graph, params)
    
    opti <- neldermead::fminbnd(cfunc, x0 = x0, xmin = rep(0, length(x0)),
                                xmax = rep(1, length(x0)),
                                options = optimisation_options)
    
    opti
  }
}
  