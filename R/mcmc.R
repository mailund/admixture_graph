
#' Collect the information about a graph and a data set needed to run an MCMC on it.
#' 
#' @param graph  The admixture graph to analyse
#' @param data   The data set to compute the posterior over
#' @return A model object wrapping functions and data needed to sample from the MCMC.
#' @export
make_mcmc_model <- function(graph, data) {
  
  f <- data$D
  concentration <- calculate_concentration(data, Z.value = TRUE) # FIXME: use the empirical covariance matrix
  params <- extract_graph_parameters(graph)
  matrix <- build_edge_optimisation_matrix(data, graph, params)$column_reduced
  
  admixture_parameters = params$admix_prop
  edge_parameters <- colnames(matrix)
  parameter_names <- c(admixture_parameters, edge_parameters)
    
  n_admix <- length(admixture_parameters)
  n_edges <- length(edge_parameters)
  admix_idx <- 1:n_admix
  edges_idx <- (n_admix+1):length(parameter_names)
                                                     
  logL <- log_likelihood(f, concentration, matrix, graph, params)
  
  transform_to_mcmc_space <- function(state) {
    admix <- state[admix_idx]
    edges <- state[edges_idx]
    c(qnorm(admix), log(edges))
  }
  
  transform_to_graph_space <- function(state) {
    admix <- state[admix_idx]
    edges <- state[edges_idx]
    c(pnorm(admix), exp(edges))
  }
  
  log_prior <- function(state) {
    return(sum(log(dnorm(state)))) # FIXME
    
    # here I'm just using a uniform prior on [0,1] for all parameters. Can always change that later.
    if (any(state < 0 || state > 1)) {
      prior <- 0
    } else {
      prior <- 1
    }
    log(prior)
  }
  
  log_likelihood <- function(state) {
    graph_space_state <- transform_to_graph_space(state)
    admix <- graph_space_state[admix_idx]
    edges <- graph_space_state[edges_idx]
    edges <- matrix(edges, ncol=1)
    tryCatch(logL(admix, edges), finally = -Inf)
  }
  
  log_posterior <- function(state) {
    log_prior(state) + log_likelihood(state)
  }

  proposal <- function(state) {
    rnorm(length(state), mean = state)
  }
  
  list(log_posterior = log_posterior, 
       admixture_parameters = admixture_parameters,
       edge_parmeters = edge_parameters, 
       parameter_names = parameter_names,
       proposal = proposal,
       transform_to_graph_space = transform_to_graph_space,
       transform_to_mcmc_space = transform_to_mcmc_space)
}

#' Run a Metropolis-Hasting MCMC to sample graph parameters.
#' 
#' 
#' 
#' @param model          Object constructed with \code{\link{make_mcmc_model}}.
#' @param initial_state  The initial set of graph parameters.
#' @param iterations     Number of iterations to sample.
#' @return A matrix containing the trace of the MCMC run.
#' @export
run_metropolis_hasting <- function(model, initial_state, iterations) {
  
  if (length(initial_state) != length(model$parameter_names)) {
    stop(paste0("The length of the initial state, (", length(initial_state), 
                ") does not match the number of graph parameters (",
                length(model$parameter_names), ")."))
  }
  
  mcmc_initial <- model$transform_to_mcmc_space(initial_state)
  trace <- matrix(nrow = iterations, ncol = length(model$parameter_names) + 1)
  colnames(trace) <- c(model$parameter_names, "posterior")
  
  current_state <- mcmc_initial
  current_posterior <- model$log_posterior(current_state)
  
  for (i in 1:iterations) {
    trace[i,] <- c(model$transform_to_graph_space(current_state), current_posterior)
    proposal_state <- model$proposal(current_state)
    proposal_posterior <- model$log_posterior(proposal_state)
    log_accept_prob <- proposal_posterior - current_posterior
    if (log(runif(1)) < log_accept_prob) {
      current_state <- proposal_state
      current_posterior <- proposal_posterior
    }
  }
  
  trace 
}

#' Removes the first k rows from a trace.
#' 
#' @param trace  A trace from an MCMC run.
#' @param k      Number of rows to discard as burn-in.
#' @export
burn_in <- function(trace, k) {
  trace[-seq(1,k), ]
}

#' Thins out an MCMC trace.
#' 
#' @param trace  A trace from an MCMC run.
#' @param k      The number of lines to skip over per retained sample.
#' @export
thinning <- function(trace, k) {
  trace[seq(1,nrow(trace), by=k),]
}

