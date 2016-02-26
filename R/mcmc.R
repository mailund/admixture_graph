
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
  if (n_admix > 0) {
    admix_idx <- 1:n_admix
  } else {
    admix_idx <- c()
  }
    
  edges_idx <- (n_admix+1):length(parameter_names)
                                                     
  logL <- log_likelihood(f, concentration, matrix, graph, params)
  
  transform_to_mcmc_space <- function(state) {
    admix <- state[admix_idx]
    edges <- state[edges_idx]
    # this is necessary to avoid -Inf in the calculations
    admix[admix == 0] <- 1e-16
    edges[edges == 0] <- 1e-16
    c(qnorm(admix), log(edges))
  }
  
  transform_to_graph_space <- function(state) {
    admix <- state[admix_idx]
    edges <- state[edges_idx]
    c(pnorm(admix), exp(edges))
  }
  
  log_prior <- function(state) {
    # just a reasonably wide normal dist in log space...
    sum(log(dnorm(state, sd=1)))
  }
  
  log_likelihood <- function(state) {
    graph_space_state <- transform_to_graph_space(state)
    admix <- graph_space_state[admix_idx]
    edges <- graph_space_state[edges_idx]
    edges <- matrix(edges, ncol=1)
    tryCatch(-logL(admix, edges), finally = -Inf)
  }
  
  proposal <- function(state) {
    rnorm(length(state), mean = state, sd = 0.001)
  }
  
  list(log_prior = log_prior, log_likelihood = log_likelihood,
       
       admixture_parameters = admixture_parameters,
       edge_parmeters = edge_parameters, 
       parameter_names = parameter_names,
       
       proposal = proposal,
       transform_to_graph_space = transform_to_graph_space,
       transform_to_mcmc_space = transform_to_mcmc_space)
}

#' Run a Metropolis-Hasting MCMC to sample graph parameters.
#' 
#' The MCMC performs a random walk in transformed parameter space (edge lengths are log transformed
#' and admixture proportions inverse Normal distribution transformed) and from this
#' explores the posterior distribution of graph parameters.
#' 
#' Using the posterior distribution of parameters is one approach to getting parameter estimates
#' and a sense of their variability. Credibility intervals can directly be obtained from sampled
#' parameters; to get confidence intervals from the likelihood maximisation approach requires
#' either estimating the Hessian matrix for the likelihood or a boot-strapping approach
#' to the data.
#' 
#' From sampling the likelihood values for each sample from the posterior we can also
#' compute the model likelihood (the probability of the data when integrating over all the model
#' parameters). This gives us a direct way of comparing graphs since the ratio of likelihoods
#' is the Bayes factor between the models. Comparing models using maximum likelihood estimtates
#' is more problematic since usually graphs are not nested models.
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
  trace <- matrix(nrow = iterations, ncol = length(model$parameter_names) + 3)
  colnames(trace) <- c(model$parameter_names, "prior", "likelihood", "posterior")
  
  current_state <- mcmc_initial
  current_prior = model$log_prior(current_state)
  current_likelihood = model$log_likelihood(current_state)
  current_posterior <- current_prior + current_likelihood
  
  pb <- utils::txtProgressBar(min = 1, max = iterations, style=3)
  for (i in 1:iterations) {
    trace[i,] <- c(model$transform_to_graph_space(current_state), current_prior, current_likelihood, current_posterior)
    
    proposal_state <- model$proposal(current_state)
    proposal_prior = model$log_prior(proposal_state)
    proposal_likelihood = model$log_likelihood(proposal_state)
    proposal_posterior <- proposal_prior + proposal_likelihood
    
    log_accept_prob <- proposal_posterior - current_posterior
    if (log(runif(1)) < log_accept_prob) {
      current_state <- proposal_state
      current_prior <- proposal_prior
      current_likelihood <- proposal_likelihood
      current_posterior <- proposal_posterior
    }
    
    utils::setTxtProgressBar(pb, i)
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

