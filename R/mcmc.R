
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
    c(stats::qnorm(admix), log(edges))
  }
  
  transform_to_graph_space <- function(state) {
    admix <- state[admix_idx]
    edges <- state[edges_idx]
    c(stats::pnorm(admix), exp(edges))
  }
  
  log_prior <- function(state) {
    # just a reasonably wide normal dist in log space...
    sum(log(stats::dnorm(state, sd=1)))
  }
  
  log_likelihood <- function(state) {
    graph_space_state <- transform_to_graph_space(state)
    admix <- graph_space_state[admix_idx]
    edges <- graph_space_state[edges_idx]
    edges <- matrix(edges, ncol=1)
    tryCatch(logL(admix, edges), finally = -Inf)
  }
  
  proposal <- function(state) {
    stats::rnorm(length(state), mean = state, sd = 0.001)
  }
  
  multnorm_proposal <- function(state, sigma){
    #     print(paste("state",state))
    #     print(paste("sigma",sigma))
    MASS::mvrnorm(1, mu = state, Sigma = sigma)
  }
  
  list(log_prior = log_prior, log_likelihood = log_likelihood,
       
       admixture_parameters = admixture_parameters,
       edge_parmeters = edge_parameters, 
       parameter_names = parameter_names,
       
       multnorm_proposal=multnorm_proposal,
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
#' @param model            Object constructed with \code{\link{make_mcmc_model}}.
#' @param initial_state    The initial set of graph parameters.
#' @param iterations       Number of iterations to sample.
#' @param no_temperatures  Number of chains in the MC3 procedure
#' @param cores            Number of cores to spread the chains across. Best performance is when cores=no_temperatures
#' @param no_flips         Mean number of times a flip between two chains should be proposed after each step
#' @param max_tmp          The highest temperature
#' 
#' @return A matrix containing the trace of the chain with temperature 1.
#' 
#' @export
run_metropolis_hasting <- function(model, initial_state, iterations, 
                                   no_temperatures = 1, cores = 1, no_flips = 1, max_tmp = 100) {
  
  if (!requireNamespace("parallel", quietly = TRUE)) {
    stop("The MCMC functionality requires that the parallel packate is installed.")
  }
  
  # tlast=proc.time()
  if(no_temperatures>1){
    temperatures <- max_tmp^(0:(no_temperatures-1)/(no_temperatures-1))
  }
  else{
    temperatures <- 1
  }
  
  
  if (length(initial_state) != length(model$parameter_names)) {
    stop(paste0("The length of the initial state, (", length(initial_state), 
                ") does not match the number of graph parameters (",
                length(model$parameter_names), ")."))
  }
  
  mcmc_initial <- model$transform_to_mcmc_space(initial_state)
  trace <- matrix(nrow = iterations, ncol = length(model$parameter_names) + 3)
  colnames(trace) <- c(model$parameter_names, "prior", "likelihood", "posterior")
  
  current_states <- t(replicate(no_temperatures,mcmc_initial)) #
  current_prior <- model$log_prior(mcmc_initial)
  current_priors <- rep(current_prior,no_temperatures)
  current_likelihood <- model$log_likelihood(mcmc_initial)
  current_likelihoods <- rep(current_likelihood,no_temperatures)
  current_posterior <- current_prior + current_likelihood
  current_posteriors <- rep(current_posterior,no_temperatures)
  
  onestep <- function(l){
    current_state <- l$current_state
    current_posterior <- l$current_posterior
    temperature <- l$temperature
    sigma <- l$sigma
    proposal_state <- model$multnorm_proposal(current_state,sigma=sigma)
    proposal_prior <- model$log_prior(proposal_state)
    proposal_likelihood <- model$log_likelihood(proposal_state)
    proposal_posterior <- proposal_prior + proposal_likelihood
    
    log_accept_prob <- proposal_posterior/temperature-current_posterior/temperature
    if (log(stats::runif(1)) < log_accept_prob) {
      current_state <- proposal_state
      current_prior <- proposal_prior
      current_likelihood <- proposal_likelihood
      current_posterior <- proposal_posterior
    }
    return(list(current_state = current_state,
                current_prior = current_prior,
                current_likelihood = current_likelihood, 
                current_posterior = current_posterior, 
                alpha = min(1, exp(log_accept_prob))))
  }
  
  #initialise adaption parameters
  ad_params <- rep(0.02,no_temperatures)
  sigmas <- replicate(no_temperatures, diag(length(current_states[1,]))*0.1, simplify = F)
  common_mean <- current_states[1,]
  
  
  pb <- utils::txtProgressBar(min = 1, max = iterations, style=3)
  
  #variable used for monitoring
  avg_temp_update_probability <- 0
  
  for (i in 1:iterations) {
    
    trace[i,] <- c(model$transform_to_graph_space(current_states[1,]), current_priors[1], current_likelihoods[1], current_posteriors[1])
    
    #making a list of lists of arguments for each chain
    listOfStepInformation=list()
    for(j in 1:no_temperatures){
      stepOfInformation <- list(current_posterior = current_posteriors[j], current_state = current_states[j,], 
                                temperature = temperatures[j], sigma = ad_params[j]*sigmas[[j]])
      listOfStepInformation[[j]] <- stepOfInformation
    }
    
    #making each step take a step
    reses <- parallel::mclapply(listOfStepInformation, onestep, mc.cores=cores)
    
    #Here adaption takes place and the new step is saved
    gamma <- 0.1/sqrt(i)
    for(j in 1:no_temperatures){
      
      #the adaption parameter is updated
      ad_params[j] <- ad_params[j]*exp(gamma*(reses[[j]]$alpha-0.234))
      
      #the new step is saved
      current_states[j,] <- reses[[j]]$current_state
      current_priors[j] <- reses[[j]]$current_prior
      current_likelihoods[j] <- reses[[j]]$current_likelihood
      current_posteriors[j] <- reses[[j]]$current_posterior
      
      #update of the covariance matrix
      common_mean=common_mean+gamma*(current_states[j,]-common_mean)/no_temperatures
      xbar <- current_states[j,]-common_mean
      sigmas[[j]] <- sigmas[[j]]+gamma*(xbar %o% xbar - sigmas[[j]])
    }
    
    #Here flips between temperatures are proposed
    flips <- stats::rpois(1,lambda=no_flips)
    if (no_temperatures > 1) {
      for (p in numeric(flips)) {
        r <- sample(no_temperatures,2)
        m <- r[1]
        j <- r[2]
        current <- current_posteriors[j]/temperatures[j]+current_posteriors[m]/temperatures[m]
        new <- current_posteriors[m]/temperatures[j]+current_posteriors[j]/temperatures[m]
        a <- exp(new-current)
        if (stats::runif(1) < a){
          tmp_post <- current_posteriors[j]
          tmp_lik <- current_likelihoods[j]
          tmp_pri <- current_priors[j]
          tmp_sta <- current_states[j,]
          current_states[j,] <- current_states[m,]
          current_priors[j] <- current_priors[m]
          current_likelihoods[j] <- current_likelihoods[m]
          current_posteriors[j] <- current_posteriors[m]
          current_states[m,] <- tmp_sta
          current_priors[m] <- tmp_pri
          current_likelihoods[m] <- tmp_lik
          current_posteriors[m] <- tmp_post
        }
        if (p==1) {
          avg_temp_update_probability <- ((i-1)*avg_temp_update_probability+a)/i
        }
      }
      
    }
    
    #     t_flip=t_flip+proc.time()-tlast
    #     tlast=proc.time()
    
    utils::setTxtProgressBar(pb, i)
  }
  
  trace 
}

#' Removes the first k rows from a trace.
#' 
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
#' Thins out an MCMC trace.
#' 
#' @param trace  A trace from an MCMC run.
#' @param k      The number of lines to skip over per retained sample.
#' @export
thinning <- function(trace, k) {
  trace[seq(1,nrow(trace), by=k),]
}



## Model likelihoods ####

#' Computes the log of a sum of numbers all given in log-space.
#' 
#' Given a sequence of numbers \eqn{[\log(x_1), \log(x_2), ..., \log(x_n)]}, computes \eqn{\log(\sum_{i=1}^n x_i)}.
#' For adding two numbers that are given in log space we use the expression \code{max(x, y) + log1p(exp( -abs(x - y) ))}
#' which is a good approximation if \code{x} and \code{y} are of the same order of magnitude, but if they
#' are of very different sizes just returns the maximum of the two. To prevent adding numbers of very different
#' magnitude we iteratively add the numbers pairwise. Because of nummerical issues with doing this, the order
#' of the input values can affect the result.
#' 
#' @param log_values    Sequence of numbers in log space \eqn{[\log(x_1), \log(x_2), ..., \log(x_n)]}
#' @return              \eqn{\log(\sum_{i=1}^n x_i)}
#' 
#' @export
log_sum_of_logs <- function(log_values) {

  # computes log(x),log(y) |-> log(x+y)
  log_add <- function(x, y) {
    if (x == -Inf)
      return(y)
    if (y == -Inf)
      return(x)
    max(x, y) + log1p(exp( -abs(x - y) ))
  }
  
  # We do this for short enough vectors to make the rest work
  # It won't quite work for long sequences because we will end up adding
  # very small numbers to very large numbers. There we want to add them pairwise
  # so we keep the magnitude of the numbers roughly equal.
  if (length(log_values) < 3)
    return(Reduce(log_add, log_values))
  
  # get an even number of values
  if (length(log_values) %% 2 == 1) {
    log_values <- 
      c(log_add(log_values[1], log_values[2]), log_values[3:length(log_values)])
  }
  
  # permute to avoid problems with low and high values clustering
  log_values <- sample(log_values, replace = FALSE)
  while(length(log_values) > 1) {
    #cat(log_values[1:min(10,length(log_values))], "\n")
    indices <- 2*seq_len(length(log_values) / 2) - 1
    log_values <- 
      unlist(Map(function(idx) log_add(log_values[idx], log_values[idx+1]), indices))  
  }
  log_values
}

#' Computes the likelihood of a model from samples from its posterior distribution.
#' 
#' The likelihood of a graph can be computed by integrating over all the graph parameters (with appropriate priors).
#' Doing this by sampling from priors is very inefficient, so we use samples from the posteriors to importance
#' sample the likelihood.
#' 
#' @param log_likelihoods    Samples of log likelihoods from the posterior distribution of the graph.
#' @return                   The likelihood of a graph where graph parameters are integrated out.
#' @export
model_likelihood <- function(log_likelihoods) {
  log_mean_inverse_log <- log_sum_of_logs(-log_likelihoods) - log(length(log_likelihoods))
  -log_mean_inverse_log
}


#' Computes the likelihood of a model from samples from its posterior distribution.
#' 
#' The likelihood of a graph can be computed by integrating over all the graph parameters (with appropriate priors).
#' Doing this by sampling from priors is very inefficient, so we use samples from the posteriors to importance
#' sample the likelihood.
#' 
#' The numerical issues with adding a lot of numbers in log space is unstable
#' so we get a better estimate by doing it several times on different permutations
#' of the data.This function calculates the mean of the likelihoods over different permutations of the
#' input and estimates the standard devition.
#'
#' @param log_likelihoods    Samples of log likelihoods from the posterior distribution of the graph.
#' @param no_samples         Number of permutations to sample when computing the result.
#' @return                   The likelihood of a graph where graph parameters are integrated out given as the mean and standard
#'                           deviation over \code{no_samples} different permutations of the input.
#' @export
model_likelihood_n <- function(log_likelihoods, no_samples = 100) {
  samples <- replicate(no_samples, model_likelihood(log_likelihoods))
  cbind(mean = mean(samples), sd = stats::sd(samples))
}

#' Computes the Bayes factor between two models from samples from their posterior distributions.
#' 
#' The likelihood of a graph can be computed by integrating over all the graph parameters (with appropriate priors).
#' Doing this by sampling from priors is very inefficient, so we use samples from the posteriors to importance
#' sample the likelihood. Given two graphs, and samples from their posteriors, we can estimate the Bayes factor
#' between them.
#' 
#' The numerical issues with adding a lot of numbers in log space is unstable
#' so we get a better estimate by doing it several times on different permutations
#' of the data. This function calculates the mean of the Bayes factors over different permutations of the
#' input and estimates the standard devition.
#'
#' @param logL1              Samples of log likelihoods from the posterior distribution of the first graph.
#' @param logL2              Samples of log likelihoods from the posterior distribution of the second graph.
#' @param no_samples         Number of permutations to sample when computing the result.
#' @return                   The Bayes factor between the two graphs given as the mean and standard
#'                           deviation over \code{no_samples} different permutations of the input.
#' @export
model_bayes_factor_n <- function(logL1, logL2, no_samples = 100) {
  samples <- replicate(no_samples, model_likelihood(logL1) - model_likelihood(logL2))
  cbind(mean = mean(samples), sd = stats::sd(samples))
}

