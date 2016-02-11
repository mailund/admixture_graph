
make_mcmc_model <- function(graph, data) {
  
  
  
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
  
  log_likelihood <- function(state, data) {
    # FIXME: this is where we need Kalle's likelihood function
    
    # this is a mock up
    sum(log(dnorm(state)))
  }
  
  log_posterior <- function(state, data) {
    log_prior(state) + log_likelihood(state, data)
  }

  list(log_posterior = log_posterior)
}

run_metropolis_hasting <- function(model, initial_state, iterations, ...) {
  if (!requireNamespace("MHadaptive", quietly = TRUE)) {
    stop("This function requires MHadaptive to be installed.")
  }
  
  MHadaptive::Metro_Hastings(li_func = model$log_posterior, 
                             pars = initial_state, 
                             par_names = letters[1:5], # FIXME
                             iterations = iterations,
                             ...)
}


model <- make_model()
mcmc_r <- run_metropolis_hasting(model, rep(0.5, 5), 10000) %>% mcmc_thin(thin = 100)
plotMH(mcmc_r)
