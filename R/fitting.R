# Kalle's suggestions to improve the model:
#
# 1) Now the bear example for instance basically has 2 equations and 2 edge
#    variables for any pair of admix variables. This means that the square sum
#    error is usually the same, unless the admix variables are such that optimal
#    edges would be negative (in principle we could have also a singularity with
#    like probability 0). So we at best get some bounds for the admix proportions,
#    the suggestions nelder mead gives are arbitrary. Thos bounds correspond to
#    cases where some edges are zero, which I assume is often nonsensical, so one
#    improvement could be setting lower bounds to edge lenghts based on someting
#    we know about drift?
# 2) Actually we're just missing some data. The system would be trurly 
#    underdetermined if we had f4(BLK, PB; A, B/C) statistics too. I either case,
#    trurly underdetermined system or not, we want to make sure the program
#    recognizes when the admix solution is really optimal and when many, possibly
#    all admixes perform equally well. A simple contour plot could be in order.
# 3) We know that the expected values of f4:s can be calculated from the graph
#    overlaps, but what about the standard deviations? They should really be
#    taken into account as weights for each equation in the system.
# 4) Really the whole model could in fact be a Bayesian model giving equal a 
#    priori weight to all admix proportions and then conditioning with reality
#    to pick the favourite candidate. Could be also that we do not know enough
#    about the random variables, right now only expected values.

build_edge_optimisation_matrix <- function(data, graph, parameters = extract_graph_parameters(graph)) {
  m <- nrow(data) # No. of equations is the no. of f4-statistics.
  n <- length(parameters$edges) # Variables are the edges.
  edge_optimisation_matrix <- matrix("0", m, n)
  # Let's fill the matrix with polynomials of admix proportions.
  for (i in seq(1, m)) {
    statistic <- f4(graph, data[i, 1], data[i, 2], data[i, 3], data[i, 4])
    for (j in seq(1, length(statistic))) {
      if (length(statistic[[j]]$prob) != 0) {
        admix_product <- ""
        for (k in seq(1, length(statistic[[j]]$prob))) {
          admix_product <- paste(admix_product, statistic[[j]]$prob[k], sep = "*")
        }
        admix_product <- substring(admix_product, 2)
        # Yeah I know this is a bit silly but the matrix is only created once.
        if (nrow(statistic[[j]]$positive) > 0) { # Insert the positive stuff
          for (k in seq(1, nrow(statistic[[j]]$positive))) {
            edge_name <- paste("edge", statistic[[j]]$positive[k, 1], statistic[[j]]$positive[k, 2], sep = "_")
            edge_number <- match(edge_name, parameters$edges)
            edge_optimisation_matrix[i, edge_number] <- paste(edge_optimisation_matrix[i, edge_number], admix_product, sep = "+")
          }
        }
        if (nrow(statistic[[j]]$negative) > 0) { # Insert the negative stuff
          for (k in seq(1, nrow(statistic[[j]]$negative))) {
            edge_name <- paste("edge", statistic[[j]]$negative[k, 1], statistic[[j]]$negative[k, 2], sep = "_")
            edge_number <- match(edge_name, parameters$edges)
            edge_optimisation_matrix[i, edge_number] <- paste(edge_optimisation_matrix[i, edge_number], admix_product, sep = "-")
          }
        }
      }
    }
  }
  # To prevent the unnecessary "0+":s and "0":s from slowing down evaluation later 
  # (don't know if serious actually), the final cleaning:
  for (i in seq(1,m)) {
    for (j in seq(1,n)) {
      if (edge_optimisation_matrix[i, j] != "0") {
        edge_optimisation_matrix[i, j] <- substring(edge_optimisation_matrix[i, j], 2)
        if (substring(edge_optimisation_matrix[i, j], 1, 1) == "+") {
          edge_optimisation_matrix[i, j] <- substring(edge_optimisation_matrix[i, j], 2)
        }
      }
    }
  }
  edge_optimisation_matrix
}

edge_optimisation_function <- function(data, matrix, graph, parameters = extract_graph_parameters(graph)) {
  if (!requireNamespace("pracma", quietly = TRUE)) {
    stop("This function requires pracma to be installed.")
  }
  force(data)
  force(matrix) # Don't know what this does but it sounds pretty cool.
  force(parameters)
  goal <- data$D
  function(x) {
    # First we evaluate the edge otimisation matrix at (admix-proportions) x
    evaluated_matrix <- matrix(0, nrow(matrix), ncol(matrix))
    for (i in seq(1, length(parameters$admix_prop))) {
      assign(parameters$admix_prop[i], x[i])
    }
    for (i in seq(1, nrow(matrix))) {
      for (j in seq(1, ncol(matrix))) {
        evaluated_matrix[i,j] <- eval(parse(text = matrix[i, j]))
      }
    }
    # It's useful to know which edge lengths are in fact free variables and
    # what edge lengths depend on one another, as the following least square
    # function just gives one exaple of an optimal solution. This takes O(n^3)
    # like the current implementation of the least squares algorithm, so making
    # this extra effort doesn't matter much. Except it's done every time and we
    # only need it for the final answer. TODO FIX
    homogeneous_matrix <- rref(evaluated_matrix)
    # Now just use a ready-made function to find the best non-negative solution
    # in the Euclidian norm. Apparently this is "slow" in the sense it takes 
    # (n^3) steps and not ~O(n^2.3) steps as it coud in principle.
    lsq_solution <- lsqnonneg(evaluated_matrix, goal)
    list(lsq_norm = lsq_solution$resid.norm, lsq_x = lsq_solution$x, homogeneous = homogeneous_matrix)
  }
}

cost_function <- function(data, matrix, graph, parameters = extract_graph_parameters(graph)) {
  function(x) {
    edge_optimisation_function(data, matrix, graph, parameters)(x)$lsq_norm
  }
}

# For testing purposes, should work whenever there's exactly 2 admix proportions.
# With more than 2 admix proportions drawing a graph becomes hard anyways.
contour_plot <- function(data, matrix, graph, parameters = extract_graph_parameters(graph)) {
  resolution <- 10
  x <- seq(0, resolution)
  y <- seq(0, resolution)
  z <- matrix(0, resolution, resolution)
  for (i in seq (1, resolution)) {
    for (j in seq(1, resolution)) {
      z[i, j] <- cost_function(data, matrix, graph, parameters)(c(i/resolution, j/resolution))
    }  
  }
  x <- 10*1:nrow(z)
  y <- 10*1:ncol(z)
  require(grDevices) # for colours
  filled.contour(x, y, z, color = heat.colors)
}

## Graph fitting ##################################################################

#' Fit the graph parameters to a data set.
#' 
#' Tries to minimize the squared distance between statistics in \code{data} and 
#' statistics given by the graph.
#' 
#' The data frame, \code{data}, must contain columns \code{W}, \code{X}, 
#' \code{Y}, and \code{Z}. The function then computes the \eqn{f_4(W,X;Y,Z)}
#' statistics for all rows from these to obtain the prediction made by the
#' graph.
#' 
#' The data frame must also contain a column, \code{D}, containing the 
#' statistics observed in the data. The fitting algorithm attempts to minimize 
#' the distance from this column and the predictions made by the graph.
#' 
#' @param data  The data set.
#' @param graph The admixture graph.
#' @param parameters In case one wants to tweak something in the graph.
#' @param optimisation_options Options to the optimisation algorithm.
#'   
#' @return A list containing the best fitted values (in an environment) and the 
#'   data extended with a column containing the graph predictions.
#'   
#' @seealso \code{\link[neldermead]{optimset}}
#'   
#' @export
fit_graph <- function(data, graph, parameters = extract_graph_parameters(graph), optimisation_options = NULL) {
  if (!requireNamespace("neldermead", quietly = TRUE)) {
    stop("This function requires neldermead to be installed.")
  }
  x0 <- rep(0.5, length(parameters$admix_prop))
  matrix <- build_edge_optimisation_matrix(data, graph, parameters)
  cfunc <- cost_function(data, matrix, graph, parameters)
  opti <- neldermead::fminbnd(cfunc, x0 = x0, xmin = rep(0, length(x0)),
                               xmax = rep(1, length(x0)),
                               options = optimisation_options)
  # The value opti is a class "neldermead" object.
  best_fit <- neldermead::neldermead.get(opti, "xopt") # Optimal parameters.
# The thing is, I'm confused with the environment object adn how to use it so at
# the moment the function just returns a list. I promise to fix back to the
# environment system if that is the true way.
  list(
    data = data,
    graph = graph,
    params = parameters,
    best_fit = best_fit,
    best_edge_fit = edge_optimisation_function(data, matrix, graph, parameters)(best_fit)$lsq_x,
    best_error = edge_optimisation_function(data, matrix, graph, parameters)(best_fit)$lsq_norm,
    homogeneous = edge_optimisation_function(data, matrix, graph, parameters)(best_fit)$homogeneous
  )  
}

# Did not intentionally break anything below this, of course nothing will work
# since I destroyed the environment.

## Interface for accessing fitted data ############################################

#' Print function for a fitted graph.
#' 
#' Print summary of the result of a fit.
#' 
#' @param x  The fitted object.
#' @param ...     Additional arguments.
#'  
#' @export
print.agraph_fit <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat("Sum of squared error:", x$error, "\n")
  invisible(x)
}

fitted_parameters <- function(object) {
  edges <- unlist(Map(function(param) get(param, object$fit_env), object$param$edges))
  admixture_proportions <- unlist(Map(function(param) get(param, object$fit_env), object$param$admix_prop))
  list(edges = edges, admixture_proportions = admixture_proportions)
}

#' Get fitted parameters for a fitted graph.
#' 
#' Extract the graph parameters for a graph fitted to data.
#' 
#' @param object  The fitted object.
#' @param ...     Additional arguments.
#' 
#' @export
coef.agraph_fit <- function(object, ...) {
  parameters <- fitted_parameters(object)
  c(parameters$edges, parameters$admixture_proportions)
}

#' Print function for a fitted graph.
#' 
#' Print summary of the result of a fit.
#' 
#' @param object  The fitted object.
#' @param ...     Additional arguments.
#' 
#' @export
summary.agraph_fit <- function(object, ...) {
  result <- fitted_parameters(object)
  print(result)
  invisible(result)
}


#' Extract fitted data for a fitted graph.
#' 
#' Get the predicted f4 statistics for a fitted graph.
#' 
#' @param object  The fitted object.
#' @param full    Should the fitted values include the full data used for fitting?
#' @param ...     Additional arguments.
#' 
#' @export
fitted.agraph_fit <- function(object, full = TRUE, ...) {
  if (full) object$fit_data
  else      object$fit_data$graph_f4
}

#' Extract the individual errors in a fitted graph.
#' 
#' Get D - graph_f4 for each data point used in the fit.
#' 
#' @param object  The fitted object.
#' @param ...     Additional arguments.
#' 
#' @export
residuals.agraph_fit <- function(object, ...) {
  with(object$fit_data, D - graph_f4)
}

#' Predict statistics on new data.
#' 
#' Predict expected f4 statistics. If \code{newdata} is not specified this function just
#' returns the predicted values on the original data.
#' 
#' @param object  The fitted object.
#' @param newdata New data frame to predict values for.
#' @param ...     Additional arguments.
#' 
#' @export
predict.agraph_fit <- function(object, newdata = NULL, ...) {
  if (is.null(newdata)) {
    fitted(object, full = TRUE)
  } else {
    add_graph_f4(newdata, object$graph, object$fit_env)
  }
}

