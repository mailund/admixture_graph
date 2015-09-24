# Kalle's suggestions to improve the model:
#
# 1) Take the standard deviations into account. One way is imposing weight
#    coefficients into the distance function. Maybe even think what is the effect
#    of assuming the data is normally distributed.
# 4) Really the whole model could perhaps be a Bayesian model giving equal a 
#    priori weight to all admix proportions and then conditioning with reality
#    to pick the favourite candidate. Could be also that we do not know enough
#    about the random variables?

## Graph fitting #################################################################

#' Build a matrix coding the linear system of edges once the admix variables
#' have been fixed.
#' 
#' The elements are characters containing numerals, admix variable names,
#' parenthesis and arithmetical operations. (Transform into expressions with
#' \code{parse} and then evaluate with \code{eval}). The column names are the
#' edge names from \code{extract_graph_parameter$edges}, the rows have no names.
#' 
#' If the essential number of equations is not higher than the essential number of
#' edge variables, the quality of edge optimisation will not depend on the admix
#' variables (expect in a very special cases), and a complaint will be given.
#' 
#' @param data        The data set.
#' @param graph       The admixture graph.
#' @param parameters  In case one wants to tweak something in the graph.
#'   
#' @return A list containing the full matrix (\code{$full}), a version with zero
#'         columns removed (\code{$column_reduced}), a version with zero rows and
#'         repeated rows also removed (\code{$double_reduced}), and an indicator
#'         of warning (\code{$complaint}).
build_edge_optimisation_matrix <- function(data, graph, parameters 
                                           = extract_graph_parameters(graph)) {
  m <- nrow(data) # Number of equations is the number of f4-statistics.
  n <- length(parameters$edges) # Variables are the edges.
  edge_optimisation_matrix <- matrix("0", m, n)
  colnames(edge_optimisation_matrix) <- parameters$edges
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
            edge_name <- paste("edge", statistic[[j]]$positive[k, 1],
                               statistic[[j]]$positive[k, 2], sep = "_")
            edge_optimisation_matrix[i, edge_name] <- 
              paste(edge_optimisation_matrix[i, edge_name],
                     admix_product, sep = " + ")
          }
        }
        if (nrow(statistic[[j]]$negative) > 0) { # Insert the negative stuff
          for (k in seq(1, nrow(statistic[[j]]$negative))) {
            edge_name <- paste("edge", statistic[[j]]$negative[k, 1], 
                               statistic[[j]]$negative[k, 2], sep = "_")
            edge_optimisation_matrix[i, edge_name] <- 
              paste(edge_optimisation_matrix[i, edge_name], 
                     admix_product, sep = " - ")
          }
        }
      }
    }
  }
  # To prevent the unnecessary "0 +":s and "0":s from slowing down evaluation
  # later (don't know if serious actually), perform some cleaning:
  for (i in seq(1,m)) {
    for (j in seq(1,n)) {
      if (edge_optimisation_matrix[i, j] != "0") {
        edge_optimisation_matrix[i, j] <- 
          substring(edge_optimisation_matrix[i, j], 3)
        if (substring(edge_optimisation_matrix[i, j], 1, 1) == "+") {
          edge_optimisation_matrix[i, j] <- 
            substring(edge_optimisation_matrix[i, j], 3)
        }
      }
    }
  }
  # Make a version with zero columns removed.
  column_reduced <- edge_optimisation_matrix
  j <- 1
  while (j <= ncol(column_reduced)) {
    for (i in seq(1, m)) {
      if (column_reduced[i, j] != "0") {
        j <- j + 1
        break
      }
      if (i == m) {
        column_reduced <- column_reduced[, -j]
      }
    }
  }
  # Make a version with repeated rows and zero rows removed.
  double_reduced <- column_reduced
  i <- 1
  while (i <= nrow(double_reduced)) {
    for (j in seq(ncol(double_reduced))) {
      if (double_reduced[i, j] != "0") {
        i <- i + 1
        break
      }
      if (j == ncol(double_reduced)) {
        double_reduced <- double_reduced[-i,]
      }
    }
  }
  double_reduced <- unique(double_reduced)
  # Make a complaint if the double reduced matrix is not high.
  complaint <- FALSE
  if (nrow(double_reduced) <= ncol(double_reduced)) {
    complaint <- TRUE
  }
  list(full = edge_optimisation_matrix, column_reduced = column_reduced, 
       double_reduced = double_reduced, complaint = complaint)
}

#' The cost function fed to nelder mead.
#' 
#' We want nelder mead to run fast so the cost function operates with the column
#' rduced edge optimisation matrix and does not give any extar information about 
#' the fit. For the details, use \code{edge_optimisation_function} instead.
#' 
#' @param data  The data set.
#' @param matrix  A column reduced edge optimisation matrix (typically given by 
#'                the function \code{edge_optimisation_matrix$column_reduced}).
#' @param graph  The admixture graph.
#' @param parameters  In case one wants to tweak something in the graph.
#'   
#' @return  Given an input vector of admix variables, returns the smallest error 
#'          regarding the edge variables.
cost_function <- function(data, matrix, graph, 
                          parameters = extract_graph_parameters(graph)) {
  if (!requireNamespace("pracma", quietly = TRUE)) {
    stop("This function requires pracma to be installed.")
  }
  goal <- data$D
  function(x) {
    # Evaluate the column reduced edge optimisation matrix at x.
    evaluated_matrix <- matrix(0, nrow(matrix), ncol(matrix))
    for (i in seq(1, length(parameters$admix_prop))) {
      assign(parameters$admix_prop[i], x[i])
    }
    for (i in seq(1, nrow(matrix))) {
      for (j in seq(1, ncol(matrix))) {
        evaluated_matrix[i, j] <- eval(parse(text = matrix[i, j]))
      }
    }
    # Now just use a ready-made function to find the best non-negative solution
    # in the Euclidian norm. Apparently this is "slow" in the sense it takes 
    # O(n^3) steps and not O(n^2.3) steps as it could in principle.
    lsq_solution <- pracma::lsqnonneg(evaluated_matrix, goal)
    lsq_solution$resid.norm
  }
}

#' More detailed edge fitting than mere \code{cost_function}.
#' 
#' Returning the cost, an example edge solution of an optimal fit, and linear 
#' relations describing the set of all edge solutions. Operating with the full
#' edge optimisation matrix.
#' 
#' @param data  The data set.
#' @param matrix  A full  edge optimisation matrix (typically given by the 
#'                function \code{edge_optimisation_matrix$full}).
#' @param graph  The admixture graph.
#' @param parameters  In case one wants to tweak something in the graph.
#'   
#' @return  Given an input vector of admix variables, returns a list \code{x} containing
#'          the minimal error (\code{x$cost}), the graph-f4-statistics 
#'          (\code{x$approximation}), an example solution (\code{x$edge_fit}), linear
#'          relations describing all the solutions (\code{x$homogeneous}) and one 
#'          way to choose the free (\code{x$free_edges}) and bounded 
#'          (\code{x$bounded_edges}) edge variables.
edge_optimisation_function <- function(data, matrix, graph, 
                              parameters = extract_graph_parameters(graph)) {
  if (!requireNamespace("pracma", quietly = TRUE)) {
    stop("This function requires pracma to be installed.")
  }
  goal <- data$D
  function(x) {
    # Evaluate the full edge otimisation matrix at x.
    evaluated_matrix <- matrix(0, nrow(matrix), ncol(matrix))
    for (i in seq(1, length(parameters$admix_prop))) {
      assign(parameters$admix_prop[i], x[i])
    }
    for (i in seq(1, nrow(matrix))) {
      for (j in seq(1, ncol(matrix))) {
        evaluated_matrix[i, j] <- eval(parse(text = matrix[i, j]))
      }
    }
    # Record the (or an example of an) optimal solution and error.
    lsq_solution <- pracma::lsqnonneg(evaluated_matrix, goal)
    edge_fit <- lsq_solution$x
    names(edge_fit) <- parameters$edges
    approximation <- evaluated_matrix %*% edge_fit
    approximation <- approximation[, 1]
    # It's useful to know which edge lengths are in fact free variables and what
    # edge lengths depend on one another, as the least square function only gave
    # one exaple of an optimal solution. This information is visible after
    # manipulating the optimisation matrix into reduced row echelon form.
    homogeneous_matrix <- pracma::rref(evaluated_matrix)
    # Make a list of (one choice of) free edges.
    free_edges <- c()
    i <- 1
    j <- 1
    while (j <= ncol(matrix)) {
      if (homogeneous_matrix[i, j] != 0) {
        if (i == nrow(matrix) && j < ncol(matrix)) {
          for (k in seq(j + 1, ncol(matrix))) {
            free_edges <- c(free_edges, parameters$edges[k])
          }
          break
        }
        i <- i + 1
        j <- j + 1
      }
      else {
        free_edges <- c(free_edges, parameters$edges[j])
        j <- j + 1
      }
    }
    # Explain the relationship between the remaining edges and the free edges.
    h <- nrow(unique(rbind(homogeneous_matrix, 0))) - 1
    bounded_edges <- rep("", h)
    for (i in seq(1, h)) {
      for (j in seq(1, ncol(matrix))) {
        if (homogeneous_matrix[i, j] != 0) {
          if (bounded_edges[i] == "") {
            bounded_edges[i] <- paste(bounded_edges[i], parameters$edges[j], " =",
                                      sep = "")
          }
          else {
            if (homogeneous_matrix[i, j] > 0) {
              bounded_edges[i] <- paste(bounded_edges[i], " - ",
                                        homogeneous_matrix[i, j], "*",
                                        parameters$edges[j], sep = "")
            }
            if (homogeneous_matrix[i, j] < 0) {
              if (substring(bounded_edges[i], nchar(bounded_edges[i])) == "=") {
                bounded_edges[i] <- paste(bounded_edges[i], " ", sep = "")
              }
              else {
                bounded_edges[i] <- paste(bounded_edges[i], " + ", sep = "")
              }
              bounded_edges[i] <- paste(bounded_edges[i], homogeneous_matrix[i, j],
                                        "*", parameters$edges[j], sep = "")
            }
          }
        }
      }
      if (substring(bounded_edges[i], nchar(bounded_edges[i])) == "=") {
        bounded_edges[i] <- paste(bounded_edges[i], " 0", sep = "")
      }
    }
    list(cost = lsq_solution$resid.norm, approximation = approximation,
         edge_fit = edge_fit, homogeneous = homogeneous_matrix, 
         free_edges = free_edges, bounded_edges = bounded_edges)
  }
}

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
#' @param graph  The admixture graph.
#' @param optimisation_options  Options to the optimisation algorithm.
#' @param parameters  In case one wants to tweak something in the graph.
#'   
#' @return A list containing everything about the fit.
#'   
#' @seealso \code{\link[neldermead]{optimset}}
#'   
#' @export
fit_graph <- function(data, graph, optimisation_options = NULL,
                      parameters = extract_graph_parameters(graph)) {
  if (!requireNamespace("neldermead", quietly = TRUE)) {
    stop("This function requires neldermead to be installed.")
  }
  x0 <- rep(0.5, length(parameters$admix_prop))
  matrix <- build_edge_optimisation_matrix(data, graph, parameters)
  full_matrix <- matrix$full
  reduced_matrix <- matrix$column_reduced
  cfunc <- cost_function(data, reduced_matrix, graph, parameters)
  opti <- neldermead::fminbnd(cfunc, x0 = x0, xmin = rep(0, length(x0)),
                              xmax = rep(1, length(x0)),
                              options = optimisation_options)
  # The value opti is a class "neldermead" object.
  best_fit <- neldermead::neldermead.get(opti, "xopt") # Optimal admix values.
  best_fit <- best_fit[, 1]
  names(best_fit) <- parameters$admix_prop
  detailed_fit <- 
    edge_optimisation_function(data, full_matrix, graph, parameters)(best_fit)
  data$graph_f4 <- detailed_fit$approximation
  # The output is a list with "agraph_fit" -mystery property.
  structure(list(
      call = sys.call(),
      data = data,
      graph = graph,
      matrix = matrix,
      complaint = matrix$complaint,
      best_fit = best_fit,
      best_edge_fit = detailed_fit$edge_fit,
      homogeneous = detailed_fit$homogeneous,
      free_edges = detailed_fit$free_edges,
      bounded_edges = detailed_fit$bounded_edges,
      best_error = detailed_fit$cost,
      approximation = detailed_fit$approximation,
      parameters = parameters
    ),  
    class = "agraph_fit"
  )
}

## Interface for accessing fitted data ############################################

#' Print function for a fitted graph.
#' 
#' Print summary of the result of a fit.
#' 
#' @param x       The fitted object.
#' @param ...     Additional parameters.
#'  
#' @export
print.agraph_fit <- function(x, ...) {
  cat("\n")
  cat("Call:")
  cat("\n")
  print(x$call)
  
  if (x$complaint == TRUE) {
    cat("\n")
    cat("The data is not sufficient to give a meaningful fit for this topology!")
    cat("\n")
  }
  
  cat("Minimal error:", x$best_error)
}

#' Get fitted parameters for a fitted graph.
#' 
#' Extract the graph parameters for a graph fitted to data.
#' 
#' @param object  The fitted object.
#' @param ...     Additional parameters.
#' 
#' @export
coef.agraph_fit <- function(object, ...) {
  c(object$best_edge_fit, object$best_fit)
}

#' Print function for a fitted graph.
#'   
#' Print summary of the result of a fit.
#' 
#' @param object  The fitted object.
#' @param ...     Additional parameters.
#' 
#' @export
summary.agraph_fit <- function(object, ...) {
  cat("\n")
  cat("Call:")
  cat("\n")
  print(object$call)
  if (object$complaint == TRUE) {
    cat("\n")
    cat("The data is not sufficient to give a meaningful fit for this topology!")
    cat("\n")
  }
  cat("\n")
  cat("Optimal admix variables:")
  cat("\n")
  print(object$best_fit)
  cat("\n")
  cat("Optimal edge variables:")
  cat("\n")
  print(object$best_edge_fit)
  cat("\n")
  cat("Solution to a homogeneous system of edges with the optimal admix variables:")
  cat("\n")
  cat("(Adding any such solution to the optimal one will not affect the error.)")
  cat("\n")
  cat("\n")
  cat("Free edge variables:")
  cat("\n")
  cat(object$free_edges, sep = "\n")
  cat("\n")
  cat("Bounded edge variables:")
  cat("\n")
  cat(object$bounded_edges, sep = "\n")
  cat("\n")
  cat("Minimal error:")
  cat("\n")
  cat(object$best_error)
}

#' Extract fitted data for a fitted graph.
#' 
#' Get the predicted f4 statistics for a fitted graph.
#' 
#' @param object  The fitted object.
#' @param ...     Additional parameters.
#' 
#' @export
fitted.agraph_fit <- function(object, ...) {
  object$data
}

#' Extract the individual errors in a fitted graph.
#' 
#' Get D - graph_f4 for each data point used in the fit.
#' 
#' @param object  The fitted object.
#' @param ...     Additional parameters.
#' 
#' @export
residuals.agraph_fit <- function(object, ...) {
  object$data$D - object$data$graph_f4
}
