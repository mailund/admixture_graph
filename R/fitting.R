# TODO:
#
# 1) Take the covariance matrix into account. When covariance matrix not given,
#    assume independence and build one from data$Z.
# 2) Find out how the f-statistics depend on one another. Find out what subsets of
#    the set of all statistics have no implied weight on some data. Preferably make
#    a program to remove extra data from already sufficient data in some fair way.

## Graph fitting #################################################################

#' Used to recognize similar expressions and to possibly simplify them.
#'
#' This is not pretty but let's see if it speeds up the program.
#'
#' @param x   Input is assumed to be a char consisting of one or more terms. Each term
#'            starts with either \code{+} or \code{-} and after that contains one or
#'            more factors separated by \code{*}. Each factor is either an admix variable,
#'            a number or a clause \code{(1 - x)} (mind the spaces, this is how the
#'            f4-function outputs), where x is again either an admix variable or a
#'            number.
#'            Everything is pretty much ruined if variable names are numbers or contain
#'            forbidden symbols \code{+, -, *, (, )}.
#'
#' @return   A polynomial in a canonical form with no parenthesis or spaces and the
#'           monomials in lexicographical order. If everything is cancelled out then "+0".
#'
#' @export
canonise_expression <- function(x) {
  # First remove the symbols 1 and - from inside the parenthesis to make things easier.
  l <- nchar(x)
  for (i in seq(0, l - 1)) {
   if (substring(x, l - i, l - i) == "(") {
      x <- paste(substring(x, 1, l - i), substring(x, l - i + 5), sep = "")
    }
  }
  # Then copy each term to a list of terms.
  terms_list <- list()
  start <- 1
  for (i in seq(2, nchar(x))) {
    if (substring(x, i, i) == "+" || substring(x, i, i) == "-") {
      terms_list <- c(terms_list, substring(x, start, i - 1))
      start <- i
    }
  }
  terms_list <- c(terms_list, substring(x, start, nchar(x)))
  # Now recursively multiply the parentheses open collecting the resulting monomials
  # to a list.
  monomials_list <- list()
  while (length(terms_list) > 0) {
    term <- terms_list[1]
    terms_list[1] <- NULL
    left <- 0
    right <- 0
    for (i in seq(1, nchar(term))) {
      if (substring(term, i, i) == "(") {
        left <- i
      }
      if (substring(term, i, i) == ")") {
        right <- i
      }
    }
    if (left == 0) {
      monomials_list <- c(monomials_list, term)
    }
    else {
      newterm1 <- paste(substring(term, 1, left - 1), "1",
                        substring(term, right + 1, nchar(term)), sep = "")
      if (substring(term, 1, 1) == "+") {
        newterm2 <- paste("-", substring(term, 2, left - 1),
                          substring(term, left + 1, right - 1),
                          substring(term, right + 1, nchar(term)), sep = "")
      }
      if (substring(term, 1, 1) == "-") {
        newterm2 <- paste("+", substring(term, 2, left - 1),
                          substring(term, left + 1, right - 1),
                          substring(term, right + 1, nchar(term)), sep = "")
      }
      terms_list <- c(terms_list, newterm1, newterm2)
    }
  }
  # It is time to recognize numbers and do the appropriate arithmetics. Save the
  # monomials as pairs of a numerical coefficient and a product of variables. These
  # can then be added together.
  math_list <- list()
  for (i in seq(1, length(monomials_list))) {
    coefficient <- 1
    variables <- ""
    temp_vector <- numeric(0)
    monomial <- monomials_list[i]
    if (substring(monomial, 1, 1) == "-") {
      coefficient <- - 1
    }
    start <- 2
    for (j in seq(2, nchar(monomial))) {
      if (substring(monomial, j, j) == "*") {
        word <- substring(monomial, start, j - 1)
        if (suppressWarnings(!is.na(as.numeric(word))) == TRUE) {
          coefficient <- coefficient * suppressWarnings(as.numeric(word))
        }
        else {
          temp_vector<- c(temp_vector, word)
        }
        start <- j + 1
      }
    }
    word <- substring(monomial, start, nchar(monomial))
    if (suppressWarnings(!is.na(as.numeric(word))) == TRUE) {
      coefficient <- coefficient * suppressWarnings(as.numeric(word))
    }
    else {
      temp_vector <- c(temp_vector, word)
    }
    if (length(temp_vector) > 0) {
      temp_vector <- sort(temp_vector)
      for (j in seq(1, length(temp_vector))) {
        variables <- paste(variables, temp_vector[j], sep = "*")
      }
      variables <- substring(variables, 2)
    }
    else {
      variables <- "1"
    }
    pair <- list(num = coefficient, let = variables)
    math_list[[i]] <- pair
  }
  i <- 1
  while (i < length(math_list)) {
    j <- i + 1
    while (j <= length(math_list)) {
      if (math_list[[i]]$let == math_list[[j]]$let) {
        math_list[[i]]$num <- math_list[[i]]$num + math_list[[j]]$num
        math_list[[j]] <- NULL
      }
      else {
        j <- j + 1
      }
    }
    i <- i + 1
  }
  l <- length(math_list)
  for (i in seq(0, l - 1)) {
    if (math_list[[l - i]]$num == 0) {
      math_list [[l - i]] <- NULL
    }
  }
  # Put elements in lexicographic order, unless the list in empty in which case we want
  # to return "+0" for later use (because as.numeric() does not interpret "" as zero).
  if (length(math_list) == 0) {
    result <- "+0"
  } else {
    result <- ""
    nums <- sapply(math_list,"[[","num")
    lets <- sapply(math_list,"[[","let")
    numerals <- nums[order(lets)]
    letters <- lets[order(lets)]
    # And finally put everything together in a long string.
    for (i in seq(1, length(math_list))) {
      if (numerals[i] < 0) {
        result <- paste(result, numerals[i], "*", letters[i], sep = "")
      }
      if (numerals[i] > 0) {
        result <- paste(result, "+", numerals[i], "*", letters[i], sep = "")
      }
    }
  }
  result
}

#' Build a matrix coding the linear system of edges once the admix variables
#' have been fixed.
#'
#' The elements are characters containing numerals, admix variable names,
#' parenthesis and arithmetical operations. (Transform into expressions with
#' parse and then evaluate with eval). The column names are the edge names from
#' extract_graph_parameter$edges, the rows have no names.
#'
#' If the essential number of equations is not higher than the essential number of
#' edge variables, the quality of edge optimisation will not depend on the admix
#' variables (expect in a very special cases), and a complaint will be given.
#'
#' @param data        The data set.
#' @param graph       The admixture graph.
#' @param parameters  In case one wants to tweak something in the graph.
#'
#' @return A list containing the full matrix ($full), a version with zero columns
#'         removed ($column_reduced) and an indicator of warning ($complaint).
#'
#' @export
build_edge_optimisation_matrix <- function(data, graph, parameters
                                           = extract_graph_parameters(graph)) {
  if (!requireNamespace("pracma", quietly = TRUE)) {
    stop("This function requires pracma to be installed.")
  }
  m <- NROW(data) # Number of equations is the number of f4-statistics.
  n <- length(parameters$edges) # Variables are the edges.
  edge_optimisation_matrix <- matrix("+0", m, n)
  colnames(edge_optimisation_matrix) <- parameters$edges
  # Let's fill the matrix with polynomials of admix proportions.
  for (i in seq(1, m)) {
    statistic <- f4(graph, data$W[i], data$X[i], data$Y[i], data$Z[i])
    if (length(statistic) != 0) {
      for (j in seq(1, length(statistic))) {
        if (length(statistic[[j]]$prob) == 0) {
          statistic[[j]]$prob <- "1"
        }
        admix_product <- ""
        for (k in seq(1, length(statistic[[j]]$prob))) {
          admix_product <- paste(admix_product, statistic[[j]]$prob[k], sep = "*")
        }
        admix_product <- substring(admix_product, 2)
        # Yeah I know this is a bit silly but the matrix is only created once.
        if (NROW(statistic[[j]]$positive) > 0) { # Insert the positive stuff
          for (k in seq(1, NROW(statistic[[j]]$positive))) {
            edge_name1 <- paste("edge", statistic[[j]]$positive[k, 1],
                                statistic[[j]]$positive[k, 2], sep = "_")
            edge_name2 <- paste("edge", statistic[[j]]$positive[k, 2],
                                statistic[[j]]$positive[k, 1], sep = "_")
            if (edge_name1 %in% parameters$edges) {
              edge_name <- edge_name1
            } else {
              edge_name <- edge_name2
            }
            edge_optimisation_matrix[i, edge_name] <-
              paste(edge_optimisation_matrix[i, edge_name],
                    admix_product, sep = "+")
          }
        }
        if (NROW(statistic[[j]]$negative) > 0) { # Insert the negative stuff
          for (k in seq(1, NROW(statistic[[j]]$negative))) {
            edge_name1 <- paste("edge", statistic[[j]]$negative[k, 1],
                                statistic[[j]]$negative[k, 2], sep = "_")
            edge_name2 <- paste("edge", statistic[[j]]$negative[k, 2],
                                statistic[[j]]$negative[k, 1], sep = "_")
            if (edge_name1 %in% parameters$edges) {
              edge_name <- edge_name1
            } else {
              edge_name <- edge_name2
            }
            edge_optimisation_matrix[i, edge_name] <-
              paste(edge_optimisation_matrix[i, edge_name],
                    admix_product, sep = "-")
          }
        }
      }
    }
  }
  # Simplify by putting each element to a canonical form. We will simplify for
  # real later but doing this before any further analysis might speed things up.
  for (i in seq(1, m)) {
    for (j in seq(1, n)) {
      edge_optimisation_matrix[i,j] <-
        canonise_expression(edge_optimisation_matrix[i, j])
    }
  }
  # In order to indentify which admix variables are trurly fitted and which are
  # not, we assign a random value to some of them, not all but maybe none. We can
  # accurately detect if none of the admix variables affect the fit, but if some
  # do and some do not we have a probability zero chance of false accusation.
  # More specifically, we treat some of the admix variables as not variables at all,
  # but constants represented by a random number from the unit interval. With
  # extremely bad luck it's theoretically possible that the random constant chosen
  # decreases the rank of the edge optimisation matrix more than a typical constant
  # would.
  complaint <- integer()
  # First thing to do is to take into account the possibility of no admix variables.
  if (length(parameters$admix_prop) == 0) {
    # Make a version with zero columns removed.
    column_reduced <- edge_optimisation_matrix
    j <- 1
    while (j <= NCOL(column_reduced)) {
      for (i in seq(1, m)) {
        if (column_reduced[i, j] != "+0") {
          j <- j + 1
          break
        }
        if (i == m) {
          column_reduced <- column_reduced[, -j, drop = FALSE]
        }
      }
    }
  } else {
    R <- 2^(length(parameters$admix_prop)) - 1
    weights <- rep(0, R)
    for (r in seq(1, R)) {
      for (a in seq(1, length(parameters$admix_prop))) {
        if ((r/(2^a)) %% 1 < 0.5) {
          weights[r] <- weights[r] + 1
        }
      }
    }
    check_order <- order(weights)
    skip <- rep(FALSE, R)
    for (r in check_order) {
      if (skip[r] == TRUE) {
        for (a in seq(1, length(parameters$admix_prop))) {
          if ((r/(2^a)) %% 1 >= 0.5) {
            skip[r - 2^(a - 1)] <- TRUE
          }
        }
        break
      } 
      A <- rep(NaN, length(parameters$admix_prop))
      for (a in seq(1, length(parameters$admix_prop))) {
        if ((r/(2^a)) %% 1 < 0.5) {
          A[a] <- runif(1) # Note that the last case r = R is not assigning anything.
        }
      }
      evaluated_matrix <- edge_optimisation_matrix # Evaluated but still chars.
      # This might be a silly way for assigning the values but I can't use assign,
      # eval and parse as some variables stay symbolic. The problematic thing about
      # this is that we are no longer allowed to have an admix variable name that is
      # a subword of another.
      for (i in seq(1, m)) {
        for (j in seq(1, n)) {
          for (a in seq(1, length(parameters$admix_prop))) {
            if (is.nan(A[a]) == FALSE) {
              evaluated_matrix[i, j] <-
                gsub(parameters$admix_prop[a], A[a], evaluated_matrix[i, j])
            }
          }
        }
      }
      # Simplify again.
      for (i in seq(1, m)) {
        for (j in seq(1, n)) {
          evaluated_matrix[i,j] <- canonise_expression(evaluated_matrix[i, j])
        }
      }
      # Make a version with zero columns removed.
      column_reduced_temp <- evaluated_matrix
      j <- 1
      while (j <= NCOL(column_reduced_temp)) {
        for (i in seq(1, m)) {
          if (column_reduced_temp[i, j] != "+0") {
            j <- j + 1
            break
          }
          if (i == m) {
            column_reduced_temp <- column_reduced_temp[, -j, drop = FALSE]
          }
        }
      }
      if (r == R) {
        column_reduced <- column_reduced_temp
      }
      # Calculate how many of the polynomials of edge variables and remaining admix
      # variables are linearly independent. (Each different product of the edge variables
      # and the remaining admix variables is treated as a new variable, the already
      # fitted admix variables are treated as real coefficients.) If the number of
      # linearly independent polynomials is not strictly higher than the number of edge
      # variables, we know that no matter what values the remaining admix values obtain,
      # the edges can be chosen in such a way that each independent polynomial gets any
      # value as we please. Thus, the remaining admix variables did not affect the fit and
      # we make a complaint.
      if (NCOL(column_reduced_temp) > 0) {
        big_matrix <- matrix(0, m, 0)
        for (i in seq(1, m)) {
          for (j in seq(1, NCOL(column_reduced_temp))) {
            word <- column_reduced_temp[i, j]
            if (word != "+0") {
              start <- 1
              for (k in seq(1, nchar(word))) {
                if (k == nchar(word) || substring(word, k + 1, k + 1) == "+" ||
                    substring(word, k + 1, k + 1) == "-") {
                  mono <- substring(word, start, k)
                  start <- k + 1
                  star <- which(strsplit(mono, "")[[1]] == "*")[1]
                  number <- as.numeric(substring(mono, 1, star - 1))
                  label <- paste(colnames(column_reduced_temp)[j], substring(mono, star + 1),
                                 sep = "*")
                  if (is.na(match(label, colnames(big_matrix))) == TRUE) {
                    v <- rep(0, m)
                    big_matrix <- cbind(big_matrix, v)
                    colnames(big_matrix)[NCOL(big_matrix)] <- label
                  }
                  big_matrix[i, label] <- number
                }
              }
            }
          }
        }
        h <- qr(big_matrix, tol = 1e-10)$rank
      } else {
        h <- 0
      }
      if (h <= NCOL(column_reduced_temp)) {
        complaint <- c(complaint, r)
        for (a in seq(1, length(parameters$admix_prop))) {
          if ((r/(2^a)) %% 1 >= 0.5) {
            skip[r - 2^(a - 1)] <- TRUE
          }
        }
      }
    }
  }
  list(full = edge_optimisation_matrix, column_reduced = column_reduced,
       complaint = complaint)
}

#' Non negative least square solution.
#'
#' This is lsqnonneg-function from the package pracma, just changed qr.solve into
#' using Moore-Penrose inverse instead (ginv from MASS) as qr.solve crashes for
#' some singular inputs. Now it won't crash but it's sometimes running for very long
#' time (forever?), presumably with those problematic inputs. After too many steps
#' the function halts and lies that the fit was terrible. I don't think this will
#' cause problems.
#'
#' @param C  The matrix.
#' @param d  The vector.
#'
#' @return A vector ($x) and the error ($resid.norm).
#'
#' @export
mynonneg <- function(C, d, iteration_multiplier = 3) {
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("This function requires MASS to be installed.")
  }
  stopifnot(is.numeric(C), is.numeric(d))
  if (!is.matrix(C) || !is.vector(d))
    stop("Argument 'C' must be a matrix, 'd' a vector.")
  m <- NROW(C); n <- NCOL(C)
  if (m != length(d))
    stop("Arguments 'C' and 'd' have nonconformable dimensions.")
  tol <- 10 * 2.220446e-16 * norm(C, "F") * (max(n, m) + 1)
  x  <- rep(0, n)             # initial point
  P  <- logical(n); Z <- !P   # non-active / active columns
  resid <- d - C %*% x
  w <- t(C) %*% resid
  wz <- numeric(n)
  # iteration parameters
  outeriter <- 0; it <- 0
  itmax <- iteration_multiplier * n; exitflag <- 1
  while (any(Z) && any(w[Z] > tol)) {
    outeriter <- outeriter + 1
    z <- numeric(n)
    wz <- rep(-Inf, n)
    wz[Z] <- w[Z]
    im <- which.max(wz)
    P[im] <- TRUE; Z[im] <- FALSE
    z[P] <- MASS::ginv(C[, P]) %*% d
    while (any(z[P] <= 0)) {
      it <- it + 1
      if (it > itmax) {
        warning("Iteration count exceeded.")
        return(list(x = rep(0, n), resid.norm = Inf))
      }
      Q <- (z <= 0) & P
      alpha <- min(x[Q] / (x[Q] - z[Q]))
      x <- x + alpha*(z - x)
      Z <- ((abs(x) < tol) & P) | Z
      P <- !Z
      z <- numeric(n)
      z[P] <- MASS::ginv(C[, P]) %*% d
    }
    x <- z
    resid <- d - C %*% x
    w <- t(C) %*% resid
  }
  return(list(x = x, resid.norm = sum(resid*resid)))
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
#'
#' @export
cost_function <- function(data, matrix, graph,
                          parameters = extract_graph_parameters(graph),
                          iteration_multiplier = 3) {
  if (!requireNamespace("pracma", quietly = TRUE)) {
    stop("This function requires pracma to be installed.")
  }
  goal <- data$D
  function(x) {
    # Evaluate the column reduced edge optimisation matrix at x.
    evaluated_matrix <- matrix(0, NROW(matrix), NCOL(matrix))
    for (i in seq(1, length(parameters$admix_prop))) {
      assign(parameters$admix_prop[i], x[i])
    }
    if (NCOL(matrix) > 0) {
      for (i in seq(1, NROW(matrix))) {
        for (j in seq(1, NCOL(matrix))) {
          evaluated_matrix[i, j] <- eval(parse(text = matrix[i, j]))
        }
      }
    }
    # Now just use a ready-made function to find the best non-negative solution
    # in the Euclidian norm. Apparently this is "slow" in the sense it takes
    # O(n^3) steps and not O(n^2.3) steps as it could in principle.
    if (NCOL(matrix) > 0) {
      lsq_solution <- mynonneg(evaluated_matrix, goal, iteration_multiplier)
      cost <- lsq_solution$resid.norm
    } else {
      cost <- sum(goal^2)
    }
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
#'
#' @export
edge_optimisation_function <- function(data, matrix, graph,
                                       parameters = extract_graph_parameters(graph),
                                       iteration_multiplier = 3) {
  if (!requireNamespace("pracma", quietly = TRUE)) {
    stop("This function requires pracma to be installed.")
  }
  goal <- data$D
  function(x) {
    # Evaluate the full edge otimisation matrix at x, if we even have admix variables.
    evaluated_matrix <- matrix(0, NROW(matrix), NCOL(matrix))
    if (length(parameters$admix_prop) != 0) {
      for (i in seq(1, length(parameters$admix_prop))) {
        assign(parameters$admix_prop[i], x[i])
      }
    }
    for (i in seq(1, NROW(matrix))) {
      for (j in seq(1, NCOL(matrix))) {
        evaluated_matrix[i, j] <- eval(parse(text = matrix[i, j]))
      }
    }
    # Record the (or an example of an) optimal solution and error.
    lsq_solution <- mynonneg(evaluated_matrix, goal, iteration_multiplier)
    edge_fit <- lsq_solution$x
    names(edge_fit) <- parameters$edges
    approximation <- evaluated_matrix %*% edge_fit
    approximation <- approximation[, 1]
    # It's useful to know which edge lengths are in fact free variables and what
    # edge lengths depend on one another, as the least square function only gave
    # one exaple of an optimal solution. This information is visible after
    # manipulating the optimisation matrix into reduced row echelon form.
    if (NROW(evaluated_matrix) != 1) {
      homogeneous_matrix <- pracma::rref(evaluated_matrix)
    } else {
      homogeneous_matrix <- evaluated_matrix
    }
    # Make a list of (one choice of) free edges.
    free_edges <- c()
    i <- 1
    j <- 1
    while (j <= NCOL(matrix)) {
      if (homogeneous_matrix[i, j] != 0) {
        if (i == NROW(matrix) && j < NCOL(matrix)) {
          for (k in seq(j + 1, NCOL(matrix))) {
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
    h <- NROW(unique(rbind(homogeneous_matrix, 0))) - 1
    bounded_edges <- rep("", h)
    if (h != 0) {
      for (i in seq(1, h)) {
        for (j in seq(1, NCOL(matrix))) {
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
                      parameters = extract_graph_parameters(graph),
                      iteration_multiplier = 3) {
  if (!requireNamespace("neldermead", quietly = TRUE)) {
    stop("This function requires neldermead to be installed.")
  }
  matrix <- build_edge_optimisation_matrix(data, graph, parameters)
  full_matrix <- matrix$full
  reduced_matrix <- matrix$column_reduced
  if (length(parameters$admix_prop) == 0) {
    # I want to create "named numeric(0)" as the optimal admix vector,
    # just for the sake of consistency.
    temp <- c(1)
    names(temp) <- c(1)
    best_fit <- temp[!1]
  } else {
    x0 <- rep(0.5, length(parameters$admix_prop))
    cfunc <- cost_function(data, reduced_matrix, graph, parameters, iteration_multiplier)
    opti <- neldermead::fminbnd(cfunc, x0 = x0, xmin = rep(0, length(x0)),
                                xmax = rep(1, length(x0)),
                                options = optimisation_options)
    # The value opti is a class "neldermead" object.
    best_fit <- neldermead::neldermead.get(opti, "xopt") # Optimal admix values.
    best_fit <- best_fit[, 1]
    names(best_fit) <- parameters$admix_prop
  }
  detailed_fit <-
    edge_optimisation_function(data, full_matrix, graph, parameters, iteration_multiplier)(best_fit)
  data$graph_f4 <- detailed_fit$approximation
  # The output is a list with class "agraph_fit"
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
  cat("Call: ")
  print(x$call)
  cat("\n")
  R <- max(2^(length(x$best_fit)) - 1, 1)
  if (length(x$complaint) > 0) {
    cat("None of the admix variables are properly fitted!")
    cat("\n")
  }
  if (R %in% x$complaint) {
    cat("Not all the admix variables are properly fitted!")
    cat("\n")
    cat("See summary() for a more detailed analysis.")
    cat("\n")
    cat("\n")
  }
  cat("Minimal error:", x$best_error)
  cat("\n")
  cat("\n")
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
  cat("Call: ")
  print(object$call)
  cat("\n")
  parameters <- extract_graph_parameters(object$graph)
  R <- 2^(length(parameters$admix_prop)) - 1
  if (R != 0) {
    for (r in object$complaint) {
      if (r != R) {
        cat("After fixing ")
      }
      fixed <- ""
      complement <- ""
      for (a in seq(1, length(parameters$admix_prop))) {
        if ((r / (2^a)) %% 1 < 0.5) {
          fixed <- paste(fixed, parameters$admix_prop[a], sep = ", ")
        } else {
          complement <- paste(complement, parameters$admix_prop[a], sep = ", ")
        }
      }
      fixed <- paste("{", substring(fixed, 3), "}", sep = "")
      complement <- paste("{", substring(complement, 3), "}", sep = "")
      if (r != R) {
        cat(fixed)
        cat(" none of the remaining variables ")
      } else {
        cat("None of the variables ")
      }
      cat(complement)
      cat(" affect the quality of the fit!")
      cat("\n")
    }
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
  cat("Adding any such solution to the optimal one will not affect the error.")
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














