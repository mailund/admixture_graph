#' Used to recognize similar expressions and to possibly simplify them.
#'
#' It's best to simplify algebraic expression a little before evaluating.
#'
#' @param x  The input is assumed to be a \code{\link{character}} consisting of one or
#'           more terms. Each term starts with either \code{+} or \code{-} and after that
#'           contains one or more factors separated by \code{*}. Each factor is either
#'           an admix variable, a number or a clause \code{(1 - x)} (mind the spaces,
#'           this is how the function \code{\link{f4}} outputs), where \code{x} is again
#'           either an admix variable or a number.
#'           Everything is pretty much ruined if variable names are numbers or contain
#'           forbidden symbols \code{+, -, *, (, )}.
#'
#' @return A polynomial in a canonical form with no parenthesis or spaces and the
#'         monomials in lexicographical order. If everything is cancelled out then \code{+0}.
canonise_expression <- function(x) {
  # First remove the symbols "1" and "-" from inside the parenthesis to make things easier.
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
      coefficient <- -1
    }
    start <- 2
    for (j in seq(2, nchar(monomial))) {
      if (substring(monomial, j, j) == "*") {
        word <- substring(monomial, start, j - 1)
        if (suppressWarnings(!is.na(as.numeric(word))) == TRUE) {
          coefficient <- coefficient * suppressWarnings(as.numeric(word))
        }
        else {
          temp_vector <- c(temp_vector, word)
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
  # to return "+0" for later use (because as.numeric does not interpret "" as zero).
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
#' The elements are characters containing numbers, admix variable names,
#' parenthesis and arithmetical operations. (Transform into expressions with
#' \code{\link{parse}} and then evaluate with \code{\link{eval}}). The default
#' column names are the edge names from \code{\link{extract_graph_parameters}},
#' the rows have no names.
#'
#' @param data        The data set.
#' @param graph       The admixture graph.
#' @param parameters  In case one wants to tweak something in the graph.
#'
#' @return A list containing the full matrix (\code{full}), a version with zero
#'         columns removed (\code{column_reduced}) and parameters to pass forward
#'        (\code{parameters}).
build_edge_optimisation_matrix <- function(data, graph, parameters
                                           = extract_graph_parameters(graph)) {
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
        if (NROW(statistic[[j]]$positive) > 0) { # Insert the positive stuff.
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
        if (NROW(statistic[[j]]$negative) > 0) { # Insert the negative stuff.
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
  # Simplify by putting each element to a canonical form.
  for (i in seq(1, m)) {
    for (j in seq(1, n)) {
      edge_optimisation_matrix[i,j] <-
        canonise_expression(edge_optimisation_matrix[i, j])
    }
  }
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
  list(full = edge_optimisation_matrix, column_reduced = column_reduced,
       parameters = parameters)
}

#' Examine the edge optimisation matrix to detect unfitted admix variables.
#'
#' If the essential number of equations is not higher than the essential number of
#' edge variables, the quality of edge optimisation will not depend on the admix
#' variables (expect possibly in isolated special cases where the quality can be worse),
#' and a complaint will be given.
#' Note: The admix variable not being fitted does not mean that there is no evidence of
#' an admix event! Isolated values of the admix variables, possibly \eqn{0} or \eqn{1},
#' might give significantly worse fit than a typical value (but not the other way around).
#'
#' @param matrix  Not really a matrix but two (should be an output of
#'                \code{\link{build_edge_optimisation_matrix}}).
#' @param tol     Calulating the rank with \code{\link{qr.solve}} sometimes crashes.
#'                Default \eqn{10^{-8}}{10^(-8)}.
#'
#' @return An indicator of warning (\code{complaint}), coding all the possibilities in
#'         a way that is interpreted elsewhere (in \code{\link{summary.agraph_fit}}).
#'
#' @seealso \code{\link{qr.solve}}
examine_edge_optimisation_matrix <- function(matrix, tol = 1e-8) {
  # In order to indentify which admix variables are trurly fitted and which are
  # not, we assign a random value to some of them, not all but maybe none. We can
  # accurately detect if none of the admix variables affect the fit, but if some
  # do and some do not we have a probability zero chance of false accusation.
  # More specifically, we treat some of the admix variables as not variables at all,
  # but constants represented by a random number from the unit interval. With
  # extremely bad luck it's theoretically possible that the random constant chosen
  # decreases the rank of the edge optimisation matrix more than a typical constant
  # would.
  edge_optimisation_matrix <- matrix$full
  m <- NROW(edge_optimisation_matrix)
  n <- NCOL(edge_optimisation_matrix)
  parameters <- matrix$parameters
  complaint <- integer()
  if (length(parameters$admix_prop) > 0) {
    # The essential amount of variables is obtained by treating all the admix variables as
    # constants and calculating the rank.
    A <- rep(NaN, length(parameters$admix_prop))
    for (a in seq(1, length(parameters$admix_prop))) {
      A[a] <- stats::runif(1)
    }
    evaluated_matrix <- edge_optimisation_matrix
    for (i in seq(1, m)) {
      for (j in seq(1, n)) {
        for (a in seq(1, length(parameters$admix_prop))) {
          evaluated_matrix[i, j] <- gsub(parameters$admix_prop[a], A[a], evaluated_matrix[i, j])
        }
        evaluated_matrix[i, j] <- eval(parse(text = evaluated_matrix[i, j]))
      }
    }
    free_vars <- qr(evaluated_matrix, tol = tol)$rank
    # That was a bit of repetition but I don't want to break anything that works ;-)
    # Now we go through all the non-full different subsets of the admix variables.
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
          A[a] <- stats::runif(1) # Note that the last case r = R is not assigning anything.
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
        h <- qr(big_matrix, tol = tol)$rank
      } else {
        h <- 0
      }
      if (h <= free_vars) {
        complaint <- c(complaint, r)
        for (a in seq(1, length(parameters$admix_prop))) {
          if ((r/(2^a)) %% 1 >= 0.5) {
            skip[r - 2^(a - 1)] <- TRUE
          }
        }
      }
    }
  }
  complaint
}

#' Non negative least square solution.
#'
#' This is the function \code{\link{lsqnonneg}} from the package \code{pracma}, 
#' I just changed \code{\link{qr.solve}} into using Moore-Penrose inverse instead
#' (\code{\link{ginv}} from \code{MASS}) as \code{\link{qr.solve}} crashes for
#' some singular inputs. Now it won't crash but it's sometimes running for very long
#' time (forever?), presumably with those problematic inputs. After too many steps
#' the function halts and lies that the fit was terrible. I don't think this will
#' cause problems.
#'
#' @param C                     The matrix.
#' @param d                     The vector.
#' @param iteration_multiplier  The definition of "too many steps". Default value is \eqn{3}
#'                              (times \eqn{10} times the matrix height).
#'
#' @return A vector (\code{x}) and the error (\code{resid.norm}).
#'
#' @seealso \code{\link[pracma]{lsqnonneg}}
#' @seealso \code{\link{qr.solve}}
#' @seealso \code{\link[MASS]{ginv}}
mynonneg <- function(C, d, iteration_multiplier = 3) {
  m <- NROW(C)
  n <- NCOL(C)
  tol <- 10 * 2.220446e-16 * norm(C, "F") * (max(n, m) + 1)
  x  <- rep(0, n)             # Initial point.
  P  <- logical(n); Z <- !P   # Non-active/active columns.
  resid <- d - C %*% x
  w <- t(C) %*% resid
  wz <- numeric(n)
  # Iteration parameters.
  outeriter <- 0; it <- 0
  itmax <- 10 * iteration_multiplier * n; exitflag <- 1
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
      z[P] <- MASS::ginv(C[, P], tol = 2.22e-16) %*% d
    }
    x <- z
    resid <- d - C %*% x
    w <- t(C) %*% resid
  }
  return(list(x = x, resid.norm = sum(resid*resid)))
}

#' The cost function fed to Nelder-Mead.
#'
#' We want Nelder-Mead to run fast so the cost function operates with the column
#' reduced edge optimisation matrix and does not give any extra information about
#' the fit. For the details, use \code{\link{edge_optimisation_function}} instead.
#'
#' @param data                  The data set.
#' @param concentration         The Cholesky decomposition of the inverted covariance
#'                              matrix.
#' @param matrix                A column reduced edge optimisation matrix (typically given
#'                              by the function \code{\link{build_edge_optimisation_matrix}}).
#' @param graph                 The admixture graph.
#' @param parameters            In case one wants to tweak something in the graph.
#' @param iteration_multiplier  Given to \code{\link{mynonneg}}.
#'
#' @return  Given an input vector of admix variables, returns the smallest error
#'          regarding the edge variables.
#'
#' @seealso \code{\link{mynonneg}}
#' @seealso \code{\link{edge_optimisation_function}}
#' @seealso \code{\link{log_likelihood}}
cost_function <- function(data, concentration, matrix, graph,
                          parameters = extract_graph_parameters(graph),
                          iteration_multiplier = 3) {
  goal <- cbind(data$D)
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
    # The matrix A is taking variances and covariances into account.
    if (NCOL(matrix) > 0) {
      C <- concentration %*% evaluated_matrix
      d <- as.vector(concentration %*% goal)
      lsq_solution <- mynonneg(C, d, iteration_multiplier)
      cost <- lsq_solution$resid.norm
    } else {
      cost <- sum(as.vector(concentration %*% goal)^2)
    }
    return(cost)
  }
}

#' More detailed edge fitting than mere cost_function.
#'
#' Returning the cost, an example edge solution of an optimal fit, and linear
#' relations describing the set of all edge solutions. Operating with the full
#' edge optimisation matrix, not the column reduced one.
#'
#' @param data                  The data set.
#' @param concentration         The Cholesky decomposition of the inverted covariance matrix.
#' @param matrix                A full edge optimisation matrix (typically given by the
#'                              function \code{\link{build_edge_optimisation_matrix}}).
#' @param graph                 The admixture graph.
#' @param parameters            In case one wants to tweak something in the graph.
#' @param iteration_multiplier  Given to \code{\link{mynonneg}}.
#'
#' @return  Given an input vector of admix variables, returns a list containing
#'          the minimal error (\code{cost}), the graph-\eqn{f4} statistics
#'          (\code{approximation}), an example solution (\code{edge_fit}), linear
#'          relations describing all the solutions (\code{homogeneous}) and one
#'          way to choose the free (\code{free_edges}) and bounded
#'          (\code{bounded_edges}) edge variables.
#'
#' @seealso \code{\link{mynonneg}}
#' @seealso \code{\link{cost_function}}
#' @seealso \code{\link{log_likelihood}}
edge_optimisation_function <- function(data, concentration, matrix, graph,
                                       parameters = extract_graph_parameters(graph),
                                       iteration_multiplier = 3) {
  goal <- cbind(data$D)
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
    C <- concentration %*% evaluated_matrix
    d <- as.vector(concentration %*% goal)
    lsq_solution <- mynonneg(C, d, iteration_multiplier)
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

#' Building a proxy concentration matrix.
#' 
#' If we don't have the true concentration matrix of the data rows calculated,
#' but at least have the \eqn{Z} scores of individual rows, (unrealistically) assuming
#' independence and calculating the concentration matrix from those is still better
#' than nothing (\emph{i. e.} the identity matrix).
#' 
#' @param data     The data containing at least the expected values of \eqn{f} statistics
#'                 (column \code{D}) and possibly also products of expected values and \eqn{f}
#'                 statistics divided by standard deviations of (the \eqn{Z} scores,
#'                 column \code{Z.value}).
#' @param Z.value  Tells whether the \eqn{Z} scores are available or should we just use the
#'                 identity matrix.
#' 
#' @return The Cholesky decomposition of the inverted covariance matrix.
calculate_concentration <- function(data, Z.value) {
  concentration <- matrix(0, NROW(data), NROW(data))
  if (Z.value == TRUE) {
    for (j in seq(1, NROW(data))) {
      concentration[j, j] <-  as.numeric(data[j, "Z.value"])/as.numeric(data[j, "D"])
    }
  }
  else {
    for (j in seq(1, NROW(data))) {
      concentration[j, j] <-  1
    }
  }
  concentration
}

#' A fast version of graph fitting.
#'
#' Given a table of observed \eqn{f} statistics and a graph, uses Nelder-Mead algorithm to
#' find the graph parameters (edge lengths and admixture proportions) that minimize the value
#' of \code{\link{cost_function}}, \emph{i. e.} maximizes the likelihood of a graph with
#' parameters given the observed data.
#' Like \code{\link{fit_graph}} but dropping most of the analysis on the result.
#' Intended for use in big iteration loops.
#'
#' @param data                  The data table, must contain columns \code{W}, \code{X},
#'                              \code{Y}, \code{Z} for sample names and \code{D} for the
#'                              observed \eqn{f_4(W, X; Y, Z)}. May contain an optional
#'                              column \code{Z.value} for the \eqn{Z} scores (the \eqn{f}
#'                              statistics divided by the standard deviations).
#' @param graph                 The admixture graph (an \code{\link{agraph}} object).
#' @param point                 If the user wants to restrict the admixture proportions somehow,
#'                              like to fix some of them. A list of two vectors: the lower and the
#'                              upper bounds. As a default the bounds are just it little bit more
#'                              than zero and less than one; this is because sometimes the
#'                              infimum of the values of cost function is at a point of
#'                              non-continuity, and zero and one have reasons to be problematic
#'                              values in this respect.
#' @param Z.value               Whether we calculate the default concentration from \eqn{Z} scores
#'                              (the default option \code{TRUE}) or just use the identity matrix.
#' @param concentration         The Cholesky decomposition of the inverted covariance matrix.
#'                              Default matrix determined by the parameter \code{Z.value}.
#' @param optimisation_options  Options to the Nelder-Mead algorithm.
#' @param parameters            In case one wants to tweak something in the graph.
#' @param iteration_multiplier  Given to \code{\link{mynonneg}}.
#'
#' @return A list containing only the essentials about the fit:
#'         \code{graph} is the graph input,
#'         \code{best_error} is the minimal value of \code{\link{cost_function}},
#'         obtained when the admixture proportions are \code{best_fit}.
#'
#' @seealso \code{\link{cost_function}}
#' @seealso \code{\link{agraph}}
#' @seealso \code{\link{calculate_concentration}}
#' @seealso \code{\link[neldermead]{optimset}}
#' @seealso \code{\link{fit_graph}}
#'
#' @examples
#' \donttest{
#' # For example, let's fit the following two admixture graph to an example data on bears:
#' 
#' data(bears)
#' print(bears)
#' 
#' leaves <- c("BLK", "PB", "Bar", "Chi1", "Chi2", "Adm1", "Adm2", "Denali", "Kenai", "Sweden") 
#' inner_nodes <- c("R", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z", "M", "N")
#' edges <- parent_edges(c(edge("BLK", "R"),
#'                         edge("PB", "v"),
#'                         edge("Bar", "x"),
#'                         edge("Chi1", "y"),
#'                         edge("Chi2", "y"),
#'                         edge("Adm1", "z"),
#'                         edge("Adm2", "z"),
#'                         edge("Denali", "t"),
#'                         edge("Kenai", "s"),
#'                         edge("Sweden", "r"),
#'                         edge("q", "R"),
#'                         edge("r", "q"),
#'                         edge("s", "r"),
#'                         edge("t", "s"),
#'                         edge("u", "q"),
#'                         edge("v", "u"),
#'                         edge("w", "M"),
#'                         edge("x", "N"),
#'                         edge("y", "x"),
#'                         edge("z", "w"),
#'                         admixture_edge("M", "u", "t"),
#'                         admixture_edge("N", "v", "w")))
#' admixtures <- admixture_proportions(c(admix_props("M", "u", "t", "a"),
#'                                       admix_props("N", "v", "w", "b")))
#' bears_graph <- agraph(leaves, inner_nodes, edges, admixtures)
#' plot(bears_graph, show_admixture_labels = TRUE)
#'
#' fit <- fast_fit(bears, bears_graph)
#' print(fit$best_error)
#' 
#' # The result is just the minimal value of the cost function and the values of admixture proportions
#' # where it's obtained, no deeper analysis of the fit.
#' }
#' 
#' @export
fast_fit <- function(data, graph,
                     point = list(rep(1e-5, length(extract_graph_parameters(graph)$admix_prop)),
                                  rep(1 - 1e-5, length(extract_graph_parameters(graph)$admix_prop))),
                     Z.value = TRUE,
                     concentration = calculate_concentration(data, Z.value),
                     optimisation_options = NULL,
                     parameters = extract_graph_parameters(graph),
                     iteration_multiplier = 3) {
  withCallingHandlers({
    inner_fast_fit(data, graph, point, Z.value, concentration, optimisation_options,
                   parameters, iteration_multiplier)
  }, error = function(e) {
    message("Something went wrong, trying again.")
    invokeRestart("try_again")
  })
}
inner_fast_fit <- function(data, graph, point, Z.value, concentration, optimisation_options,
                           parameters, iteration_multiplier) {
  withRestarts({
    matrix <- build_edge_optimisation_matrix(data, graph, parameters)
    reduced_matrix <- matrix$column_reduced
    if (length(parameters$admix_prop) == 0) {
      # I want to create "named numeric(0)" as the optimal admix vector,
      # just for the sake of consistency.
      temp <- c(1)
      names(temp) <- c(1)
      best_fit <- temp[!1]
      full_matrix <- matrix$full
      detailed_fit <-
        edge_optimisation_function(data, concentration, full_matrix, graph, parameters,
                                   iteration_multiplier)(best_fit)
      best_error <- detailed_fit$cost
    } else {
      x0 <- 0.5*(point[[1]] + point[[2]])
      cfunc <- cost_function(data, concentration, reduced_matrix, graph, parameters, iteration_multiplier)
      opti <- neldermead::fminbnd(cfunc, x0 = x0, xmin = point[[1]], xmax = point[[2]],
                                  options = optimisation_options)
      # The value opti is a class "neldermead" object.
      best_error <- neldermead::neldermead.get(opti, "fopt") # Optimal error.
      best_fit <- neldermead::neldermead.get(opti, "xopt") # Optimal admix values.
      best_fit <- best_fit[, 1]
      names(best_fit) <- parameters$admix_prop
    }
    # The output is a list containing the square sum error, admix variable fit and the graph.
    list(best_error = best_error, best_fit = best_fit, graph = graph)
  }, try_again = function() {
    inner_fast_fit(data, graph, point, Z.value, concentration, optimisation_options, parameters,
                   iteration_multiplier)
  })
}

#' Fit the graph parameters to a data set.
#'
#' Given a table of observed \eqn{f} statistics and a graph, uses Nelder-Mead algorithm to
#' find the graph parameters (edge lengths and admixture proportions) that minimize the value
#' of \code{\link{cost_function}}, \emph{i. e.} maximizes the likelihood of a graph with
#' parameters given the observed data.
#' Like \code{\link{fast_fit}} but outputs a more detailed analysis on the results.
#'
#' @param data                  The data table, must contain columns \code{W}, \code{X},
#'                              \code{Y}, \code{Z} for sample names and \code{D} for the
#'                              observed \eqn{f_4(W, X; Y, Z)}. May contain an optional
#'                              column \code{Z.value} for the \eqn{Z} scores (the \eqn{f}
#'                              statistics divided by the standard deviations).
#' @param graph                 The admixture graph (an \code{\link{agraph}} object).
#' @param point                 If the user wants to restrict the admixture proportions somehow,
#'                              like to fix some of them. A list of two vectors: the lower and the
#'                              upper bounds. As a default the bounds are just it little bit more
#'                              than zero and less than one; this is because sometimes the
#'                              infimum of the values of cost function is at a point of
#'                              non-continuity, and zero and one have reasons to be problematic
#'                              values in this respect.
#' @param Z.value               Whether we calculate the default concentration from \eqn{Z} scores
#'                              (the default option \code{TRUE}) or just use the identity matrix.
#' @param concentration         The Cholesky decomposition of the inverted covariance matrix.
#'                              Default matrix determined by the parameter \code{Z.value}.
#' @param optimisation_options  Options to the Nelder-Mead algorithm.
#' @param parameters            In case one wants to tweak something in the graph.
#' @param iteration_multiplier  Given to \code{\link{mynonneg}}.
#' @param qr_tol                Given to \code{\link{examine_edge_optimisation_matrix}}.
#'
#' @return A class \code{agraph_fit} list containing a lot of information about the fit: \cr
#'         \code{data} is the input data, \cr
#'         \code{graph} is the input graph, \cr
#'         \code{matrix} is the output of \code{\link{build_edge_optimisation_matrix}},
#'         containing the \code{full} matrix, the \code{column_reduced} matrix without zero
#'         columns, and graph \code{parameters}, \cr
#'         \code{complaint} coding wchich subsets of admixture proportions are trurly fitted, \cr
#'         \code{best_fit} is the optimal admixture proportions (might not be unique if they
#'         are not trurly fitted), \cr
#'         \code{best_edge_fit} is an example of optimal edge lengths, \cr
#'         \code{homogeneous} is the reduced row echelon form of the matrix describing when
#'         a vector of edge lengths have no effect on the prediced statistics \eqn{F}, \cr
#'         \code{free_edges} is one way to choose a subset of edge lengths in such a vector as
#'         free variables, \cr
#'         \code{bounded_edges} is how we calculate the reamining edge lengths from the free ones, \cr
#'         \code{best_error} is the minimum value of the \code{\link{cost_function}}, \cr
#'         \code{approximation} is the predicted statistics \eqn{F} with the optimal graph parameters, \cr
#'         \code{parameters} is jsut a shortcut for the graph parameters. \cr
#'         See \code{\link{summary.agraph_fit}} for the interpretation of some of these results.
#'
#' @seealso \code{\link{cost_function}}
#' @seealso \code{\link{agraph}}
#' @seealso \code{\link{calculate_concentration}}
#' @seealso \code{\link[neldermead]{optimset}}
#' @seealso \code{\link{fast_fit}}
#'
#' @examples
#' \donttest{
#' # For example, let's fit the following two admixture graph to an example data on bears:
#' 
#' data(bears)
#' print(bears)
#' 
#' leaves <- c("BLK", "PB", "Bar", "Chi1", "Chi2", "Adm1", "Adm2", "Denali", "Kenai", "Sweden") 
#' inner_nodes <- c("R", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z", "M", "N")
#' edges <- parent_edges(c(edge("BLK", "R"),
#'                         edge("PB", "v"),
#'                         edge("Bar", "x"),
#'                         edge("Chi1", "y"),
#'                         edge("Chi2", "y"),
#'                         edge("Adm1", "z"),
#'                         edge("Adm2", "z"),
#'                         edge("Denali", "t"),
#'                         edge("Kenai", "s"),
#'                         edge("Sweden", "r"),
#'                         edge("q", "R"),
#'                         edge("r", "q"),
#'                         edge("s", "r"),
#'                         edge("t", "s"),
#'                         edge("u", "q"),
#'                         edge("v", "u"),
#'                         edge("w", "M"),
#'                         edge("x", "N"),
#'                         edge("y", "x"),
#'                         edge("z", "w"),
#'                         admixture_edge("M", "u", "t"),
#'                         admixture_edge("N", "v", "w")))
#' admixtures <- admixture_proportions(c(admix_props("M", "u", "t", "a"),
#'                                       admix_props("N", "v", "w", "b")))
#' bears_graph <- agraph(leaves, inner_nodes, edges, admixtures)
#' plot(bears_graph, show_admixture_labels = TRUE)
#'
#' fit <- fit_graph(bears, bears_graph)
#' summary(fit)
#' 
#' # It turned out the values of admixture proportions had no effect on the cost function. This is not
#' # too surprising because the huge graph contains a lot of edge variables compared to the tiny 
#' # amount of data we used! Note however that the mere existence of the admixture event with non- 
#' # trivial (not zero or one) admixture proportion might still decrease the cost function.
#' }
#'
#' @export
fit_graph <- function(data, graph,
                      point = list(rep(1e-5, length(extract_graph_parameters(graph)$admix_prop)),
                                   rep(1 - 1e-5, length(extract_graph_parameters(graph)$admix_prop))),
                      Z.value = TRUE,
                      concentration = calculate_concentration(data, Z.value),
                      optimisation_options = NULL,
                      parameters = extract_graph_parameters(graph), 
                      iteration_multiplier = 3, qr_tol = 1e-8) {
  withCallingHandlers({
    inner_fit_graph(data, graph, point, Z.value, concentration, optimisation_options,
                    parameters, iteration_multiplier, qr_tol)
  }, error = function(e) {
    message("Something went wrong, trying again.")
    invokeRestart("try_again")
  })
}
inner_fit_graph <- function(data, graph, point, Z.value, concentration, optimisation_options,
                            parameters, iteration_multiplier, qr_tol) {
  withRestarts({
    matrix <- build_edge_optimisation_matrix(data, graph, parameters)
    full_matrix <- matrix$full
    reduced_matrix <- matrix$column_reduced
    complaint <- examine_edge_optimisation_matrix(matrix, qr_tol)
    if (length(parameters$admix_prop) == 0) {
      # I want to create "named numeric(0)" as the optimal admix vector,
      # just for the sake of consistency.
      temp <- c(1)
      names(temp) <- c(1)
      best_fit <- temp[!1]
    } else {
      x0 <- 0.5*(point[[1]] + point[[2]])
      cfunc <- cost_function(data, concentration, reduced_matrix, graph, parameters, iteration_multiplier)
      opti <- neldermead::fminbnd(cfunc, x0 = x0, xmin = point[[1]], xmax = point[[2]],
                                  options = optimisation_options)
      # The value opti is a class "neldermead" object.
      best_fit <- neldermead::neldermead.get(opti, "xopt") # Optimal admix values.
      best_fit <- best_fit[, 1]
      names(best_fit) <- parameters$admix_prop
    }
    detailed_fit <-
      edge_optimisation_function(data, concentration, full_matrix, graph, parameters,
                                 iteration_multiplier)(best_fit)
    data$graph_f4 <- detailed_fit$approximation
    # The output is a list with class "agraph_fit"
    structure(list(
      call = sys.call(),
      data = data,
      graph = graph,
      matrix = matrix,
      complaint = complaint,
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
  }, try_again = function() {
    inner_fit_graph(data, graph, point, Z.value, concentration, optimisation_options, parameters,
                    iteration_multiplier, qr_tol)
  })
}

#' Calculate (essentially) the log likelihood of a graph with parameters, given the observation.
#'
#' Or the log likelihood of the observation, given graph with parameters, depending how things are modeled.
#' Basically this is just \code{\link{cost_function}} that doesn't optimize the edge variables
#' but has them as an argument instead.
#' 
#' @param f              The observed \eqn{f} statistics (the column \code{D} from \code{data}).
#' @param concentration  The Cholesky decomposition of the inverted covariance matrix. So if \eqn{S}
#'                       is the covariance matrix, this is \eqn{C = chol(S^{-1})}{C = chol(S^(-1))} satisfying
#'                       \eqn{S^{-1} = C^t C}{S^(-1) = C^t*C}.
#' @param matrix         A column reduced edge optimisation matrix (typically given by the function
#'                       \code{\link{build_edge_optimisation_matrix}}).
#' @param graph          The admixture graph. Here to give default value for:
#' @param parameters     Just because we need to know variable names.
#'
#' @return The output is a function. Given admixture proportions \code{x} and edge lengths \code{e}, the graph
#'         topology \eqn{T} implies an estimate \eqn{F} for the statistics \eqn{f}. This output function
#'         calculates
#'         \deqn{l = (F-f)^t S^{-1}(F-f)}{l = (F-f)^t*S^(-1)*(F-f)}
#'         from \code{x} and \code{e}. Up to a constant error and multiplier that is a log likelihood function, as
#'         \deqn{\det(2 \pi S)^{-1/2} e^{-l/2}}{det(2*\pi*S)^(-1/2)*exp(-l/2)}
#'         can be seen as a likelihood of a graph with parameters, given the observation, or the other way around
#'         (possibly up to a normalization constant).
#'
#' @seealso \code{\link{cost_function}}
#' @seealso \code{\link{edge_optimisation_function}}
#' @seealso \code{\link{calculate_concentration}}
#'
#' @export
log_likelihood <- function(f, concentration, matrix, graph, parameters = extract_graph_parameters(graph)) {
  function(x, e) {
    # Evaluate the column reduced edge optimisation matrix at admix values x.
    evaluated_matrix <- matrix(0, NROW(matrix), NCOL(matrix))
    if (length(parameters$admix_prop) > 0) {
      for (i in seq(1, length(parameters$admix_prop))) {
        assign(parameters$admix_prop[i], x[i])
      }
    }
    if (NCOL(matrix) > 0) {
      for (i in seq(1, NROW(matrix))) {
        for (j in seq(1, NCOL(matrix))) {
          evaluated_matrix[i, j] <- eval(parse(text = matrix[i, j]))
        }
      }
    }
    # Now just multiply with edge values e to get F.
    if (NCOL(matrix) > 0) {
      FF <- evaluated_matrix %*% e
      vector <- concentration %*% (FF - f)
    } else {
      vector <- concentration %*% f
    }
    likelihood <- -0.5 * sum(as.vector(vector)^2)
    return(likelihood)
  }
}

#' Print function for the fitted graph.
#'
#' Prints the value of \code{\link{cost_function}} of the fitted graph, and complains
#' if some or all of the admixture proportions aren't trurly fitted.
#' Note: the admixture proportion not being trurly fitted does not necessarily mean that
#' there is no evidence of an admix event!
#'
#' @param x    The fitted object.
#' @param ...  Additional parameters.
#'
#' @seealso \code{link{summary.agraph_fit}}
#'
#' @export
print.agraph_fit <- function(x, ...) {
  cat("\n")
  cat("Call: ")
  print(x$call)
  cat("\n")
  R <- max(2^(length(x$best_fit)) - 1, 1)
  if (length(x$complaint) > 0) {
    cat("None of the admixture proportions are properly fitted!")
    cat("\n")
  }
  if (R %in% x$complaint) {
    cat("Not all of the admixture proportions are properly fitted!")
    cat("\n")
    cat("See summary.agraph_fit for a more detailed analysis.")
    cat("\n")
    cat("\n")
  }
  cat("Minimal error:", x$best_error)
  cat("\n")
  cat("\n")
}

#' Parameters for the fitted graph.
#'
#' Extracts the graph parameters for the graph fitted to data. Note that the optimal 
#' parameters are generally not unique.
#'
#' @param object  The fitted object.
#' @param ...     Additional parameters.
#' 
#' @seealso \code{link{summary.agraph_fit}}
#'
#' @export
coef.agraph_fit <- function(object, ...) {
  c(object$best_edge_fit, object$best_fit)
}

#' Summary for the fitted graph.
#'
#' Prints: \cr
#' Optimal admixture proportions and a complaint if some of them are not trurly fitted,
#' \emph{i. e.} if after fixing a (possibly empty) subset of them, the rest have typically
#' no effect on the cost function. Here typically means that some isolated values of the 
#' admixture proportions, like \eqn{0} or \eqn{1}, might actually give a significantly worse
#' fit than the constant fit given by any other values (but not better).
#' Thus, an admixture proportion not being fitted does not always mean that there is no
#' evidence of an admix event, as fixing them at \eqn{0} or \eqn{1} could make the fit
#' worse while the exact value won't matter otherwise. \cr
#' The optimal edge lengths give one of the solutions for the best fit. It is generally not
#' unique, as after fixing the admixture proportions, the best edge lengths are a non-negative
#' least square solution for a system of linear equations. To get all the solutions one has
#' to add any solution of the corresponding homogeneous system to the given exaple solution
#' (and exclude possible negative values). The solutions of the homogeneous system are given
#' as a set of free edge lengths that may obtain any non-negative value, and bounded edge
#' lengths that linearly depend on the free ones. \cr
#' Minimal error is the value of the \code{\link{cost_function}} at the fit.
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
        cat(" none of the remaining proportions ")
      } else {
        cat("None of the proportions ")
      }
      cat(complement)
      cat(" affect the quality of the fit!")
      cat("\n")
    }
  }
  cat("\n")
  cat("Optimal admixture proportions:")
  cat("\n")
  print(object$best_fit)
  cat("\n")
  cat("Optimal edge lengths:")
  cat("\n")
  print(object$best_edge_fit)
  cat("\n")
  cat("Solution to a homogeneous system of edge lengths with the optimal admixture proportions:")
  cat("\n")
  cat("Adding any such solution to the optimal one will not affect the error.")
  cat("\n")
  cat("\n")
  cat("Free edge lengths:")
  cat("\n")
  cat(object$free_edges, sep = "\n")
  cat("\n")
  cat("Bounded edge lengths:")
  cat("\n")
  cat(object$bounded_edges, sep = "\n")
  cat("\n")
  cat("Minimal error:")
  cat("\n")
  cat(object$best_error)
}

#' Predicted f statistics for the fitted graph.
#'
#' Gets the predicted \eqn{f} statistics \eqn{F} for the fitted graph.
#'
#' @param object  The fitted object.
#' @param ...     Additional parameters.
#' 
#' @seealso \code{link{summary.agraph_fit}}
#'
#' @export
fitted.agraph_fit <- function(object, ...) {
  object$data
}

#' Errors of prediction in the fitted graph
#'
#' Gets \eqn{f - F}, the difference between predicted and observed  statistics,
#' for each data point used in the fit.
#'
#' @param object  The fitted object.
#' @param ...     Additional parameters.
#' 
#' @seealso \code{link{summary.agraph_fit}}
#'
#' @export
residuals.agraph_fit <- function(object, ...) {
  object$data$D - object$data$graph_f4
}
