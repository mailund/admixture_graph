#' Plot an admixture graph.
#' 
#' This is a basic drawing routine for visualising the graph. For publication
#' quality graphs a lot more tweaking is probably needed.
#' 
#' @param x The admixture graph.
#' @param ordered_leaves The leaf-nodes in the left to right order they should
#'   be drawn. I don't have a good algorithm for figuring out that order so for
#'   now it is required as a function argument.
#' @param show_admixture_labels A flag determining if the plot should include
#'   the names of admixture proportions.
#' @param show_inner_node_labels A flat determining if the plot should include
#'   the names of inner nodes.
#'   
#' @param ... Additional plotting options
#'   
#' @export
plot.agraph <- function(x, 
                        ordered_leaves = NULL,
                        show_admixture_labels = FALSE,
                        show_inner_node_labels = FALSE,
                        ...) {
  
  graph <- x
  
  if (is.null(ordered_leaves))
    ordered_leaves <- graph$leaves
  
  dfs <- function(node, basis, step) {
    result <- rep(NA, length(graph$nodes))
    names(result) <- graph$nodes
    
    dfs_ <- function(node) {
      children <- which(graph$children[node,])
      if (length(children) == 0) {
        result[node] <<- basis(node)
      } else {
        result[node] <<- step(vapply(children, dfs_, numeric(1)))
      }
    }
    dfs_(node)
    result
  }
  
  no_parents <- function(node) length(which(graph$parents[node, ]))
  roots <- which(Map(no_parents, graph$nodes) == 0)
  if (length(roots) > 1) stop("Don't know how to handle more than one root")
  root <- roots[1]
  ypos <- dfs(root, basis = function(x) 0.0, step = function(x) max(x) + 1.0)
  
  leaf_index <- function(n) {
    result <- which(graph$nodes[n] == ordered_leaves)
    if (length(result) != 1) stop("Unexpected number of matching nodes")
    result
  }
  left_x  <- dfs(root, basis = leaf_index, step = min)
  right_x <- dfs(root, basis = leaf_index, step = max)
  xpos <- left_x + (right_x - left_x) / 2.0
  
  # Start the actual drawing of the graph...
  plot(xpos, ypos, type = "n", axes = FALSE, frame.plot = FALSE,
       xlab = "", ylab = "", ylim = c(-1, max(ypos) + 0.5), ...)
  
  for (node in graph$nodes) {
    parents <- graph$nodes[graph$parents[node, ]]
    if (length(parents) == 1) {
      lines(c(xpos[node],xpos[parents]), c(ypos[node], ypos[parents]))
      
    } else if (length(parents) == 2) {
      break_y <- ypos[node]
      break_x_left <- xpos[node] - 0.3
      break_x_right <- xpos[node] + 0.3
      
      if (xpos[parents[1]] < xpos[parents[2]]) {
        lines(c(xpos[parents[1]], break_x_left), c(ypos[parents[1]], break_y))
        lines(c(xpos[parents[2]], break_x_right), c(ypos[parents[2]], break_y))  
      } else {
        lines(c(xpos[parents[2]], break_x_left), c(ypos[parents[2]], break_y))
        lines(c(xpos[parents[1]], break_x_right), c(ypos[parents[1]], break_y))
      }
      
      
      segments(break_x_left, break_y, xpos[node], ypos[node], col = "red")
      segments(break_x_right, break_y, xpos[node], ypos[node], col = "red")
      
      if (show_admixture_labels) {
        if (xpos[parents[1]] < xpos[parents[2]]) {
          text(break_x_left, break_y, graph$probs[parents[[1]], node],
               cex = 0.5, pos = 1, col = "red", offset = 0.1)
          text(break_x_right, break_y, graph$probs[parents[[2]], node],
               cex = 0.5, pos = 1, col = "red", offset = 0.1)
        } else {
          text(break_x_left, break_y, graph$probs[parents[[2]], node],
               cex = 0.5, pos = 1, col = "red", offset = 0.1)
          text(break_x_right, break_y, graph$probs[parents[[1]], node],
               cex = 0.5, pos = 1, col = "red", offset = 0.1)          
        }
      }
    }
  }
  
  is_inner <- Vectorize(function(n) sum(graph$children[n, ]) > 0)
  inner_nodes <- which(is_inner(graph$nodes))
  leaves <- which(!is_inner(graph$nodes))
  
  if (show_inner_node_labels) {
    text(xpos[inner_nodes], ypos[inner_nodes], 
         labels = graph$nodes[inner_nodes], cex = 0.6, col = "blue", pos = 3)
  }
  text(xpos[leaves], ypos[leaves], labels = graph$nodes[leaves], 
       cex = 0.7, col = "black", pos = 1)
  
  invisible()
}

# Determining the crossing numbers of a graph is a difficult problem, and requiring
# edges to be straight lines makes it even worse. The following algorithm is by no
# means meant to be optimal. We only perform some clean-up on the cycles of the graph,
# using a heuristic of dealing with the worst cycles first.

# We already have the graph environment but we will need some more details.
# Some stuff we just borrow straight from the agraph envionment, repeated here for
# consistency's sake.

refined_graph <- function(graph) {
  leaves <- graph$leaves
  inner_nodes <- graph$inner_nodes
  nodes <- graph$nodes
  probs <- graph$probs
  parents <- graph$parents
  parent_of <- rep(NA, length(nodes))
  names(parent_of) <- nodes
  children <- graph$children
  child_order <- list()
  admix_nodes <- character()
  parent_order <- list()
  for (i in seq(1, length(nodes))) {
    child_list <- character()
    parent_list <- character()
    parent_count <- 0
    for (j in seq(1, length(nodes))) {
      if (children[i, j] == TRUE) {
        child_list <- c(child_list, nodes[j])
      }
      if (parents[i, j] == TRUE) {
        parent_of[i] <- nodes[j]
        parent_count <- parent_count + 1
        parent_list <- c(parent_list, nodes[j])
      }
    }
    child_order[[i]] <- child_list
    names(child_order)[i] <- nodes[i]
    if (parent_count == 0) {
      root <- nodes[i]
    }
    if (parent_count == 1) {
      parent_list <- character()
    }
    if (parent_count > 1) {
      parent_of[i] <- NA
      admix_nodes <- c(admix_nodes, nodes[i])
    }
    parent_order[[i]] <- parent_list
    names(parent_order)[i] <- nodes[i]
  }
  structure(list(nodes = nodes,
                 leaves = leaves,
                 inner_nodes = inner_nodes,
                 admix_nodes = admix_nodes,
                 root = root,
                 probs = probs,
                 parents = parents,
                 parent_of = parent_of,
                 parent_order = parent_order,
                 children = children,
                 child_order = child_order),
            class = "refined_graph")
}

#' @export
arrange_graph <- function(graph) {
  graph <- refined_graph(graph)
  print("INPUT ORDER:")
  print(graph$child_order)
  print(graph$parent_order)
  # First of all let's make a convention to only call this function if there is one ore more
  # admix events, trees can always be drawn no matter what.
  # Determine cycles corresponding to admix nodes.  We follow both branches until a common
  # acestor is found. Unfortunately this is not unique if we hit more admix nodes on the
  # way. Then our first tie-breaker is the age of the ancestor, if one ancestor is an
  # ancestor of the other, choose the younger one, and the second tie-breaker is the branch
  # count of the cycle, the number of branches deviating from it, up or down, not counting
  # the admix node itself or the ancestor (doesn' matter when breaking ties of course but
  # we will use the branch count later also).
  cycles <- list()
  for (mix in seq(1, length(graph$admix_nodes))) {
    mix <- graph$admix_nodes[mix]
    cycle_candidates <- list()
    for (i in seq(1, length(all_paths_to_root(graph, mix)) - 1)) {
      for (j in seq(i + 1, length(all_paths_to_root(graph, mix)))) {
        # This is ok because mix is an admix variable and so there is at least 2 paths to root.
        everything_OK <- TRUE
        v1 <- all_paths_to_root(graph, mix)[[i]][-1]
        v2 <- all_paths_to_root(graph, mix)[[j]][-1]
        if (v1[1] == v2[1]) {
          everything_OK <- FALSE
        }
        collision <- first_collision(v1, v2)
        # Truncate for convenience.
        v1 <- v1[1:match(collision, v1) - 1]
        v2 <- v2[1:match(collision, v2) - 1]
        branch_count <- branch_count(graph, v1, v2)
        if (everything_OK == TRUE) {
          temp <- list()
          temp[[all_paths_to_root(graph, mix)[[i]][2]]] <- v1
          temp[[all_paths_to_root(graph, mix)[[j]][2]]] <- v2
          temp[["collision"]] <- collision
          temp[["branch_count"]] <- branch_count
          cycle_candidates[[length(cycle_candidates) + 1]] <- temp
        }
      }
    }
    # Now we choose the cycle candidate that has the least branch count among those cycles whose
    # collision node is not a proper ancestor of the collision node of another candidate.
    least_branch_count <- Inf
    for (i in seq(1, length(cycle_candidates))) {
      could_this_be_it <- TRUE
      for (j in seq(1, length(cycle_candidates))) {
        candidate <- cycle_candidates[[i]]$collision
        comparison <- cycle_candidates[[j]]$collision
        if (is_descendant_of(graph, comparison, candidate) == TRUE) {
          could_this_be_it <- FALSE
        }
        if (cycle_candidates[[i]]$branch_count < least_branch_count && could_this_be_it == TRUE) {
          cycles[[mix]] <- cycle_candidates[[i]]
          least_branch_count <- cycle_candidates[[i]]$branch_count
        }
      }
    }
  }
  # Before clearing the cycles we need to agree on the order on which to clear them, which is quite
  # relevant. As clearing also might change the parent ordering of an older admixture node, we should
  # at least make sure the cycles originating from younger admixture nodes are cleared first. It also
  # looks like the cycles with large branch count should be cleared first. Again, this does not remove
  # all the issues.
  cleaning_order <- character()
  # First a topological sort on the partial order of ancestry among admix nodes.
  arrow_list <- list()
  S <- graph$admix_nodes
  for (i in seq(1, length(graph$admix_nodes))) {
    for (j in seq(1, length(graph$admix_nodes))) {
      if (is_descendant_of(graph, graph$admix_nodes[i], graph$admix_nodes[j]) == TRUE) {
        key <- paste(graph$admix_nodes[i], graph$admix_nodes[j], sep = "_")
        arrow_list[[key]] <- c(graph$admix_nodes[i], graph$admix_nodes[j])
        if (is.na(match(graph$admix_nodes[j], S)) == FALSE) {
          S <- S[-match(graph$admix_nodes[j], S)]
        }
      }
    }
  }
  # S is non-empty.
  while (length(S) > 0) {
    n <- S[1]
    S <- S[-1]
    cleaning_order <- c(cleaning_order, n)
    for (i in arrow_list) {
      if (i[1] == n) {
        key <- paste(i[1], i[2], sep = "_")
        arrow_list[[key]] <- NULL # I don't want to talk about it.
        should_add <- TRUE
        for (j in arrow_list) {
          if (j[2] == i[2]) {
            should_add <- FALSE
          }
        }
        if (should_add == TRUE) {
          S <- c(S, i[2])
        }
      }
    }
  }
  # Then the admix nodes with strictly larger brach count cut in line if not forbidden by ancestry.
  if (length(graph$admix_nodes) > 1) {
    for (i in seq(2, length(graph$admix_nodes))) {
      j <- i
      while (j > 1 && is_descendant_of(graph, graph$admix_nodes[j - 1], graph$admix_nodes[j]) == FALSE
             && cycles[[graph$admix_nodes[j]]]$branch_count > cycles[[graph$admix_nodes[j - 1]]]$branch_count) {
        tempo <- cleaning_order[j - 1]
        cleaning_order[j - 1] <- cleaning_order[j]
        cleaning_order[j] <- tempo
        j <- j - 1
      }
    }
  }
  # Then accorning to the cleaning order, turn all branches away from the cycle, the orientation of which is
  # determined by the parents of the admix node.
  for (mix in cleaning_order) { # Hey why aren't all my loops like this?
    # First take the collision node and place the edges that belong to the cycle next to each other, in the
    # correct order, without skewing much to the left or right. Then clean the branches outside the cycle,
    # making an extra effort to make sure any encountered admix node with the common collision node will
    # be on the correct side of the cycle.
    left_start <- graph$parent_order[[mix]][1]
    right_start <- graph$parent_order[[mix]][2]
    if (length(cycles[[mix]][[left_start]]) > 0) {
      left_last <- cycles[[mix]][[left_start]][length(cycles[[mix]][[left_start]])]
    } else {
      left_last <- mix
    }
    if (length(cycles[[mix]][[right_start]]) > 0) {
      right_last <- cycles[[mix]][[right_start]][length(cycles[[mix]][[right_start]])]
    } else {
      right_last <- mix
    }
    left_index <- match(left_last, graph$child_order[[cycles[[mix]]$collision]])
    right_index <- match(right_last, graph$child_order[[cycles[[mix]]$collision]])
    u <- right_index - left_index - sign(right_index - left_index)
    left_move <- sign(u)*floor(0.5*abs(u))
    right_move <- -sign(u)*ceiling(0.5*abs(u)) + (1 - sign(right_index - left_index))/2
    graph$child_order[[cycles[[mix]]$collision]] <- move(graph$child_order[[cycles[[mix]]$collision]],
                                                         left_index, left_move)
    graph$child_order[[cycles[[mix]]$collision]] <- move(graph$child_order[[cycles[[mix]]$collision]],
                                                         right_index, right_move)
    for (l in cycles[[mix]][[left_start]]) {
      if (l == left_start) {
        previous <- mix
      } else {
        l_index <- match(l, cycles[[mix]][[left_start]])
        previous <- cycles[[mix]][[left_start]][l_index - 1]
      }
      graph$child_order[[l]] <- move(graph$child_order[[l]],
                                     match(previous, graph$child_order[[l]]),
                                     length(graph$child_order[[l]]) - match(previous, graph$child_order[[l]]))
      if (l %in% graph$admix_nodes) {
        if (graph$parent_order[[l]][1] %in% cycles[[mix]][[left_start]] ||
            graph$parent_order[[l]][1] == cycles[[mix]]$collision) {
          graph$parent_order[[l]] <- move(graph$parent_order[[l]], 1, 1)
        }
        if (cycles[[mix]]$collision == cycles[[l]]$collision) {
          for (v in graph$child_order[[cycles[[mix]]$collision]]) {
            if (v %in% cycles[[l]][[graph$parent_order[[l]][1]]]) {
              target <- match(left_last, graph$child_order[[cycles[[mix]]$collision]])
              own_position <- match(v, graph$child_order[[cycles[[mix]]$collision]])
              graph$child_order[[cycles[[mix]]$collision]] <-
                move(graph$child_order[[cycles[[mix]]$collision]], own_position, target - own_position - 1)
            }
          }
        }
      }
    } 
    for (r in cycles[[mix]][[right_start]]) {
      if (r == right_start) {
        previous <- mix
      } else {
        r_index <- match(r, cycles[[mix]][[right_start]])
        previous <- cycles[[mix]][[right_start]][r_index - 1]
      }
      graph$child_order[[r]] <- move(graph$child_order[[r]],
                                     match(previous, graph$child_order[[r]]),
                                     -(match(previous, graph$child_order[[r]]) - 1))
      if (r %in% graph$admix_nodes) {
        if (graph$parent_order[[r]][2] %in% cycles[[mix]][[right_start]] ||
            graph$parent_order[[r]][2] == cycles[[mix]]$collision) {
          graph$parent_order[[r]] <- move(graph$parent_order[[r]], 1, 1)
        }
        if (cycles[[mix]]$collision == cycles[[r]]$collision) {
          for (v in graph$child_order[[cycles[[mix]]$collision]]) {
            if (v %in% cycles[[r]][[graph$parent_order[[r]][2]]]) {
              target <- match(right_last, graph$child_order[[cycles[[mix]]$collision]])
              own_position <- match(v, graph$child_order[[cycles[[mix]]$collision]])
              graph$child_order[[cycles[[mix]]$collision]] <-
                move(graph$child_order[[cycles[[mix]]$collision]], own_position, target - own_position + 1)
            }
          }
        }
      }
    }
  }
  print("NEW ORDER:")
  print(graph$child_order)
  print(graph$parent_order)
}

all_paths_to_root <- function(graph, node) {
  path_list <- list()
  if (node == graph$root) {
    path_list[[1]] <- c(node)
  } else if (node %in% graph$admix_nodes) {
    previous <- c(all_paths_to_root(graph, graph$parent_order[[node]][1]),
                  all_paths_to_root(graph, graph$parent_order[[node]][2]))
    for (j in seq(1, length(previous))) {
      path_list[[j]] <- c(node, previous[[j]])
    }
  } else {
    previous <- all_paths_to_root(graph, graph$parent_of[node])
    for (j in seq(1, length(previous))) {
      path_list[[j]] <- c(node, previous[[j]])
    }
  }
  path_list
}

is_descendant_of <- function(graph, offspring, ancestor) {
  answer <- FALSE
  all_ancestors <- all_paths_to_root(graph, offspring)
  for (i in seq(1, length(all_ancestors))) {
    for (j in seq(1, length(all_ancestors[[i]]))) {
      if (all_ancestors[[i]][j] == ancestor) {
        answer <- TRUE
      }
    }
  }
  # However, offspring does not count as a descendant of itself.
  if (offspring == ancestor) {
    answer <- FALSE
  }
  answer
}

first_collision <- function(v1, v2) {
  # Assuming no-one is one's own ancestor.
  answer <- NA
  for (i in seq(1, length(v1))) {
    for (j in seq(1, length(v2))) {
      if (v1[i] == v2[j] && is.na(answer) == TRUE) {
        answer <- v1[i]
      }
    }
  }
  answer
}

branch_count <- function(graph, v1, v2) {
  # Assuming v1 and v2 don't share endpoints.
  count <- 0
  if (length(v1) > 0) {
    for (i in seq(1, length(v1))) {
      count <- count + length(graph$child_order[[v1[i]]]) - 1
      # Not sure if we should but also count the admix branches.
      if (v1[i] %in% graph$admix_nodes) {
        count <- count + 1
      }
    }
  }
  if (length(v2) > 0) {
    for (j in seq(1, length(v2))) {
      count <- count + length(graph$child_order[[v2[j]]]) - 1
      # Not sure if we should but also count the admix branches.
      if (v2[j] %in% graph$admix_nodes) {
        count <- count + 1
      }
    }
  }
  count
}

# In vector v, move n:th element x steps to the right. Will cause problems if used wrong.
move <- function(v, n, x) {
  for (j in seq(1, abs(x))) {
    temp <- v[n + sign(x)*j]
    v[n + sign(x)*j] <- v[n + sign(x)*(j - 1)]
    v[n + sign(x)*(j - 1)] <- temp
  }
  v
}