#' Fast version of graph plotting.
#' 
#' This is a fast, deterministic and stand-alone function for visualizing the
#' admixture graph. Has the bad habit if sometimes drawing several nodes at the
#' exact same coordinates; for clearer reasults try \code{\link{plot.agraph}}
#' (which, on the other hand, relies on numerical optimising of a compicated cost
#' function and might be unpredictable).
#' 
#' @param x                       The admixture graph.
#' @param ordered_leaves          The leaf-nodes in the left to right order they
#'                                should be drawn.
#' @param show_admixture_labels   A flag determining if the plot should include
#'                                the names of admixture proportions.
#' @param show_inner_node_labels  A flag determining if the plot should include
#'                                the names of inner nodes.
#' @param ...                     Additional plotting options.
#'
#' @return A plot.
#'
#' @seealso \code{\link{plot.agraph}}
#' 
#' @examples
#' # Here is an example of a graph that doesn't behave right under fast_plot(),
#' # taken from the collection of all the admixture graphs with four leaves and at
#' # most two admixture events:
#' 
#' fast_plot(four_leaves_graphs[[24]](c("A", "B", "C", "D")))
#' 
#' # To be fair, here is a graph that looks all right:
#' 
#' fast_plot(four_leaves_graphs[[25]](c("A", "B", "C", "D")))
#' 
#' @export
fast_plot <- function(x, 
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
  graphics::plot(xpos, ypos, type = "n", axes = FALSE, frame.plot = FALSE,
                 xlab = "", ylab = "", ylim = c(-1, max(ypos) + 0.5), ...)
  for (node in graph$nodes) {
    parents <- graph$nodes[graph$parents[node, ]]
    if (length(parents) == 1) {
      graphics::lines(c(xpos[node],xpos[parents]), c(ypos[node], ypos[parents]))
    } else if (length(parents) == 2) {
      break_y <- ypos[node]
      break_x_left <- xpos[node] - 0.3
      break_x_right <- xpos[node] + 0.3
      if (xpos[parents[1]] < xpos[parents[2]]) {
        graphics::lines(c(xpos[parents[1]], break_x_left), c(ypos[parents[1]], break_y))
        graphics::lines(c(xpos[parents[2]], break_x_right), c(ypos[parents[2]], break_y))  
      } else {
        graphics::lines(c(xpos[parents[2]], break_x_left), c(ypos[parents[2]], break_y))
        graphics::lines(c(xpos[parents[1]], break_x_right), c(ypos[parents[1]], break_y))
      }
      graphics::segments(break_x_left, break_y, xpos[node], ypos[node], col = "red")
      graphics::segments(break_x_right, break_y, xpos[node], ypos[node], col = "red")
      if (show_admixture_labels) {
        if (xpos[parents[1]] < xpos[parents[2]]) {
          graphics::text(break_x_left, break_y, graph$probs[parents[[1]], node],
                         cex = 0.5, pos = 1, col = "red", offset = 0.1)
          graphics::text(break_x_right, break_y, graph$probs[parents[[2]], node],
                         cex = 0.5, pos = 1, col = "red", offset = 0.1)
        } else {
          graphics::text(break_x_left, break_y, graph$probs[parents[[2]], node],
                         cex = 0.5, pos = 1, col = "red", offset = 0.1)
          graphics::text(break_x_right, break_y, graph$probs[parents[[1]], node],
                         cex = 0.5, pos = 1, col = "red", offset = 0.1)          
        }
      }
    }
  }
  is_inner <- Vectorize(function(n) sum(graph$children[n, ]) > 0)
  inner_nodes <- which(is_inner(graph$nodes))
  leaves <- which(!is_inner(graph$nodes))
  if (show_inner_node_labels) {
    graphics::text(xpos[inner_nodes], ypos[inner_nodes], 
                   labels = graph$nodes[inner_nodes], cex = 0.6, col = "blue", pos = 3)
  }
  graphics::text(xpos[leaves], ypos[leaves], labels = graph$nodes[leaves], 
                 cex = 0.7, col = "black", pos = 1)
  invisible()
}

#' Plot an admixture graph.
#' 
#' This is a basic drawing routine for visualising the graph. Uses Nelder-Mead
#' algorithm and complicated heuristic approach to find aestethic node coordinates,
#' but due to bad luck or an oversight in the heuristics, especially with larger
#' graphs, might sometimes produce a weird looking result. Usually plotting again
#' helps and if not, use the optional parameters to help the algorithm or try the
#' faster and deterministic function \code{\link{fast_plot}} (which unfortunately is
#' not very good at handling multiple admixture events).
#'
#' @param x                       The admixture graph.
#' @param show_leaf_labels        A flag determining if leaf names are shown.
#' @param draw_leaves             A flag determining if leaf nodes are drawn as
#'                                little circles.
#' @param color                   Color of all drawn nodes unless overriden by
#'                                \code{inner_node_color}.
#' @param show_inner_node_labels  A flag determining if the plot should include
#'                                the names of inner nodes.
#' @param draw_inner_nodes        A flag determining if inner nodes are drawn as
#'                                little circles.
#' @param inner_node_color        Color of inner node circles, if drawn.
#' @param show_admixture_labels   A flag determining if the plot should include
#'                                the names of admixture proportions.
#' @param parent_order            An optional list of instuctions on which order
#'                                from left to right to draw the parents of nodes.
#'                                The list should contain character vectors of parents
#'                                with the name of the child, \emph{e.g.}
#'                                \code{child = c("left_parent", "right_parent")}.
#'                                Using automated heuristics for nodes not specified.
#' @param child_order             An optional list of instuctions on which order
#'                                from left to right to draw the children of nodes.
#'                                The list should contain character vectors of children
#'                                with the name of the parent, \emph{e.g.}
#'                                \code{parent = c("left_child", "right_child")}.
#'                                Using automated heuristics for nodes not specified.
#' @param leaf_order              An optional vector describing in which order should
#'                                the leaves be drawn. Using automated heuristic
#'                                depending on \code{parent_order} and \code{child_order}
#'                                if not specified. Accepts both a character vector of
#'                                the leaves or a numeric vector interpreted as a
#'                                permutation of the default order.
#' @param fix                     If nothing else helps, the list \code{fix} can be used to
#'                                correct the inner node coordinates given by the heuristics.
#'                                Should contain numeric vectors of length 2 with the name
#'                                of an inner node, \emph{e.g.} \code{inner_node = c(0, 10)},
#'                                moving \code{inner_node} to the right 10 units where 100 is
#'                                the plot width. Non-specified inner nodes are left in peace.
#' @param platform                By default admixture nodes are drawn with a horizontal
#'                                platform for proportion labels, the width of which is
#'                                half the distance between any two leaves. The number
#'                                \code{platform} tells how many default platform widths
#'                                should the platforms be wide, \emph{i. e.} zero means no
#'                                platform.
#' @param title                   Optional title for the plot.
#' @param ...                     Additional plotting options.
#'
#' @return A plot.
#'
#' @seealso \code{\link{agraph}}
#' @seealso \code{\link{fast_plot}}
#'
#' @examples
#' \donttest{
#' leaves <- c("salmon", "sea horse", "mermaid", "horse", "human", "lobster")
#' inner_nodes <- c("R", "s", "t", "u", "v", "w", "x", "y", "z")
#' edges <- parent_edges(c(edge("salmon", "t"),
#'                         edge("sea horse", "y"),
#'                         edge("mermaid", "z"),
#'                         edge("horse", "w"),
#'                         edge("human", "x"),
#'                         edge("lobster", "R"),
#'                         edge("s", "R"),
#'                         edge("t", "s"),
#'                         edge("u", "t"),
#'                         edge("v", "s"),
#'                         edge("w", "v"),
#'                         edge("x", "v"),
#'                         admixture_edge("y", "u", "w"),
#'                         admixture_edge("z", "u", "x")))
#' admixtures <- admixture_proportions(c(admix_props("y", "u", "w", "a"),
#'                                       admix_props("z", "u", "x", "b")))
#' graph <- agraph(leaves, inner_nodes, edges, admixtures)
#' plot(graph, show_inner_node_labels = TRUE)
#' 
#' # Suppose now that we prefer to have the outgroup "lobster" on the right side.
#' # This is achieved by telling the algorithm that the children of "R" should be in
#' # the order ("s", "lobster"), from left to right:
#' 
#' plot(graph, child_order = list(R = c("s", "lobster")))
#' 
#' # Suppose further that we prefer to have "mermaid" and "human" next to each other,
#' # as well as "sea horse" and "horse". This is easily achieved by rearranging the leaf
#' # order proposed by the algorithm. We can also fine-tune by moving "y" a little bit
#' # to the right, make the admixture platforms a bit shorter, color the nodes, show the
#' # admixture proportions and give the plot a title:
#'
#' plot(graph, child_order = list(R = c("s", "lobster")), leaf_order = c(1, 2, 4, 3, 5, 6),
#'      fix = list(y = c(5, 0)), platform = 0.8, color = "deepskyblue",
#'      inner_node_color = "skyblue", show_admixture_labels = TRUE,
#'      title = "Evolution of fish/mammal hybrids")
#' }
#' 
#' @export
plot.agraph <- function(x,
                        show_leaf_labels = TRUE,
                        draw_leaves = TRUE,
                        color = "yellowgreen",
                        show_inner_node_labels = FALSE,
                        draw_inner_nodes = draw_leaves,
                        inner_node_color = color,
                        show_admixture_labels = FALSE,
                        parent_order = list(),
                        child_order = list(),
                        leaf_order = NULL,
                        fix = list(),
                        platform = 1,
                        title = NULL,
                        ...) {
  # Combine the user instructions and automated heuristics about the graph orderings.
  graph <- x
  arranged <- arrange_graph(graph)
  parent_order <- utils::modifyList(arranged$parent_order, parent_order)
  child_order <- utils::modifyList(arranged$child_order, child_order)
  if (is.null(leaf_order) == TRUE) {
    leaf_order <- leaf_order(graph, parent_order, child_order)
  } else if (typeof(leaf_order) != "character") {
    machine_order <- leaf_order(graph, parent_order, child_order)
    human_order <- leaf_order
    leaf_order <- rep("", length(leaf_order))
    for (i in seq(1, length(human_order))) {
      leaf_order[i] <- machine_order[human_order[i]]
    }
  }
  # Assign initial coordinates for all nodes.
  leaves <- list()
  for (i in seq(1, length(leaf_order))) {
    leaves[[i]] <- c(100*(i - 1)/(length(leaf_order) - 1), 0)
  }
  names(leaves) <- leaf_order
  parents <- graph$parents
  for (i in seq(1, length(graph$inner_nodes))) {
    candidate <- graph$inner_nodes[i]
    abandon <- FALSE
    for (j in seq(1, NCOL(parents))) {
      if (parents[candidate, j] == TRUE) {
        abandon <- TRUE
      }
    }
    if (abandon == FALSE) {
      root <- list(c(50, 100))
      names(root) <- c(candidate)
      delete <- i
    }
  }
  root_removed <- graph$inner_nodes[-delete]
  inner <- list()
  if (length(root_removed) > 0) {
    for (i in seq(1, length(root_removed))) {
      inner[[i]] <- c(0, 0)
    }
  }
  names(inner) <- root_removed
  # Assign the y-coordinates according to my arbitrary principles.
  refined_graph <- refined_graph(graph)
  heights <- rep(0, length(inner))
  names(heights) <- names(inner)
  global_longest <- 0
  for (inner_node in names(inner)) {
    paths <- all_paths_to_leaves(refined_graph, inner_node)
    longest <- 0
    for (path in paths) {
      if (length(path) > longest) {
        longest <- length(path)
      }
    }
    heights[inner_node] <- longest - 1
    if (longest > global_longest) {
      global_longest <- longest
    }
  }
  for (inner_node in names(inner)) {
    heights[inner_node] <- global_longest - heights[inner_node]
    inner[[inner_node]][2] <- 100*(1 - heights[inner_node]/global_longest)
  }
  # Perform Nelder-Mead to optimize the x-coordinates of the non-root inner nodes.
  if (length(inner) > 0) {
    x0 <- rep(50, length(inner))
    min <- rep(0, length(inner))
    max <- rep(100, length(inner))
    cfunc <- drawing_cost(graph, leaves, root, inner, child_order, parent_order, platform)
    opti <- suppressWarnings(neldermead::fminbnd(cfunc, x0 = x0, xmin = min, xmax = max))
    x <- neldermead::neldermead.get(opti, "xopt")
    for (i in seq(1, length(inner))) {
      inner[[i]][1] <- x[i]
    }
  }
  # Plot everything asked for.
  xpd <- graphics::par()$xpd
  graphics::par(xpd = NA)
  level <- platform*25/(length(leaves) - 1)
  for (inner_node in names(inner)) {
    inner[[inner_node]][2] <- 100*(1 - heights[inner_node]/global_longest)
  }
  for (fixed in names(fix)) {
    inner[[fixed]] <- inner[[fixed]] + fix[[fixed]]
  }
  nodes <- graph$nodes
  coordinates <- c(leaves, root, inner)
  graphics::plot(c(-level, 100 + level), c(0, 100), type = "n", axes = FALSE,
                 frame.plot = FALSE, xlab = "", ylab = "", main = title, ...)
  for (i in nodes) {
    for (j in nodes) {
      if (parents[i, j] == TRUE) {
        i_thing <- 0
        if (length(parent_order[[i]]) == 2) {
          if (parent_order[[i]][1] == j) {
            i_thing <- -level
            graphics::segments(coordinates[[i]][1] + i_thing, coordinates[[i]][2], coordinates[[i]][1] - i_thing,
                               coordinates[[i]][2], col = "black", lwd = 2)
          }
          if (parent_order[[i]][2] == j) {
            i_thing <- level
          }
          if (show_admixture_labels == TRUE) {
            label <- graph$probs[i, j]
            if (substr(label, 1, 1) == "(") {
              label <- substr(label, 2, nchar(label) - 1)
            }
            graphics::text(coordinates[[i]][1] + 0.75*i_thing, coordinates[[i]][2], label,
                           adj = c(0.5, 1.6), cex = 0.8)
          }
        }
        graphics::segments(coordinates[[i]][1] + i_thing, coordinates[[i]][2], coordinates[[j]][1],
                           coordinates[[j]][2], col = "black", lwd = 2)
      }
    }
  }
  for (i in seq(1, length(leaves))) {
    leaf <- leaves[[i]]
    if (draw_leaves == TRUE) {
      graphics::points(leaf[1], leaf[2], lwd = 2, pch = 21, col = "black", bg = color, cex = 2)
    }
    if (show_leaf_labels == TRUE) {
      graphics::text(leaf[1], leaf[2], names(leaves)[i], adj = c(0.5, 2.6), cex = 0.8)
    }
  }
  if (length(inner) > 0) {
    for (i in seq(1, length(inner))) {
      vertex <- inner[[i]]
      if (draw_inner_nodes == TRUE) {
        graphics::points(vertex[1], vertex[2], lwd = 2, pch = 21, col = "black", bg = inner_node_color, cex = 2)
      }
      if (show_inner_node_labels == TRUE) {
        graphics::text(vertex[1], vertex[2], names(inner)[i], adj = c(0.5, -1.6), cex = 0.8)
      }
    }
  }
  juuri <- root[[1]]
  if (draw_inner_nodes == TRUE) {
    graphics::points(juuri[1], juuri[2], lwd = 2, pch = 21, col = "black", bg = inner_node_color, cex = 2)
  }
  if (show_inner_node_labels == TRUE) {
    graphics::text(juuri[1], juuri[2], names(root)[1], adj = c(0.5, -1.6), cex = 0.8)
  }
  graphics::par(xpd = xpd)
}

drawing_cost <- function(graph, leaves, root, inner, child_order, parent_order, platform) {
  function(x) {
    # Put the variable x in place.
    for (i in seq(1, length(inner))) {
      inner[[i]][1] <- x[i]
    }
    all <- c(leaves, root, inner)
    # Calculate the cost of the input and default coordinates.
    cost <- 0
    new <- function(u, v) {
      w <- u + v
      U <- sqrt(u[1]^2 + u[2]^2)
      V <- sqrt(v[1]^2 + v[2]^2)
      W <- sqrt(w[1]^2 + w[2]^2)
      kos <- u%*%v/(U*V + 1e-5)
      if (v[1]*abs(u[2]) > u[1]*abs(v[2])) {
        angle <- abs(kos)^4 # ARBITRARY CONSTANTS HERE!
        if (angle > 0.95) {
          angle <- angle + 5
        }
      } else {
        angle <- 10 - kos
      }
      verticality <- abs(w%*%c(0, 1)/(W + 1e-5))  
      return((4 + 4*angle - 2*verticality)*(U + V)) # ARBITRARY CONSTANTS HERE!
    }
    new2 <- function(u) {
      U <- sqrt(u[1]^2 + u[2]^2)
      verticality <- abs(u%*%c(0, 1)/(U + 1e-5))  
      return((4 - 2*verticality)*U) # ARBITRARY CONSTANTS HERE!
    }
    level <- platform*25/(length(leaves) - 1)
    for (i in seq(1, length(inner))) {
      p <- parent_order[[names(inner)[i]]]
      c <- child_order[[names(inner)[i]]]
      if (length(p) == 2) {
        u <- c(all[[p[1]]][1] - inner[[i]][1] - level, all[[p[1]]][2] - inner[[i]][2])
        v <- c(all[[p[2]]][1] - inner[[i]][1] + level, all[[p[2]]][2] - inner[[i]][2])
        cost <- cost + new(u, v)
        if (length(parent_order[[c[1]]]) == 2) {
          if (parent_order[[c[1]]][1] == names(inner)[[i]]) {
            thing <- -level
          }
          if (parent_order[[c[1]]][2] == names(inner)[[i]]) {
            thing <- level
          }
        } else {
          thing <- 0
        }
        u <- c(all[[c[1]]][1] - inner[[i]][1] + thing, all[[c[1]]][2] - inner[[i]][2])
        cost <- cost + new2(u)
      }
      if (length(c) == 2) {
        if (length(parent_order[[c[1]]]) == 2) {
          if (parent_order[[c[1]]][1] == names(inner)[[i]]) {
            thing <- -level
          }
          if (parent_order[[c[1]]][2] == names(inner)[[i]]) {
            thing <- level
          }
        } else {
          thing <- 0
        }
        u <- c(all[[c[1]]][1] - inner[[i]][1] + thing, all[[c[1]]][2] - inner[[i]][2])
        if (length(parent_order[[c[2]]]) == 2) {
          if (parent_order[[c[2]]][1] == names(inner)[[i]]) {
            thing <- -level
          }
          if (parent_order[[c[2]]][2] == names(inner)[[i]]) {
            thing <- level
          }
        } else {
          thing <- 0
        }
        v <- c(all[[c[2]]][1] - inner[[i]][1] + thing, all[[c[2]]][2] - inner[[i]][2])
        cost <- cost + new(u, v)
        u <- c(all[[p[1]]][1] - inner[[i]][1], all[[p[1]]][2] - inner[[i]][2])
        cost <- cost + new2(u)
      }
    }
    # Root also.
    c <- child_order[[names(root)[1]]]
    if (length(parent_order[[c[1]]]) == 2) {
      if (parent_order[[c[1]]][1] == names(inner)[[i]]) {
        thing <- -level
      }
      if (parent_order[[c[1]]][2] == names(inner)[[i]]) {
        thing <- level
      }
    } else {
      thing <- 0
    }
    u <- c(all[[c[1]]][1] - root[[1]][1] + thing, all[[c[1]]][2] - root[[1]][2])
    if (length(parent_order[[c[2]]]) == 2) {
      if (parent_order[[c[2]]][1] == names(inner)[[i]]) {
        thing <- -level
      }
      if (parent_order[[c[2]]][2] == names(inner)[[i]]) {
        thing <- level
      }
    } else {
      thing <- 0
    }
    v <- c(all[[c[2]]][1] - root[[1]][1] + thing, all[[c[2]]][2] - root[[1]][2])
    cost <- cost + new(u, v)
    return(cost)
  } 
}

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

arrange_graph <- function(graph) {
  graph <- refined_graph(graph)
  # Determining the crossing numbers of a graph is a difficult problem, and requiring
  # edges to be straight lines makes it even worse. The following algorithm is by no
  # means meant to be optimal. We only perform some clean-up on the cycles of the graph,
  # using a heuristic of dealing with the worst cycles first.
  
  # First determine cycles corresponding to admix nodes.  We follow both branches until 
  # a common acestor is found. Unfortunately this is not unique if we hit more admix nodes 
  # on the way. Then our first tie-breaker is the age of the ancestor, if one ancestor is an
  # ancestor of the other, choose the younger one, and the second tie-breaker is the branch
  # count of the cycle, the number of branches deviating from it, up or down, not counting
  # the admix node itself or the ancestor (doesn' matter when breaking ties of course but
  # we will use the branch count later also).
  cycles <- list()
  if (length(graph$admix_nodes) > 0) {
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
            temp[["left"]] <- c(mix, v1, collision)
            temp[["right"]] <- c(mix, v2, collision)
            temp[["collision"]] <- collision
            temp[["branch_count"]] <- branch_count
            temp[["admix"]] <- mix
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
        }
        if (cycle_candidates[[i]]$branch_count < least_branch_count && could_this_be_it == TRUE) {
          cycles[[mix]] <- cycle_candidates[[i]]
          least_branch_count <- cycle_candidates[[i]]$branch_count
        }
      }
    }
    # Now that the cycles are determined, we need to make sure they really are cycles and not "eights",
    # that is, change the parent order of the admix node to match the child order of the collision node.
    # As several cycles might share the same collision node but not the admix node, the initial orientation
    # has to be determined by the collision node.
    # Always keep to left side of the cycle at [[1]] and the right side at [[2]].
    for (i in seq(1, length(cycles))) {
      cycle <- cycles[[i]]
      order <- graph$child_order[[cycle$collision]]
      if (cycle$left[length(cycle$left) - 1] != order[1]) {
        cycles[[i]]$left <- cycle$right
        cycles[[i]]$right <- cycle$left
      }
      graph$parent_order[[cycle$admix]] <- c(cycles[[i]]$left[2], cycles[[i]]$right[2])
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
    if (length(cleaning_order) > 1) {
      for (i in seq(2, length(cleaning_order))) {
        for (j in seq(1, i - 1)) {
          if (is_descendant_of(graph, cleaning_order[i - j], cleaning_order[i - j + 1]) == FALSE) {
            if (cycles[[cleaning_order[i - j + 1]]]$branch_count > cycles[[cleaning_order[i - j]]]$branch_count) {
              tempo <- cleaning_order[i - j]
              cleaning_order[i - j] <- cleaning_order[i - j + 1]
              cleaning_order[i - j + 1] <- tempo
            }
          }
        }
      }
    }
    # Then accorning to the cleaning order, turn all branches away from the cycle.
    # If we encounter other admix nodes or parts of other cycles that do not share the
    # collision node, make sure the cycles have different orientations.
    # If we encounter parts of other cycles that do share the collision node, it is
    # probably better to draw the cycles intersecting (one crossing) than it is to
    # draw one inside of the other (one or more crossings plus coordinate troubles), so
    # some turnings need to be omitted.
    for (mix in cleaning_order) { # Hey why aren't all my loops like this?
      # Left side:
      for (l in seq(1, length(cycles[[mix]]$left) - 1)) {
        vertex <- cycles[[mix]]$left[l]
        if (length(graph$child_order[[vertex]]) == 2) {
          # So let's try this stunt:
          # Only here on the left side of the cycle, if the vertex we see belongs to
          # cycle with the same collision node, turn it inside instead of outside and
          # otherwise do everything as before.
          turn_in <- FALSE
          for (other_cycle in cycles) {
            other_left <- other_cycle$left[-length(other_cycle$left)]
            other_right <- other_cycle$right[-length(other_cycle$right)]
            other <- c(other_left, other_right)
            if (vertex %in% other && other_cycle$admix != mix &&
                other_cycle$collision == cycles[[mix]]$collision) {
              if (graph$child_order[[vertex]][1] %in% other_cycle == TRUE &&
                  graph$child_order[[vertex]][1] %in% cycles[[mix]]$left == FALSE) {
                turn_in <- TRUE
              } 
              if (graph$child_order[[vertex]][2] %in% other_cycle == TRUE &&
                  graph$child_order[[vertex]][2] %in% cycles[[mix]]$left == FALSE) {
                turn_in <- TRUE
              }
            }
          }
          if (turn_in == FALSE) {
            if (graph$child_order[[vertex]][2] != cycles[[mix]]$left[l - 1]) {
              graph$child_order[[vertex]][1] <- graph$child_order[[vertex]][2]
              graph$child_order[[vertex]][2] <- cycles[[mix]]$left[l - 1]
            }
          } else {
            if (graph$child_order[[vertex]][1] != cycles[[mix]]$left[l - 1]) {
              graph$child_order[[vertex]][2] <- graph$child_order[[vertex]][1]
              graph$child_order[[vertex]][1] <- cycles[[mix]]$left[l - 1]
              # Any cycle we need to turn around because of breaking the rules like this?
              # Possibly yes, but avoid cleaning them too much.
              for (first_cycle in cycles) {
                if (first_cycle$collision == vertex) {
                  graph$child_order[[vertex]] <- c(graph$parent_order[[vertex]][2],
                                                   graph$parent_order[[vertex]][1])
                  cycles[[first_cycle$admix]]$left <- first_cycle$right
                  cycles[[first_cycle$admix]]$right <- first_cycle$left
                  rest <- character(0)
                  for (other_cycle in cycles) {
                    if (other_cycle$admix != first_cycle$admix) {
                      rest <- c(rest, other_cycle$left, other_cycle$right)
                    }
                  }
                  for (L in seq(1, length(first_cycle$left) - 1)) {
                    mess <- first_cycle$left[L]
                    if (mess %in% rest == FALSE && length(graph$child_order[[mess]] == 2)) {
                      if (graph$child_order[[mess]][2] != first_cycle$left[L - 1]) {
                        graph$child_order[[mess]][1] <- graph$child_order[[mess]][2]
                        graph$child_order[[mess]][2] <- first_cycle$left[L - 1]
                      }
                    }
                  }
                  for (R in seq(1, length(first_cycle$right) - 1)) {
                    mess <- first_cycle$right[R]
                    if (mess %in% rest == FALSE && length(graph$child_order[[mess]] == 2)) {
                      if (graph$child_order[[mess]][1] != first_cycle$right[R - 1]) {
                        graph$child_order[[mess]][2] <- graph$child_order[[mess]][1]
                        graph$child_order[[mess]][1] <- first_cycle$right[R - 1]
                      }
                    }
                  }
                }
              }  
            }
          }
        }
        if (length(graph$parent_order[[vertex]]) == 2 && l > 1) {
          # The collision node can't be shared here between the two cycles.
          if (graph$parent_order[[vertex]][2] != cycles[[mix]]$left[l + 1]) {
            graph$parent_order[[vertex]][1] <- graph$parent_order[[vertex]][2]
            graph$parent_order[[vertex]][2] <- cycles[[mix]]$left[l + 1]
            other_cycle <- cycles[[vertex]]
            cycles[[vertex]]$left <- other_cycle$right
            cycles[[vertex]]$right <- other_cycle$left
            graph$child_order[[other_cycle$collision]] <- c(other_cycle$right[length(other_cycle$right) - 1],
                                                            other_cycle$left[length(other_cycle$left) - 1])
          }
        }
        # If two cycles share some part, we should at least suggest that the common part belongs to one left and
        # one right side.
        for (i in seq(1, length(cycles))) {
          other_cycle <- cycles[[i]]
          if (cycles[[mix]]$collision != other_cycle$collision && l != 1) {
            # Also ensures the cycles are different.
            other_left <- other_cycle$left[-length(other_cycle$left)][-1]
            if (vertex %in% other_left) {
              cycles[[i]]$left <- other_cycle$right
              cycles[[i]]$right <- other_cycle$left
              graph$child_order[[other_cycle$collision]] <- c(other_cycle$right[length(other_cycle$right) - 1],
                                                              other_cycle$left[length(other_cycle$left) - 1])
              graph$parent_order[[other_cycle$admix]] <- c(other_cycle$right[2], other_cycle$left[2])
              # To avoid endless looping, we mustn't perform a full clearing on the cycle we just switched around.
              # Instead we clear only such child edges that are not part of another cycle.
              for (L in seq(1, length(cycles[[i]]$left) - 1)) {
                old <- cycles[[i]]$left[[L]]
                if (length(graph$child_order[[old]]) == 2) {
                  if (graph$child_order[[old]][2] != cycles[[i]]$left[L - 1]) {
                    young <- graph$child_order[[old]][2]
                    let_it_be <- FALSE
                    for (yet_another_cycle in cycles) {
                      if (old %in% yet_another_cycle$left && young %in% yet_another_cycle$left) {
                        let_it_be <- TRUE
                      }
                      if (old %in% yet_another_cycle$right && young %in% yet_another_cycle$right) {
                        let_it_be <- TRUE
                      }
                    }
                    if (let_it_be == FALSE) {
                      graph$child_order[[old]] <- c(young, cycles[[i]]$left[L - 1])
                    }
                  }
                }
              }
              for (R in seq(1, length(cycles[[i]]$right) - 1)) {
                old <- cycles[[i]]$right[[R]]
                if (length(graph$child_order[[old]]) == 2) {
                  if (graph$child_order[[old]][1] != cycles[[i]]$right[R - 1]) {
                    young <- graph$child_order[[old]][1]
                    let_it_be <- FALSE
                    for (yet_another_cycle in cycles) {
                      if (old %in% yet_another_cycle$left && young %in% yet_another_cycle$left) {
                        let_it_be <- TRUE
                      }
                      if (old %in% yet_another_cycle$right && young %in% yet_another_cycle$right) {
                        let_it_be <- TRUE
                      }
                    }
                    if (let_it_be == FALSE) {
                      graph$child_order[[old]] <- c(cycles[[i]]$right[R - 1], young)
                    }
                  }
                }
              }
            }
          }
        }
      }
      # Right side:
      for (r in seq(1, length(cycles[[mix]]$right) - 1)) {
        vertex <- cycles[[mix]]$right[r]
        if (length(graph$child_order[[vertex]]) == 2) {
          if (graph$child_order[[vertex]][1] != cycles[[mix]]$right[r - 1]) {
            graph$child_order[[vertex]][2] <- graph$child_order[[vertex]][1]
            graph$child_order[[vertex]][1] <- cycles[[mix]]$right[r - 1]
          }
        }
        if (length(graph$parent_order[[vertex]]) == 2 && r > 1) {
          # The collision node can't be shared here between the two cycles.
          if (graph$parent_order[[vertex]][1] != cycles[[mix]]$right[r + 1]) {
            graph$parent_order[[vertex]][2] <- graph$parent_order[[vertex]][1]
            graph$parent_order[[vertex]][1] <- cycles[[mix]]$right[r + 1]
            other_cycle <- cycles[[vertex]]
            cycles[[vertex]]$left <- other_cycle$right
            cycles[[vertex]]$right <- other_cycle$left
            graph$child_order[[other_cycle$collision]] <- c(other_cycle$right[length(other_cycle$right) - 1],
                                                            other_cycle$left[length(other_cycle$left) - 1])
          }
        }
        # If two cycles share some part, we should at least suggest that the common part belongs to one left and
        # one right side.
        for (i in seq(1, length(cycles))) {
          other_cycle <- cycles[[i]]
          if (cycles[[mix]]$collision != other_cycle$collision && r != 1) {
            # Also ensures the cycles are different.
            other_right <- other_cycle$right[-length(other_cycle$right)][-1]
            if (vertex %in% other_right) {
              cycles[[i]]$left <- other_cycle$right
              cycles[[i]]$right <- other_cycle$left
              graph$child_order[[other_cycle$collision]] <- c(other_cycle$right[length(other_cycle$right) - 1],
                                                              other_cycle$left[length(other_cycle$left) - 1])
              graph$parent_order[[other_cycle$admix]] <- c(other_cycle$right[2], other_cycle$left[2])
              # To avoid endless looping, we mustn't perform a full clearing on the cycle we just switched around.
              # Instead we clear only such child edges that are not part of another cycle.
              for (L in seq(1, length(cycles[[i]]$left) - 1)) {
                old <- cycles[[i]]$left[[L]]
                if (length(graph$child_order[[old]]) == 2) {
                  if (graph$child_order[[old]][2] != cycles[[i]]$left[L - 1]) {
                    young <- graph$child_order[[old]][2]
                    let_it_be <- FALSE
                    for (yet_another_cycle in cycles) {
                      if (old %in% yet_another_cycle$left && young %in% yet_another_cycle$left) {
                        let_it_be <- TRUE
                      }
                      if (old %in% yet_another_cycle$right && young %in% yet_another_cycle$right) {
                        let_it_be <- TRUE
                      }
                    }
                    if (let_it_be == FALSE) {
                      graph$child_order[[old]] <- c(young, cycles[[i]]$left[L - 1])
                    }
                  }
                }
              }
              for (R in seq(1, length(cycles[[i]]$right) - 1)) {
                old <- cycles[[i]]$right[[R]]
                if (length(graph$child_order[[old]]) == 2) {
                  if (graph$child_order[[old]][1] != cycles[[i]]$right[R - 1]) {
                    young <- graph$child_order[[old]][1]
                    let_it_be <- FALSE
                    for (yet_another_cycle in cycles) {
                      if (old %in% yet_another_cycle$left && young %in% yet_another_cycle$left) {
                        let_it_be <- TRUE
                      }
                      if (old %in% yet_another_cycle$right && young %in% yet_another_cycle$right) {
                        let_it_be <- TRUE
                      }
                    }
                    if (let_it_be == FALSE) {
                      graph$child_order[[old]] <- c(cycles[[i]]$right[R - 1], young)
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  # One last thing: we want to record a single parent in parent_order just likewe record a single child
  # to child_order.
  for (node in graph$nodes) {
    if (is.na(graph$parent_of[node]) == FALSE) {
      graph$parent_order[[node]] <- graph$parent_of[node]
    }
  }
  return(list(parent_order = graph$parent_order, child_order = graph$child_order, cycles = cycles))
}

leaf_order <- function(graph, parent_order, child_order) {
  # We actually still need the refined graph object, so build it and fix the orderings it has by default.
  graph <- refined_graph(graph)
  graph$parent_order <- parent_order
  graph$child_order <- child_order
  # Give some arrangement for leaves. Doesn't really follow from parent and child orders, but almost.
  # Assign the y-coordinates according to my arbitrary principles.
  heights <- rep(0, length(graph$nodes))
  names(heights) <- graph$nodes
  for (node in graph$nodes) {
    paths <- all_paths_to_root(graph, node)
    longest <- 0
    for (path in paths) {
      if (length(path) > longest) {
        longest <- length(path)
      }
    }
    heights[node] <- longest - 1
  }
  # Calculate some tentative values for the x-coordinates, just to get a goog guess for leaf ordering.
  x_order <- graph$nodes[order(heights)]
  silly_x <- rep(0, length(graph$nodes))
  names(silly_x) <- graph$nodes
  for (i in seq(2, length(x_order))) {
    vertex <- x_order[i]
    if (length(graph$parent_order[[vertex]]) < 2) {
      parent <- graph$parent_order[[vertex]][1]
      # Only child of a single parent:
      if (length(graph$child_order[[parent]]) == 1) {
        silly_x[vertex] <- silly_x[parent]
      }
      if (length(graph$child_order[[parent]]) == 2) {
        # First born of a single parent:
        if (graph$child_order[[parent]][1] == vertex) {
          silly_x[vertex] <- silly_x[parent] - 2^(-heights[vertex] + 1)
        }
        # Second child of a single parent:
        if (graph$child_order[[parent]][2] == vertex) {
          silly_x[vertex] <- silly_x[parent] + 2^(-heights[vertex] + 1)
        }
      }
    } else {
      # Child of two parents:
      left_parent <- graph$parent_order[[vertex]][1]
      right_parent <- graph$parent_order[[vertex]][2]
      if (silly_x[left_parent] < silly_x[right_parent]) {
        # Consistent parents:
        silly_x[vertex] <- 0.5*(silly_x[left_parent] + silly_x[right_parent])
        if (length(graph$child_order[[left_parent]]) == 2 &&
            length(graph$child_order[[right_parent]]) == 2) {
          if (graph$child_order[[left_parent]][1] == vertex) {
            silly_x[vertex] <- silly_x[vertex] - 2^(-heights[vertex] + 2)
          }
          if (graph$child_order[[right_parent]][1] == vertex) {
            silly_x[vertex] <- silly_x[vertex] - 2^(-heights[vertex] + 2)
          }
          if (graph$child_order[[left_parent]][2] == vertex) {
            silly_x[vertex] <- silly_x[vertex] + 2^(-heights[vertex] + 2)
          }
          if (graph$child_order[[right_parent]][2] == vertex) {
            silly_x[vertex] <- silly_x[vertex] + 2^(-heights[vertex] + 2)
          }
        }
      } else {
        # Inconsistent parents:
        decision <- 0
        if (length(graph$child_order[[left_parent]]) == 2) {
          if (graph$child_order[[left_parent]][1] == vertex) {
            decision <- decision - 1
          }
          if (graph$child_order[[left_parent]][2] == vertex) {
            decision <- decision + 1
          }
        }
        if (length(graph$child_order[[right_parent]]) == 2) {
          if (graph$child_order[[right_parent]][1] == vertex) {
            decision <- decision - 1
          }
          if (graph$child_order[[right_parent]][2] == vertex) {
            decision <- decision + 1
          }
        }
        if (decision < 0) {
          silly_x[vertex] <- silly_x[right_parent] - 2^(-heights[vertex] + 1)
        }
        if (decision > 0) {
          silly_x[vertex] <- silly_x[left_parent] + 2^(-heights[vertex] + 1)
        }
        if (decision == 0) {
          silly_x[vertex] <- 0.5*(silly_x[left_parent] + silly_x[right_parent])
        }
      }
    }
  }
  leaf_order <- graph$leaves[order(silly_x)]
  leaf_order <- leaf_order[!is.na(leaf_order)]
  return(leaf_order)
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

all_paths_to_leaves <- function(graph, node) {
  path_list <- list()
  if (node %in% graph$leaves) {
    path_list[[1]] <- c(node)
  } else if (node %in% graph$admix_nodes) {
    previous <- all_paths_to_leaves(graph, graph$child_order[[node]][1])
    for (j in seq(1, length(previous))) {
      path_list[[j]] <- c(node, previous[[j]])
    }
  } else {
    previous <- c(all_paths_to_leaves(graph, graph$child_order[[node]][1]),
                  all_paths_to_leaves(graph, graph$child_order[[node]][2]))
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