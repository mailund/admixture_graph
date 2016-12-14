#' Get the number of admixture events in a graph.
#' 
#' Get the number of admixture events in a graph.
#' 
#' @param x  The graph.
#' 
#' @return Number of admixture events in the graph.
#' 
#' @export
no_admixture_events <- 
  function(x) UseMethod("no_admixture_events")

#' Get the number of admixture events in a graph.
#' 
#' Get the number of admixture events in a graph.
#' 
#' @param x  The graph.
#' 
#' @return Number of admixture events in the graph.
#' 
#' @export
no_admixture_events.agraph <- 
  function(x) sum(rowSums(x$parents) == 2)

#' Get the number of admixture events in a fitted graph.
#' 
#' Get the number of admixture events in a fitted graph.
#' 
#' @param x  The fitted graph.
#' 
#' @return Number of admixture events in the graph.
#' 
#' @export
no_admixture_events.agraph_fit <- 
  function(x) no_admixture_events(x$graph)

#' Get the number of admixture events in a list of fitted graph.
#' 
#' Get the number of admixture events in a list of fitted graph.
#' 
#' @param x  The graphs.
#' 
#' @return Number of admixture events in the graphs.
#' 
#' @export
no_admixture_events.agraph_fit_list <- 
  function(x) unlist(Map(no_admixture_events, x))

#' Get the sum of squared errors for a fitted graph.
#' 
#' Get the sum of squared errors for a fitted graph.
#' 
#' @param x  The fitted graph.
#' 
#' @return SSE for the fit.
#' 
#' @export
sum_of_squared_errors <- 
  function(x) UseMethod("sum_of_squared_errors")

#' Get the sum of squared errors for a fitted graph.
#' 
#' Get the sum of squared errors for a fitted graph.
#' 
#' @param x  The fitted graph.
#' 
#' @return SSE for the fit.
#' 
#' @export
sum_of_squared_errors.agraph_fit <- 
  function(x) x$best_error

#' Get the sum of squared errors for a list of fitted graph.
#' 
#' Get the sum of squared errors for a list of fitted graph.
#' 
#' @param x  The fitted graphs.
#' 
#' @return SSE for the fits.
#' 
#' @export
sum_of_squared_errors.agraph_fit_list <- 
  function(x) unlist(Map(sum_of_squared_errors, x))

#' Get the tests in the fit where the predictions fall outside of the error bars.
#' 
#' Get the tests in the fit where the predictions fall outside of the error bars.
#' 
#' @param fit    The fitted graph.
#' @param sigma  The width of the error bars.
#' 
#' @return The poorly fitted tests.
#' 
#' @export
poor_fits <- 
  function(fit, sigma = 6) UseMethod("poor_fits")

#' Get the tests in the fit where the predictions fall outside of the error bars.
#' 
#' Get the tests in the fit where the predictions fall outside of the error bars.
#' 
#' @param fit    The fitted graph.
#' @param sigma  The width of the error bars.
#' 
#' @return The poorly fitted tests.
#' 
#' @export
poor_fits.agraph_fit <- function(fit, sigma=6) {
  x <- stats::fitted(fit) 
  x$stderr = x$D/x$Z.value
  x[with(x, (abs(D - graph_f4) > sigma/2*stderr)),]
}

#' Get the tests in the fit where the predictions fall outside of the error bars.
#' 
#' Get the tests in the fit where the predictions fall outside of the error bars.
#' 
#' @param fit    The fitted graphs.
#' @param sigma  The width of the error bars.
#' 
#' @return The poorly fitted tests.
#' 
#' @export
poor_fits.agraph_fit_list <- function(fit, sigma=6) {
  lapply(fit, poor_fits, sigma = sigma)
}

#' Get the number of tests in the fit where the predictions fall outside of the error bars.
#' 
#' Get the number of tests in the fit where the predictions fall outside of the error bars.
#' 
#' @param fit    The fitted graph.
#' @param sigma  The width of the error bars.
#' 
#' @return The poorly fitted tests.
#' 
#' @export
no_poor_fits <- 
  function(fit, sigma=6) UseMethod("no_poor_fits")

#' Get the number of tests in the fit where the predictions fall outside of the error bars.
#' 
#' Get the number of tests in the fit where the predictions fall outside of the error bars.
#' 
#' @param fit    The fitted graph.
#' @param sigma  The width of the error bars.
#' 
#' @return The poorly fitted tests.
#' 
#' @export
no_poor_fits.agraph_fit <- function(fit, sigma=6) {
  nrow(no_poor_fits(fit, sigma))
}

#' Get the number of tests in the fit where the predictions fall outside of the error bars.
#' 
#' Get the number of tests in the fit where the predictions fall outside of the error bars.
#' 
#' @param fit    The fitted graph.
#' @param sigma  The width of the error bars.
#' 
#' @return The poorly fitted tests.
#' 
#' @export
no_poor_fits.agraph_fit_list <- function(fit, sigma=6) {
  vapply(poor_fits(fit, sigma = sigma), nrow, 1)
}