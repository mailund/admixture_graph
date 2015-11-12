library(magrittr)

data(bears)
populations <- c("BLK", "PB",
                 "Bar", "Chi1", "Chi2", "Adm1", "Adm2",
                 "Denali", "Kenai", "Sweden") 

four <- bears %>% fit_all_graphs(c("BLK", "PB", "Kenai", "Adm1"))
five <- bears %>% fit_all_graphs(c("BLK", "PB", "Kenai", "Adm1", "Sweden"))

no_admixture_events <- function(fit) {
  sum(rowSums(fit$graph$parents) == 2)
}

sses <- unlist(Map(sum_of_squared_errors, four))

no_admixture_events <- function(x) UseMethod("no_admixture_events")
no_admixture_events.agraph <- function(x) sum(rowSums(x$parents) == 2)
no_admixture_events.agraph_fit <- function(x) no_admixture_events(x$graph)
no_admixture_events.agraph_fit_list <- function(x) unlist(Map(no_admixture_events, x))


sum_of_squared_errors <- function(x) UseMethod("sum_of_squared_errors")
sum_of_squared_errors.agraph_fit <- function(x) x$best_error
sum_of_squared_errors.agraph_fit_list <- function(x) unlist(Map(sum_of_squared_errors, x))

fit_by_admixtures <- function(x) {
  no_admixtures <- no_admixture_events(x)
  sse <- sum_of_squared_errors(x)
  
  by(sse, no_admixtures, min)
}

summary.agraph_fit_list <- function(object, ...) {
  xx <- fit_by_admixtures(object)
  xx
}

fit_by_admixtures(four)
(xx <- summary(four))
str(xx)




