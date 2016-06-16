#' admixturegraph: Visualising and analysing admixture graphs.
#' 
#' The package provides functionality to analyse and test admixture graphs 
#' against the \eqn{f} statistics described in the paper 
#' \href{http://tinyurl.com/o5a4kr4}{Ancient Admixture in Human History},
#' Patternson \emph{et al.}, Genetics, Vol. 192, 1065--1093, 2012.
#' 
#' The \eqn{f} statistics -- \eqn{f_2}, \eqn{f_3}, and \eqn{f_4} -- extract 
#' information about correlations between gene frequencies in different 
#' populations (or single diploid genome samples), which can be informative 
#' about patterns of gene flow between these populations in form of admixture 
#' events. If a graph is constructed as a hypothesis for the relationship 
#' between the populations, equations for the expected values of the \eqn{f} 
#' statistics can be extracted, as functions of edge lengths -- representing 
#' genetic drift -- and admixture proportions.
#' 
#' This package provides functions for extracting these equations and for 
#' fitting them against computed \eqn{f} statistics. It does not currently 
#' provide functions for computing the \eqn{f} statistics -- for that we refer 
#' to  the \href{https://github.com/DReichLab/AdmixTools}{ADMIXTOOLS} software 
#' package.
#' 
#' @docType package
#' @name admixturegraph
NULL