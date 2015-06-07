<!-- README.md is generated from README.Rmd. Please edit that file -->
admixturegraph: Admixture Graph Manipulation and Fitting
========================================================

The package provides functionality to analyse and test admixture graphs against the f statisticsdescribed in the paper [Ancient Admixture in Human History](http://tinyurl.com/o5a4kr4), Patternson et al., Genetics, Vol. 192, 1065--1093, 2012.

The f statistics --- \(f_2\), \(f_3\), and \(f_4\) --- extract information about correlations between gene frequencies in different populations (or single diploid genome samples), which can be informative about patterns of gene flow between these populations in form of admixture events. If a graph is constructed as a hypothesis for the relationship between the populations, equations for the expected values of the f statistics can be extracted, as functions of edge lenghs — representing genetic drift — and admixture proportions.

This package provides functions for extracting these equations and for fitting them against computed f statistics. It does not currently provide functions for computing the f statistics — for that we refer to the [ADMIXTOOLS](https://github.com/DReichLab/AdmixTools) software package.
