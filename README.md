<!-- README.md is generated from README.Rmd. Please edit that file -->
Admixture Graph Manipulation and Fitting
========================================

The package provides functionality to analyse and test admixture graphs against the *f* statistics described in the paper Ancient Admixture in Human History, Patterson *et al.*, Genetics, Vol. 192, 1065--1093, 2012.

The *f* statistics --- *f2*, *f3*, and *f4* --- extract information about correlations between gene frequencies in different populations (or single diploid genome samples), which can be informative about patterns of gene flow between these populations in form of admixture events. If a graph is constructed as a hypothesis for the relationship between the populations, equations for the expected values of the *f* statistics can be extracted, as functions of edge lenghs --- representing genetic drift --- and admixture proportions.

This package provides functions for extracting these equations and for fitting them against computed *f* statistics. It does not currently provide functions for computing the *f* statistics --- for that we refer to the [ADMIXTOOLS](https://github.com/DReichLab/AdmixTools) software package.

Example
-------

Below is a quick example of how the package can be used. The example uses data from polar bears and brown bears with a black bear as outgroup and is taken from [Genomic evidence of geographically widespread effect of gene flow from polar bears into brown bears](http://onlinelibrary.wiley.com/doi/10.1111/mec.13038/abstract).

The BLK sample is the black bear, the PB sample is a polar bear, and the rest are brown bears.

I have taken the *f* statistics from Table 1 in the paper:

``` r
data(bears)
bears
#>      W  X      Y      Z      D Z.value
#> 1  BLK PB Sweden   Adm1 0.1258    12.8
#> 2  BLK PB  Kenai   Adm1 0.0685     5.9
#> 3  BLK PB Denali   Adm1 0.0160     1.3
#> 4  BLK PB Sweden   Adm2 0.1231    12.2
#> 5  BLK PB  Kenai   Adm2 0.0669     6.1
#> 6  BLK PB Denali   Adm2 0.0139     1.1
#> 7  BLK PB Sweden    Bar 0.1613    14.7
#> 8  BLK PB  Kenai    Bar 0.1091     8.9
#> 9  BLK PB Denali    Bar 0.0573     4.3
#> 10 BLK PB Sweden   Chi1 0.1786    17.7
#> 11 BLK PB  Kenai   Chi1 0.1278    11.3
#> 12 BLK PB Denali   Chi1 0.0777     6.4
#> 13 BLK PB Sweden   Chi2 0.1819    18.3
#> 14 BLK PB  Kenai   Chi2 0.1323    12.1
#> 15 BLK PB Denali   Chi2 0.0819     6.7
#> 16 BLK PB Sweden Denali 0.1267    14.3
#> 17 BLK PB  Kenai Denali 0.0571     5.6
#> 18 BLK PB Sweden  Kenai 0.0719     9.6
```

The `D` column is the f4(W,X;Y,Z) statistic and the `Z` column is the *Z*-values obtained from a blocked jacknife (see Patterson *et al.* for details).

From the statistics we can see that the ABC bears (Adm, Bar and Chi) are closer related to the polar bears compared to the other brown bears. The paper explains this with gene flow from polar bears into the ABC bears and going further out from there, but we can also explain this by several waves of admixture from ancestral polar bears into brown bears:

``` r
leaves <- c("BLK", "PB",
            "Bar", "Chi1", "Chi2", "Adm1", "Adm2",
            "Denali", "Kenai", "Sweden") 
inner_nodes <- c("R", "PBBB",
                 "Adm", "Chi", "BC", "ABC",
                 "x", "y", "z",
                 "pb_a1", "pb_a2", "pb_a3", "pb_a4",
                 "bc_a1", "abc_a2", "x_a3", "y_a4")

edges <- parent_edges(c(edge("BLK", "R"),
                        edge("PB", "pb_a1"),
                        edge("pb_a1", "pb_a2"),
                        edge("pb_a2", "pb_a3"),
                        edge("pb_a3", "pb_a4"),
                        edge("pb_a4", "PBBB"),
                        
                        edge("Chi1", "Chi"),
                        edge("Chi2", "Chi"),
                        edge("Chi", "BC"),
                        edge("Bar", "BC"),
                        edge("BC", "bc_a1"),
                        
                        edge("Adm1", "Adm"),
                        edge("Adm2", "Adm"),
                        
                        admixture_edge("bc_a1", "pb_a1", "ABC", "a"),
                        edge("Adm", "ABC"),
                        
                        edge("ABC", "abc_a2"),
                        admixture_edge("abc_a2", "pb_a2", "x", "b"),
                        
                        edge("Denali", "x"),
                        edge("x", "x_a3"),
                        admixture_edge("x_a3", "pb_a3", "y", "c"),
                      
                        edge("Kenai", "y"),
                        edge("y", "y_a4"),                        
                        admixture_edge("y_a4", "pb_a4", "z", "d"),
                        
                        edge("Sweden", "z"),
                        
                        edge("z", "PBBB"),
                        edge("PBBB", "R")))

bears_graph <- agraph(leaves, inner_nodes, edges)
plot(bears_graph, show_admixture_labels = TRUE)
#> fminbnd:  Exiting: Maximum number of function evaluations has been exceeded
#>          - increase MaxFunEvals option.
#>          Current function value: 3027.37262644651
```

![](README-unnamed-chunk-3-1.png)

Fitting a graph to data
-----------------------

The graph makes predictions on how the *f4* statistics should look. The graph parameters can be fit to observed statistics using the `fit_graph` function:

``` r
fit <- fit_graph(bears, bears_graph)
fit
#> 
#> Call: inner_fit_graph(data, graph, point, Z.value, concentration, optimisation_options, 
#>     parameters, iteration_multiplier, qr_tol)
#> 
#> None of the admixture proportions are properly fitted!
#> Not all of the admixture proportions are properly fitted!
#> See summary.agraph_fit for a more detailed analysis.
#> 
#> Minimal error: 12.98523
```

You can get detailsabout the fit by calling the `summary.agraph_fit` function:

``` r
summary(fit)
#> 
#> Call: inner_fit_graph(data, graph, point, Z.value, concentration, optimisation_options, 
#>     parameters, iteration_multiplier, qr_tol)
#> 
#> None of the proportions {a, b, c, d} affect the quality of the fit!
#> 
#> Optimal admixture proportions:
#>         a         b         c         d 
#> 0.3666992 0.4977105 0.9565926 0.7986799 
#> 
#> Optimal edge lengths:
#>        edge_R_BLK       edge_R_PBBB       edge_PBBB_z   edge_PBBB_pb_a4 
#>        0.00000000        0.00000000        0.00000000        0.07852837 
#>     edge_Adm_Adm1     edge_Adm_Adm2     edge_Chi_Chi1     edge_Chi_Chi2 
#>        0.00000000        0.00000000        0.00000000        0.00000000 
#>       edge_BC_Bar       edge_BC_Chi      edge_ABC_Adm    edge_ABC_bc_a1 
#>        0.00000000        0.00000000        0.00000000        0.00000000 
#>     edge_x_Denali     edge_x_abc_a2      edge_y_Kenai       edge_y_x_a3 
#>        0.00000000        0.00000000        0.00000000        0.00000000 
#>     edge_z_Sweden       edge_z_y_a4     edge_pb_a1_PB  edge_pb_a1_bc_a1 
#>        0.00000000        0.00000000        0.00000000        0.00000000 
#>  edge_pb_a2_pb_a1 edge_pb_a2_abc_a2  edge_pb_a3_pb_a2   edge_pb_a3_x_a3 
#>        0.13643125        0.00000000        0.02156832        0.00000000 
#>  edge_pb_a4_pb_a3   edge_pb_a4_y_a4     edge_bc_a1_BC   edge_abc_a2_ABC 
#>        0.04010857        0.00000000        0.00000000        0.00000000 
#>       edge_x_a3_x       edge_y_a4_y 
#>        0.00000000        0.00000000 
#> 
#> Solution to a homogeneous system of edge lengths with the optimal admixture proportions:
#> Adding any such solution to the optimal one will not affect the error.
#> 
#> Free edge lengths:
#> edge_R_BLK
#> edge_R_PBBB
#> edge_PBBB_z
#> edge_Adm_Adm1
#> edge_Adm_Adm2
#> edge_Chi_Chi1
#> edge_Chi_Chi2
#> edge_BC_Bar
#> edge_BC_Chi
#> edge_ABC_Adm
#> edge_ABC_bc_a1
#> edge_x_Denali
#> edge_x_abc_a2
#> edge_y_Kenai
#> edge_y_x_a3
#> edge_z_Sweden
#> edge_z_y_a4
#> edge_pb_a1_PB
#> edge_pb_a1_bc_a1
#> edge_pb_a2_abc_a2
#> edge_pb_a3_x_a3
#> edge_pb_a4_y_a4
#> edge_bc_a1_BC
#> edge_abc_a2_ABC
#> edge_x_a3_x
#> edge_y_a4_y
#> 
#> Bounded edge lengths:
#> edge_PBBB_pb_a4 = 0
#> edge_pb_a2_pb_a1 = 0
#> edge_pb_a3_pb_a2 = 0
#> edge_pb_a4_pb_a3 = 0
#> 
#> Minimal error:
#> 12.98523
```

You can make a plot of the fit against the data by calling the `plot.agraph_fit` function:

``` r
plot(fit)
```

![](README-unnamed-chunk-6-1.png)

The plot shows the observed *f4* statistics with error bars (in black) plus the predicted values from the graph.

The result of this is a `ggplot2` object that you can modify by adding `ggplot2` commands in the usual way.

Read the vignette `admixturegraph` for more examples.
