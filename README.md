<!-- README.md is generated from README.Rmd. Please edit that file -->
Admixture Graph Manipulation and Fitting
========================================

The package provides functionality to analyse and test admixture graphs against the f statisticsdescribed in the paper [Ancient Admixture in Human History](http://tinyurl.com/o5a4kr4), Patternson *et al.*, Genetics, Vol. 192, 1065--1093, 2012.

The f statistics --- f2, f3, and f4 --- extract information about correlations between gene frequencies in different populations (or single diploid genome samples), which can be informative about patterns of gene flow between these populations in form of admixture events. If a graph is constructed as a hypothesis for the relationship between the populations, equations for the expected values of the f statistics can be extracted, as functions of edge lenghs — representing genetic drift — and admixture proportions.

This package provides functions for extracting these equations and for fitting them against computed f statistics. It does not currently provide functions for computing the f statistics — for that we refer to the [ADMIXTOOLS](https://github.com/DReichLab/AdmixTools) software package.

Example
-------

Below is a quick example of how the package can be used. The example uses data from polar bears and brown bears with a black bear as outgroup and is taken from [Genomic evidence of geographically widespread effect of gene flow from polar bears into brown bears](http://onlinelibrary.wiley.com/doi/10.1111/mec.13038/abstract).

The BLK sample is the black bear, the PB sample is a polar bear, and the rest are brown bears.

I have taken the D statistics from Table 1 in the paper and have the statistics:

``` r
data(bears)
bears
#>      W  X      Y      Z       D Z.value
#> 1  BLK PB Sweden   Adm1  0.1258    12.8
#> 2  BLK PB  Kenai   Adm1  0.0685     5.9
#> 3  BLK PB Denali   Adm1  0.0160     1.3
#> 4  BLK PB Sweden   Adm2  0.1231    12.2
#> 5  BLK PB  Kenai   Adm2  0.0669     6.1
#> 6  BLK PB Denali   Adm2  0.0139     1.1
#> 7  BLK PB Sweden    Bar  0.1613    14.7
#> 8  BLK PB  Kenai    Bar  0.1091     8.9
#> 9  BLK PB Denali    Bar  0.0573     4.3
#> 10 BLK PB Sweden   Chi1  0.1786    17.7
#> 11 BLK PB  Kenai   Chi1  0.1278    11.3
#> 12 BLK PB Denali   Chi1  0.0777     6.4
#> 13 BLK PB Sweden   Chi2  0.1819    18.3
#> 14 BLK PB  Kenai   Chi2  0.1323    12.1
#> 15 BLK PB Denali   Chi2  0.0819     6.7
#> 16 BLK PB Sweden Denali  0.1267    14.3
#> 17 BLK PB  Kenai Denali  0.0571     5.6
#> 18 BLK PB Sweden  Kenai  0.0719     9.6
#> 19 BLK PB Denali  Kenai -0.0571    -5.6
```

The D column is the f4(W,X;Y,Z) statistic and the Z column is the Z-values obtained from a blocked jacknife (see Patterson *et al.* for details).

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
                        
                        admixture_edge("bc_a1", "pb_a1", "ABC"),
                        edge("Adm", "ABC"),
                        
                        edge("ABC", "abc_a2"),
                        admixture_edge("abc_a2", "pb_a2", "x"),
                        
                        edge("Denali", "x"),
                        edge("x", "x_a3"),
                        admixture_edge("x_a3", "pb_a3", "y"),
                      
                        edge("Kenai", "y"),
                        edge("y", "y_a4"),                        
                        admixture_edge("y_a4", "pb_a4", "z"),
                        
                        edge("Sweden", "z"),
                        
                        edge("z", "PBBB"),
                        edge("PBBB", "R")))
 

admixtures <- admixture_proportions(c(admix_props("bc_a1", "pb_a1", "ABC", "a"),
                                      admix_props("abc_a2", "pb_a2", "x", "b"),
                                      admix_props("x_a3", "pb_a3", "y", "c"),
                                      admix_props("y_a4", "pb_a4", "z", "d")))
                                
bears_graph <- agraph(leaves, inner_nodes, edges, admixtures)
plot(bears_graph, show_inner_node_labels = TRUE, show_admixture_labels = TRUE)
```

![](README-graph-1.png)

The graph makes predictions on how the f4 statistics should look, in particular it allows us to predict the signs of the f4 statistics.

``` r
add_graph_f4_sign(bears, bears_graph)
#>      W  X      Y      Z       D Z.value graph_f4_sign
#> 1  BLK PB Sweden   Adm1  0.1258    12.8             1
#> 2  BLK PB  Kenai   Adm1  0.0685     5.9             1
#> 3  BLK PB Denali   Adm1  0.0160     1.3             1
#> 4  BLK PB Sweden   Adm2  0.1231    12.2             1
#> 5  BLK PB  Kenai   Adm2  0.0669     6.1             1
#> 6  BLK PB Denali   Adm2  0.0139     1.1             1
#> 7  BLK PB Sweden    Bar  0.1613    14.7             1
#> 8  BLK PB  Kenai    Bar  0.1091     8.9             1
#> 9  BLK PB Denali    Bar  0.0573     4.3             1
#> 10 BLK PB Sweden   Chi1  0.1786    17.7             1
#> 11 BLK PB  Kenai   Chi1  0.1278    11.3             1
#> 12 BLK PB Denali   Chi1  0.0777     6.4             1
#> 13 BLK PB Sweden   Chi2  0.1819    18.3             1
#> 14 BLK PB  Kenai   Chi2  0.1323    12.1             1
#> 15 BLK PB Denali   Chi2  0.0819     6.7             1
#> 16 BLK PB Sweden Denali  0.1267    14.3             1
#> 17 BLK PB  Kenai Denali  0.0571     5.6             1
#> 18 BLK PB Sweden  Kenai  0.0719     9.6             1
#> 19 BLK PB Denali  Kenai -0.0571    -5.6            -1
```

The way the signs are predicted is by extracting the equations for the f4 statistics that the graph implies: For each quartet of leaves we can extract an equation for the corresponding f4 statistics --- an equation in the edge lenghts and admixture proportions --- and if this equation only have positive values we know that the sign must be positive, if it only has negative values we know that it must be negative, and if it constant zero we know it must be zero.

In general we will not always have only positive or negative terms, in which case we cannot this simply predict the sign for f4 statistics. If this is the case we need to set the parameters of the graph --- the edge lengths and admixture proportions --- to get the sign, and in that case we can also predict the numerical value of the f4 statistics from the graph.

Fitting a graph to data
-----------------------

If you have the *neldermead* package installed you can also fit graph parameters to data. This is done using the *fit\_graph* function

``` r
fit <- fit_graph(bears, bears_graph, optimset(MaxFunEvals=10000))
#> fminbnd:  Exiting: Maximum number of function evaluations has been exceeded
#>          - increase MaxFunEvals option.
#>          Current function value: 0.00295657755132507
fit
#> Call:
#> fit_graph(bears, bears_graph, optimset(MaxFunEvals = 10000))
#> 
#> Sum of squared error: 0.002956578
```

The object it returns contains an environment that contains the fitted parameters and a data frame containing the original data together with an extra column, graph\_f4, containing the fitted values.

You can get the fitted values by calling the *summary* function.

``` r
summary(fit)
#> $edges
#>        edge_R_BLK       edge_R_PBBB       edge_PBBB_z   edge_PBBB_pb_a4 
#>        0.19057921        0.00000100        0.71193890        0.12445339 
#>     edge_Adm_Adm1     edge_Adm_Adm2     edge_Chi_Chi1     edge_Chi_Chi2 
#>        0.83561262        0.12008785        0.73746522        0.10110172 
#>       edge_BC_Bar       edge_BC_Chi      edge_ABC_Adm    edge_ABC_bc_a1 
#>        0.99999900        0.47150452        0.00000100        0.00000100 
#>     edge_x_Denali     edge_x_abc_a2      edge_y_Kenai       edge_y_x_a3 
#>        0.40288193        0.86128120        0.25411597        0.00000100 
#>     edge_z_Sweden       edge_z_y_a4     edge_pb_a1_PB  edge_pb_a1_bc_a1 
#>        0.35447815        0.98815364        0.74276859        0.73857230 
#>  edge_pb_a2_pb_a1 edge_pb_a2_abc_a2  edge_pb_a3_pb_a2   edge_pb_a3_x_a3 
#>        0.02048059        0.53348206        0.09907696        0.84540469 
#>  edge_pb_a4_pb_a3   edge_pb_a4_y_a4     edge_bc_a1_BC   edge_abc_a2_ABC 
#>        0.00000100        0.49357293        0.45251097        0.90377890 
#>       edge_x_a3_x       edge_y_a4_y 
#>        0.00000100        0.09023466 
#> 
#> $admixture_proportions
#>         a         b         c         d 
#> 0.5405489 0.1558058 0.7111461 0.4715302
```

This function also returns the fitted values as a list, so you can assign the result to an object if you need to access it later.

You can also get the fitted parameters using the generic *coef* or *coefficients* funcions

``` r
coef(fit)
#>        edge_R_BLK       edge_R_PBBB       edge_PBBB_z   edge_PBBB_pb_a4 
#>        0.19057921        0.00000100        0.71193890        0.12445339 
#>     edge_Adm_Adm1     edge_Adm_Adm2     edge_Chi_Chi1     edge_Chi_Chi2 
#>        0.83561262        0.12008785        0.73746522        0.10110172 
#>       edge_BC_Bar       edge_BC_Chi      edge_ABC_Adm    edge_ABC_bc_a1 
#>        0.99999900        0.47150452        0.00000100        0.00000100 
#>     edge_x_Denali     edge_x_abc_a2      edge_y_Kenai       edge_y_x_a3 
#>        0.40288193        0.86128120        0.25411597        0.00000100 
#>     edge_z_Sweden       edge_z_y_a4     edge_pb_a1_PB  edge_pb_a1_bc_a1 
#>        0.35447815        0.98815364        0.74276859        0.73857230 
#>  edge_pb_a2_pb_a1 edge_pb_a2_abc_a2  edge_pb_a3_pb_a2   edge_pb_a3_x_a3 
#>        0.02048059        0.53348206        0.09907696        0.84540469 
#>  edge_pb_a4_pb_a3   edge_pb_a4_y_a4     edge_bc_a1_BC   edge_abc_a2_ABC 
#>        0.00000100        0.49357293        0.45251097        0.90377890 
#>       edge_x_a3_x       edge_y_a4_y                 a                 b 
#>        0.00000100        0.09023466        0.54054887        0.15580580 
#>                 c                 d 
#>        0.71114614        0.47153019
```

To get the fitted predictions, together with the data used for fitting, use the *fitted* function.

``` r
fitted(fit)
#>      W  X      Y      Z       D Z.value    graph_f4
#> 1  BLK PB Sweden   Adm1  0.1258    12.8  0.12385302
#> 2  BLK PB  Kenai   Adm1  0.0685     5.9  0.06516948
#> 3  BLK PB Denali   Adm1  0.0160     1.3  0.01839679
#> 4  BLK PB Sweden   Adm2  0.1231    12.2  0.12385302
#> 5  BLK PB  Kenai   Adm2  0.0669     6.1  0.06516948
#> 6  BLK PB Denali   Adm2  0.0139     1.1  0.01839679
#> 7  BLK PB Sweden    Bar  0.1613    14.7  0.18880479
#> 8  BLK PB  Kenai    Bar  0.1091     8.9  0.13012126
#> 9  BLK PB Denali    Bar  0.0573     4.3  0.08334856
#> 10 BLK PB Sweden   Chi1  0.1786    17.7  0.18880479
#> 11 BLK PB  Kenai   Chi1  0.1278    11.3  0.13012126
#> 12 BLK PB Denali   Chi1  0.0777     6.4  0.08334856
#> 13 BLK PB Sweden   Chi2  0.1819    18.3  0.18880479
#> 14 BLK PB  Kenai   Chi2  0.1323    12.1  0.13012126
#> 15 BLK PB Denali   Chi2  0.0819     6.7  0.08334856
#> 16 BLK PB Sweden Denali  0.1267    14.3  0.10545623
#> 17 BLK PB  Kenai Denali  0.0571     5.6  0.04677269
#> 18 BLK PB Sweden  Kenai  0.0719     9.6  0.05868353
#> 19 BLK PB Denali  Kenai -0.0571    -5.6 -0.04677269
```

You can make a plot of the fit against the data using the *plot* function.

``` r
plot(fit)
```

![](README-fitted_data-1.png)

The plot shows the data f4 statistics with error bars (in black) plus the predicted values from the graph.

The result of this is a ggplot2 object that you can modify by adding ggplot2 commands in the usual way.
