# Version 1.0.2

 * Fixed a bug in the MCMC code.
 * Introduced admixture graph library in the data folder in the form of matrices, and new
   functions graph_to_vector() and vector_to_graph() to interpret the matrix data.
 * Added rename_nodes() for renaming graph nodes in user defined or automised manner.
 * Added canonise_graph() and remove_duplicates() for detecting isomorphism among graphs
   with labeled leaves.
 * Added all_graphs() although the data folder already contains a big admixture graph library.
 * New parallel computing fitting function fit_graph_list() that takes one list and data
   as an input instead of a list and permutations.
 * Updated the vignette.

# Version 1.0.1

 * Updated code for laying out graphs with nodes that have more than two children.
 * Added options for fitting plots to make them in grayscale (Issue #11).
 * Updated functions for plotting the cost function and number of fits as functions of admixture
   proportions (Issue #12).
 * Changed specification of admixture proportions (Iussue #9) to avoid specifying redundant
   information. The changes are still backward compatible, so old analysis scripts should
   still work.

# Version 1.0.0

 * A version ready for publication.

# Version 0.4.1

 * Fixed bug that caused crossing edges for admixture events in plotting (Issue #3).
 * Added bears to data/
 * Added syntactic suggar to graph construction
 * Fix plotting bug (Issue #4)
 * Wrote an interface to fitted data objects using the generic functions usually used
   for fitted models.
 * Added code for plotting fitted data
 * Added predict() function.
 * Added function for filtering data based on the leaves in a tree and for mapping
   between samples and populations (leaves).
 * New optimisation code based on solving linear systems for those variables that are linear.
 
