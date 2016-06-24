# Version 1.0.0.9000

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
 
