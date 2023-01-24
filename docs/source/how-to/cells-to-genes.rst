.. _Cells to Genes:

How do I connect cells and expression data?
###########################################

You want to compute a summary statistic for your cells from a subset of genes, or you want to identify genes that are differentially expressed between two subsets of cells. Because cells and genes are only connected through their expression matrix, glue genes needs to do a calculation to make these kinds of connections.

Create a summary of gene expression over cells
************************************************

#. Create a subset over some genes of interest. This is just a starting point; the subset can subsequently be updated. This subset must be defined on the genes in the single cell data set, either because it is a direct selection of these genes or because it is a selection on another dataset that is appropriately joined to the gene expression data in the single cell data set.
#. 