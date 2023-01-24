.. _Cells to Genes:

How do I connect cells and expression data?
###########################################

You want to compute a summary statistic for your cells from a subset of genes, or you want to identify genes that are differentially expressed between two subsets of cells. Because cells and genes are only connected through their expression matrix, glue genes needs to do a calculation to make these kinds of connections, and you have to use a dedicated plug-in to tell glue how to do this calculation.

.. _Gene Expression Summary:

Create a summary of gene expression over cells
************************************************

#. Create a subset over some genes of interest. This is just a starting point; the subset can subsequently be updated. This subset must be defined on the genes in the single cell data set, either because it is a direct selection of these genes or because it is a selection on another dataset that is appropriately joined to the gene expression data in the single cell data set.

   .. figure:: images/connect_create_gene_subset.png
      :align: center
      :width: 66%
    
      Create a subset over the genes in your single cell data set.


#. Choose the **Calculate Summary Over Gene Subset** menu item from the **Plugins** menubar item. In the dialog window that pops up, make sure that the correct dataset for the cells is selected, as well as the subset that you just defined over genes. Finally, choose the summary statistic you want to use. *Mean Value* is fastest if you have very large datasets, but *Module Score* or *PCA* do a better job of specifically highlighting the cells where the genes are differentially expressed.

   .. figure:: images/connect_choose_plugin.png
      :align: center
      :width: 66%
    
      Choose the Calculate Summary Over Gene Subset menubar item.


   .. figure:: images/connect_summary_dialog.png
      :align: center
      :width: 66%
    
      Choose the parameters for the summary calculation.


#. A dialog will tell you that a new component has been added to your cells dataset. This new component contains the summary statistic for the gene subset, and will update as the gene subset changes. You can now either use this new component as a Linear color-scale in a visualization of the cells or plot this new component on one axis of another Viewer. Now explore what happens to these views as you update the gene subset.

.. warning::
    Changing the name of the gene subset after creation will not update the name of the component in the cells dataset, which can be confusing. 

.. _DGE Between Cell Subsets:

Identify differentially expressed genes between two subsets of cells
*********************************************************************

#. Create the two subsets of cells you want to compare. This plug-in does not live-update the set of differentially expressed genes, so you need to decide on the subsets you want to compare first.

#. Choose the **Scanpy Differential Gene Expression** menu item from the **Plugins** menubar item. In the dialog window that pops up, make sure to choose the X matrix (gene expression matrix) for these subsets of cells, identify the component that contains the gene names or IDs, and select the two subsets of cells as Subset 1 and Subset 2. This plug-in will return a ranking for all the genes, but will automatically create a subset for the top N genes; choose a different value of N if you desire.

   .. figure:: images/connect_dge_menubar.png
      :align: center
      :width: 66%
    
      Choose the Scanpy Differential Gene Expression menubar item.

   .. figure:: images/connect_dge_dialog.png
      :align: center
      :width: 66%
    
      Set up the appropriate datasets and parameters for running differential gene expression.

#. Refine and visualize the subset and dataset created. This plugin creates a new dataset called **DGE between Subset 1 and Subset 2** (where Subset 1 and Subset 2 are the names of your subsets over cells) which is linked to the gene metadata in your single-cell data set. It also creates a subset with the same name corresponding to the top N genes. This subset will be displayed on visualizations or table of genes, and you can use the full dataset to refine this subset for futher analysis or export into other analysis programs.