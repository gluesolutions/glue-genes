![glue_genes_logo_stacked](https://user-images.githubusercontent.com/3639698/137145077-2c2c9011-68bd-4770-9d58-494bf7632a33.png)


The glue-genes meta-package
===========================

[glue genes](https://github.com/gluesolutions/glue-genes) is a custom version of the [glue software](https://glueviz.org) for the 
visual exploration of genomics data.

### Features

glue genes provides all the core features of glue:

* Interactive, linked statistical graphics
* A user interface for 'glueing' together multiple datasets with complicated relationships
* A highly scriptable and extendable environment

In addition, glue genes provides:

* Data loaders for genomics data file formats such as .bed, .bedgraph, .bedpe, .bigwig, and .hic
* Viewers including: Genome Track Viewer, 2D Heatmaps, Dendrogram/Tree Viewer, Network Viewer
* Changes to core-glue to enable easy exploration of genetics data

### Install

This package does not contain any code, but instead is a meta-package that 
installs all the packages necessary to use glue with genomics data. 

Due to some complicated binary dependencies, the recommended procedure for 
installing this package is to use [conda](https://www.anaconda.com) and intall into a new environment as follows:

```
conda create -n glue-genes python==3.8
conda activate glue-genes
conda install -c glueviz glueviz
conda install -c bioconda pairix
conda install -c bioconda tabix
pip install glue-genes@git+https://github.com/gluesolutions/glue-genes.git
```

### Documention

The documentation for glue genes currently lives at https://docs.google.com/document/d/13rQj2H5ZOq5-6dIJf2My9Lkbvov9xwHMrwMXWNZqA94/edit# and a PDF version is available with this repository. 

### Development

glue genes is currently under active development towards a 1.0 release. Please file all bug reports as github issues. 
