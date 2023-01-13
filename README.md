![glue_genes_logo_stacked](https://user-images.githubusercontent.com/3639698/137145077-2c2c9011-68bd-4770-9d58-494bf7632a33.png)


The glue-genes meta-package
===========================

[glue genes](https://github.com/gluesolutions/glue-genes) is a custom version of the [glue software](https://glueviz.org) for the 
visual exploration of genomics data, developed in partnership with [The Jackson Laboratory](https://jax.org).

### Features

glue genes provides all the core features of glue:

* Interactive, linked statistical graphics
* A user interface for 'glueing' together multiple datasets with complicated relationships
* A highly scriptable and extendable environment

In addition, glue genes provides:

* Data loaders for genomics data file formats such as .bed, .bigwig, .h5ad, and .loom
* Viewers including: 2D Heatmap, Small Multiples, QTL Viewer
* Menubar plug-ins to facilitate analysis of single-cell transcriptomics data

### Install

This package does not contain any code, but instead is a meta-package that 
installs all the packages necessary to use glue with genomics data. 

Due to some complicated binary dependencies, the recommended procedure for 
installing this package is to use [conda](https://www.anaconda.com) and install into a new environment as follows:

```
conda create -n glue-genes python==3.9
conda activate glue-genes

conda install glue-core
pip install glue-genes@git+https://github.com/gluesolutions/glue-genes.git
```

### Usage

The Scanpy differential gene expression plugin requires some custom actions to run at startup, which can be achieved by starting glue using this command:

`glue --startup=setup_anndata`



### Development

glue genes is currently under active development towards a 1.0 release. Please file all bug reports as github issues. 
