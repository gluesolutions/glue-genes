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

See [the glue genes documentation](https://glue-genes.readthedocs.io/en/latest/how-to/installation.html) for installation help.

We recommend installing into a new virtual environment. 

`pip install glue-genes[qt6]`

### Usage

After installation, glue (with the glue genes plug-ins loaded) can be invoked at the terminal:
`glue`


### Development

glue genes is under active development. Please file all bug reports as github issues. 
