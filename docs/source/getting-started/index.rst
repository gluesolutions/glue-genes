.. _Getting Started:

Getting Started with glue genes
################################

Installation
================

You can install glue genes as a Python package.

We recommend installing into a dedicated virtual environment using conda::

    conda create -n glue-genes-env
    conda activate glue-genes-env

    conda install -c conda-forge glue-core
    pip install git+https://github.com/gluesolutions/glue-genes.git

This process uses conda to get some of the difficult-to-install dependencies
for glue set up first.

Brief Introduction to glue
============================

glue genes is a custom version of the glue Python library to explore relationships within and among related datasets.
The main features of glue include:

* **Linked Statistical Graphics.** With glue, users can create scatter plots, histograms and images (2D and 3D) of their data. glue is focused on the brushing and linking paradigm, where selections in any graph propagate to all others.
* **Flexible linking across data.** glue uses the logical links that exist between different data sets to overlay visualizations of different data, and to propagate selections across data sets. These links are specified by the user, and are arbitrarily flexible.
* **Full scripting capability.** glue is written in Python, and built on top of its standard scientific libraries (i.e., Numpy, Matplotlib, Scipy). Users can easily integrate their own Python code for data input, cleaning, and analysis.

The following video gives a 2 minute overview of "the three meanings of glue".

.. raw:: html

   <center>
   <iframe width="560" height="315" src="https://www.youtube.com/embed/SFVnqlHPJ6M" title="YouTube video player" frameborder="0" allow="accelerometer; clipboard-write; encrypted-media; gyroscope; picture-in-picture;" allowfullscreen></iframe>
   </center>

The `full glue documentation <http://docs.glueviz.org/>`_ is a great resource for learning more about the basic and advanced features of glue.