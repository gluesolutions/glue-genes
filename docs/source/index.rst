.. glue-genes documentation master file, created by
   sphinx-quickstart on Wed Jan 11 16:58:00 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

glue genes Documentation
###########################

`glue genes <https://github.com/gluesolutions/glue-genes>`_ is a custom version of the `glue software <http://glueviz.org>`_
for exploratory data analysis of genomics and transcriptomics data.

Features
=========

`glue genes <https://github.com/gluesolutions/glue-genes>`_ provides all
the core features of glue:

* Interactive, linked statistical graphics
* A user interface for 'glueing' together multiple datasets with complicated relationships
* A highly scriptable and extendable environment

In addition, `glue genes <https://github.com/gluesolutions/glue-genes>`_ provides:

* Data loaders for genomics data file formats such as .bed, .bigwig, .h5ad, and .loom
* Viewers including: 2D Heatmap, Small Multiples, QTL Viewer
* Menubar plug-ins to facilitate analysis of single-cell transcriptomics data

Documentation Contents
======================

.. grid:: 1 1 2 2
    :gutter: 2

    .. grid-item-card:: Getting Started
        :img-top: images/getting_started.svg
        :link: getting-started/index
        :link-type: doc

        New to *glue genes*? Use the getting started guides for help
        with installation, a brief introduction to glue, and tutorials to
        get you started.

    .. grid-item-card:: How-To Guides
        :img-top: images/how_to.svg
        :link: how-to/index
        :link-type: doc

        Are you looking for help on how to do something specific? The
        How-To Guides will explain how to perfom specific tasks using
        *glue genes*.

    .. grid-item-card::  User Guide
        :img-top: images/user_guide.svg
        :link: user-guide/index
        :link-type: doc

        The user guide provides in-depth information on the
        key concepts of *glue* and *glue genes* with useful
        background information and explanation about the
        how data is represented and displayed.


    .. grid-item-card::  API reference
        :img-top: images/api_reference.svg
        :link: api/index
        :link-type: doc

        The API reference contains a description of the main methods
        and classes in *glue genes*, along with their parameters. This
        is most useful for developers and users who wish to interact
        with *glue genes* using the built-in terminal.



.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: For New Users

   Getting Started <getting-started/index>
   How-To Guides <how-to/index>
   User Guide <user-guide/index>

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: For Advanced Users/Developers
   
   API <api/index>


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
