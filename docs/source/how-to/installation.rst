.. _Installation:

How do I install glue genes?
############################

The best way to install glue genes depends on your use case. If you are not
familiar with Python, the easiest way to get started is using the
:ref:`standalone-applications`. If you want to develop glue genes, or
install additional glue plugins, you should follow the instructions for
:ref:`python-installation`.


.. _standalone-applications:

Standalone Applications
========================

The easiest way to get started with glue genes is to download and run these
pre-built single-file applications.

.. warning::
    Currently, these applications will prompt security warnings and apparently fail to work
    on most operating systems the first time you launch them. The core glue documentation
    shows you how to work around this problem. See `Mac instructions <http://docs.glueviz.org/en/stable/installation/standalone.html#macos-x>`_ and `Windows instructions <http://docs.glueviz.org/en/stable/installation/standalone.html#windows>`_ for detailed walkthroughs (with
    screenshots) on how to fix this.

* `Mac <https://gluesolutions.s3.amazonaws.com/installers/genes/glue+genes.dmg>`_
* `Windows <https://gluesolutions.s3.amazonaws.com/installers/genes/glue+genes.exe>`_
* `Linux <https://gluesolutions.s3.amazonaws.com/installers/genes/glue-genes>`_

These applications contain the latest stable version of glue, along with the glue genes
plug-ins and other core plug-ins such as 3D viewers, export to plotly, and WWT.

.. note::
    The standalone Mac application will probably launch with two icons in the dock;
    there are not actually two instances of glue running. We hope to fix this soon.


.. note::
    Using the standalone applications does not allow you to install additional glue plug-ins
    or upgrade/modify glue code. If you want to develop glue genes or make use of any of the 
    myriad ways you can `customize glue <https://glueviz.readthedocs.io/en/stable/customizing_guide/customization.html>`_ you should follow the 
    :ref:`python-installation` instructions.


.. _python-installation:

Python Installation
====================

Installing glue genes as a Python package allows you to modify glue genes,
install additional `glue plugins <https://glueviz.org/plugins.html>`_, and
take full advantage of the many ways you can `customize <http://docs.glueviz.org/en/stable/customizing_guide/customization.html>`_ glue/glue genes to suit
you needs.

We recommend installing into a dedicated virtual environment using `conda <https://www.anaconda.com>`_::

    conda create -n glue-genes-env
    conda activate glue-genes-env

    conda install -c conda-forge glue-core
    pip install glue-genes

This process uses conda to get some of the difficult-to-install dependencies
for glue installed first. If you prefer not to use conda (or the alternative `mamba <https://mamba.readthedocs.io/en/latest/>`_) then please consult the core glue documentation for `alternative methods for installing glue <http://docs.glueviz.org/en/stable/installation/installation.html>`_, and afterwards simply::

    pip install glue-genes
