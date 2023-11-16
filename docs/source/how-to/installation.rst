.. _Installation:

How do I install glue genes?
############################

.. _python-installation:

Python Installation
====================

Installing glue genes as a Python package allows you to modify glue genes,
install additional `glue plugins <https://glueviz.org/plugins.html>`_, and
take full advantage of the many ways you can `customize <http://docs.glueviz.org/en/stable/customizing_guide/customization.html>`_ glue/glue genes to suit
you needs.

We recommend installing glue genes into a dedicated virtual environment using pip::

    python -m venv glue-genes-env
    source glue-genes-env/bin/activate

    pip install glue-genes[qt6]

If you encounter difficulties, then please consult the core glue documentation for `alternative methods for installing glue <http://docs.glueviz.org/en/stable/installation/installation.html>`_, and afterwards simply::

    pip install glue-genes[qt6]
