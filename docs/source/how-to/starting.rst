.. _Starting:

How do I start glue genes?
####################################

If you have install glue as a Python package you need to start it from a terminal by typing::

    glue --startup=setup_anndata

.. note::

   The very first time glue genes launches after installation (or an update) it can take quite a long time (a minute or more). Subsequent launches are fast.

If you forget the ``--startup=setup_anndata``, glue genes will still run, but a few more advanced features related to loading and displaying single cell data will not work.

If something does not seem to be working correctly with glue genes (typically some of data loaders or viewers you expect are not present) try launching glue genes with the ``-v`` flag::

    glue -v --startup=setup_anndata

and checking the terminal output for error messages. 

If, instead, you have installed the :ref:`standalone-applications` then you can simply double-click
the application icon to launch glue. Note that you may have to deal with the security prompts if you
did not deal with them during installation. See `Mac security instructions <http://docs.glueviz.org/en/stable/installation/standalone.html#macos-x>`_ and `Windows security instructions <http://docs.glueviz.org/en/stable/installation/standalone.html#windows>`_.


What Next?
**********
After launching glue genes you probably want to :ref:`load some data for analysis<Get Data In>`.