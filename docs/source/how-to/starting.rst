.. _Starting:

How do I start glue genes?
####################################

Once you have :ref:`installed<Installation>` glue genes, you need to start it from
a terminal by typing::

    glue --startup=setup_anndata

.. note::

   The very first time glue genes launches after installation (or an update)
   it can take quite a long time (a minute or more). Subsequent launches are fast.

If you forget the ``--startup=setup_anndata``, glue genes will still run, but
a few more advanced features related to loading and displaying single cell data
will not work.

If something does not seem to be working correctly with glue genes (typically some of
data loaders or viewers you expect are not present) try launching glue genes with the
``-v`` flag::

    glue -v --startup=setup_anndata

and checking the terminal output for error messages. 


What Next?
**********
After launching glue genes you probably want to :ref:`load some data for analysis<Get Data In>`.