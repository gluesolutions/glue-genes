=============================
Multi-Resolution Imaging Data
=============================

.. currentmodule:: glue_genes.glue_genomics_data.multires_data

.. automodule:: glue_genes.glue_genomics_data.multires_data


Multi-resolution images are handled with a special sub-class of the glue Data object that uses the full-resolution image as a main components, and then creates ReducedResolutionData objects for all the downsampled images. The MultiResolutionData object then passes off calls to `get_data()` to the appropriate ReducedResolutionData object.

MultiResolutionData Object
----------------------------

.. autosummary::
   :toctree: api/
   
   MultiResolutionData

ReducedResolutionData Object
-----------------------------


.. autosummary::
   :toctree: api/

   ReducedResolutionData