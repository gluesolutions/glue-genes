from pathlib import Path

import numpy as np
from dask.array import from_zarr
from glue.config import data_factory
from glue.core import Data

__all__ = ["read_zarr"]


def is_zarr(filename, **kwargs):
    """
    Check if a file is a Zarr file.
    """
    return filename.endswith(".zarr")


@data_factory("Zarr data loader", is_zarr, priority=999)
def read_zarr(file_name):
    """
    Read a Zarr file into a glue Data object.

    This works on a sample Zarr from this tutorial
    https://squidpy.readthedocs.io/en/stable/auto_tutorials/tutorial_visium_hne.html
    but is probably not general.

    The OME-NGFF stuff here:
    https://ome-zarr.readthedocs.io/en/stable/python.html

    Is probably a better target for full support.

    Parameters
    ----------
    file_name: str
        The pathname to the input file.
    """

    # zarr_data = zarr.open(file_name, mode="r").image[:, :].squeeze()
    dask_array = from_zarr(file_name, component="image").squeeze()
    channel_index = 2  # Don't know if this will generally by true
    num_channels = dask_array.shape[channel_index]
    if channel_index == 2:
        channel_dict = {
            f"ch{i}": np.flipud(np.squeeze(dask_array[:, :, i]))
            for i in range(num_channels)
        }

    return Data(**channel_dict, label=Path(file_name).stem)
