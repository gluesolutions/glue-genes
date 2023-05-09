from pathlib import Path

import zarr
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
    zarr_data = zarr.open(file_name, mode="r").image[:, :].squeeze()

    ch_0 = zarr_data[:, :, 0]
    ch_1 = zarr_data[:, :, 1]
    ch_2 = zarr_data[:, :, 2]

    return Data(ch0=ch_0, ch1=ch_1, ch2=ch_2, label=Path(file_name).stem)
