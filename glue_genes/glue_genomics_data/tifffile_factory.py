from pathlib import Path

import dask.array as da
import numpy as np
import zarr
from glue.config import data_factory
from tifffile import imread
from .multires_data import MultiResolutionData
from glue.core.data import Data


def is_tifffile(filename):
    """
    Check if a file can be opened with tifffile
    """
    try:
        imread(filename, aszarr=True)
    except Exception:
        return False
    return True


@data_factory("Tifffile Loader", is_tifffile, priority=100)
def read_tifffile_as_zarr(filename):
    """
    Read TIFF-like files as Zarr using tifffile

    Parameters
    ----------
    filename: str
        The pathname to the file.

    Notes
    -----
    """

    store = imread(filename, aszarr=True)
    grp = zarr.open(store, mode="r")

    try:
        multiscales = grp.attrs["multiscales"][0]

        dask_data = [
            da.from_zarr(store, component=d["path"]) for d in multiscales["datasets"]
        ]
        data_components = []
        for dask_array in dask_data:
            if dask_array.shape[-1] == 4:
                channel_names = ["red", "green", "blue", "opacity"]
            elif dask_array.shape[-1] == 3:
                channel_names = ["red", "green", "blue"]
            else:
                channel_names = ["ch{i}" for i in range(dask_array.shape[-1])]
            channel_dict = {
                ch: np.flipud(np.squeeze(dask_array[..., i]))
                for i, ch in enumerate(channel_names)
            }
            data_components.append(channel_dict)

        return MultiResolutionData(
            **data_components[0],
            all_resolutions=data_components,
            label=Path(filename).stem,
            reduced_dims=[  # This assumes that the image is two-dimensional
                0,
                1,
            ],
        )
    except KeyError:
        dask_array = da.from_zarr(store)
        if dask_array.shape[-1] == 4:
            channel_names = ["red", "green", "blue", "opacity"]
        elif dask_array.shape[-1] == 3:
            channel_names = ["red", "green", "blue"]
        else:
            channel_names = ["ch{i}" for i in range(dask_array.shape[-1])]
        channel_dict = {
            ch: np.flipud(np.squeeze(dask_array[..., i]))
            for i, ch in enumerate(channel_names)
        }
        return Data(**channel_dict, label=Path(filename).stem)
