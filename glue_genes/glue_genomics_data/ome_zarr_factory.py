from pathlib import Path

import numpy as np
from glue.config import data_factory
from ome_zarr.io import parse_url
from ome_zarr.reader import Reader

from .multires_data import MultiResolutionData

__all__ = ["is_ome_zarr"]


def is_ome_zarr(filename, **kwargs):
    """Check if a file is a Zarr file"""
    zarr = parse_url(filename)
    if zarr:
        return True
    else:
        return False


def get_nontrivial_dimensions(im):
    """
    Get the number of nontrivial dimensions in an image.
    """
    nontrivial_dims = []
    for i, dim in enumerate(im.shape):
        if dim > 1:
            nontrivial_dims.append(i)
    return len(nontrivial_dims)


@data_factory("OME-ZARR Loader", is_ome_zarr, priority=900)
def read_ome_zarr(filename):
    """
    Read an OME-Zarr file into a glue MultiResolutionData object.

    Parameters
    ----------
    filename: str
        The pathname to the file.

    Notes
    -----
    TODO: Use metadata to set up coordinates
    This is particularly useful for giving t and z proper names in the 4/5D case

    """
    reader = Reader(parse_url(filename))
    # nodes may include images, labels etc
    nodes = list(reader())
    # first node will be the image pixel data
    image_node = nodes[0]
    meta = image_node.metadata
    axes_types = [x["name"] for x in meta["axes"]]
    print(axes_types)
    dask_data = image_node.data
    data_components = []

    # A channel axis is treated differently than the other axes
    # because we load channels as separate glue components
    # We squeeze out trivial dimensions (length 1) but we
    # need to keep track of how many this is so we know
    # which dimensions have been scaled down.
    non_trivial_dims = get_nontrivial_dimensions(dask_data[0])

    try:
        channel_index = axes_types.index("c")
        num_channels = dask_data[0].shape[channel_index]
        non_trivial_dims -= 1
    except ValueError:
        channel_index = None
        num_channels = 1

    if channel_index == 0:
        for dask_array in dask_data:
            channel_names = meta.get("name", None)
            if len(channel_names) == num_channels:
                channel_dict = {
                    ch: np.flipud(np.squeeze(dask_array[i, ...]))
                    for i, ch in enumerate(channel_names)
                }
            else:
                channel_dict = {
                    f"ch{i}": np.flipud(np.squeeze(dask_array[i, ...]))
                    for i in range(num_channels)
                }
            data_components.append(channel_dict)

    elif channel_index == 1:
        for dask_array in dask_data:
            channel_names = meta.get("name", None)
            if len(channel_names) == num_channels:
                channel_dict = {
                    ch: np.flipud(np.squeeze(dask_array[:, i, ...]))
                    for i, ch in enumerate(channel_names)
                }
            else:
                channel_dict = {
                    f"ch{i}": np.flipud(np.squeeze(dask_array[:, i, ...]))
                    for i in range(num_channels)
                }
            data_components.append(channel_dict)

    elif channel_index == 2:
        for dask_array in dask_data:
            channel_names = meta.get("name", None)
            if len(channel_names) == num_channels:
                channel_dict = {
                    ch: np.flipud(np.squeeze(dask_array[:, :, i, ...]))
                    for i, ch in enumerate(channel_names)
                }
            else:
                channel_dict = {
                    f"ch{i}": np.flipud(np.squeeze(dask_array[:, :, i, ...]))
                    for i in range(num_channels)
                }
            data_components.append(channel_dict)
    elif channel_index is None:
        for dask_array in dask_data:
            channel_dict = {}
            channel_dict["value"] = np.flipud(np.squeeze(dask_array))
            data_components.append(channel_dict)
    else:
        pass

    return MultiResolutionData(
        **data_components[0],
        all_resolutions=data_components,
        label=Path(filename).stem,
        reduced_dims=[
            non_trivial_dims - 2,
            non_trivial_dims - 1,
        ],
    )
