from pathlib import Path

import dask.array as da
import numpy as np
import zarr
from glue.config import data_factory

from .multires_data import MultiResolutionData

try:
    from napari_lazy_openslide import OpenSlideStore
    from openslide import (
        PROPERTY_NAME_COMMENT,
        OpenSlide,
        OpenSlideUnsupportedFormatError,
    )

    OPENSLIDE_INSTALLED = True

except ImportError:
    OPENSLIDE_INSTALLED = False


__all__ = ["is_openslide"]


def is_openslide(filename, **kwargs):
    """Check if a file can be read with openslide"""
    if not OPENSLIDE_INSTALLED:
        return False
    if OpenSlide.detect_format(filename) is None:
        return False
    try:
        slide = OpenSlide(filename)
    except OpenSlideUnsupportedFormatError:
        return False

    description = slide.properties.get(PROPERTY_NAME_COMMENT)
    # Don't try to handle OME-TIFF
    # https://github.com/cgohlke/tifffile/blob/b346e3bd7de81de512a6715b01124c8f6d60a707/tifffile/tifffile.py#L5781
    if description and description[-4:] == "OME>":
        return False

    if slide.level_count == 1:
        return None

    slide.close()
    return True


@data_factory("OpenSlide Loader", is_openslide, priority=999)
def read_open_slide_as_zarr(filename):
    """
    Read an OpenSlide compatible image using napari-lazy-openslide
    to make it look like a OME-Zarr

    Parameters
    ----------
    filename: str
        The pathname to the file.

    Notes
    -----
    """
    store = OpenSlideStore(filename)
    grp = zarr.open(store, mode="r")
    multiscales = grp.attrs["multiscales"][0]
    dask_data = [
        da.from_zarr(store, component=d["path"]) for d in multiscales["datasets"]
    ]
    data_components = []
    for dask_array in dask_data:
        channel_names = ["red", "green", "blue", "opacity"]
        channel_dict = {
            ch: np.flipud(np.squeeze(dask_array[..., i]))
            for i, ch in enumerate(channel_names)
        }
        data_components.append(channel_dict)

    # add_kwargs = {"name": multiscales["name"]}
    return MultiResolutionData(
        **data_components[0],
        all_resolutions=data_components,
        label=Path(filename).stem,
    )
