from pathlib import Path

import pyranges as pr
from glue.config import data_factory
from glue.core import Data

__all__ = ["read_bigwig"]


def is_bigwig(filename, **kwargs):
    """Check if a file is a BigWig file"""
    return filename.endswith(".bigwig") or filename.endswith(".bw")


@data_factory("BigWig data loader", is_bigwig, priority=999)
def read_bigwig(file_name):
    """
    Read a BigWig file into a glue Data object.

    Parameters
    ----------
    file_name: str
        The pathname to the BigWig file.

    Notes
    -----
    The ecosystem around bigwig files supports selective data loading and
    selection from the on-disk file. This will probably require creating a
    custom glue data object that supports keeping the data on disk. Currently,
    this loader loads all data into memory.

    Also of note is the relation with the bed files, i.e. that they may serve to
    identify particular areas of the BigWig files that are of interest. It may
    be necessary to create a loader that also takes an associated BED file so
    that only relevant information is loaded.

    TODO: Investigate if there is a faster library than pyranges for this.
    Perhaps pyBigWig?

    """
    bw_data = pr.read_bigwig(file_name).as_df()

    return Data(**{k: bw_data[k] for k in bw_data.columns}, label=Path(file_name).stem)
