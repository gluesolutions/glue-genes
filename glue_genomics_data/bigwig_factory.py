from glue.config import data_factory
from glue.core import Data
import pyranges as pr
from pathlib import Path

__all__ = ['is_bigwig', 'read_bigwig']


def is_bigwig(filename, **kwargs):
    return filename.endswith('.bigwig') or filename.endswith('.bw')


@data_factory('BigWig data loader', is_bigwig, priority=999)
def read_bigwig(file_name):
    """
    Read a bigwig file into glue.

    Notes
    -----
    The ecosystem around bigwig files supports selective data loading and
    selection from the on-disk file. This will probably require creating a 
    custom glue data object that supports keeping the data on disk. Currently,
    this loader loads all data into memory.

    Also of note is the relation with the bed files, i.e. that they serve to
    identify particular areas of the bigwig files that are of interest. It may
    be necessary to create a loader that also takes an associated bed file so
    that only relevant information is loaded.
    """
    # TODO: pyranges seem slow, switch to pyBigWig?
    bw_data = pr.read_bigwig(file_name).as_df()

    return Data(**{k: bw_data[k] for k in bw_data.columns}, 
                label=Path(file_name).stem)