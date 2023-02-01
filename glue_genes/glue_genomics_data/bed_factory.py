from pathlib import Path

import pandas as pd
from glue.config import data_factory
from glue.core import Data

__all__ = ["read_bed"]


def is_bed(filename, **kwargs):
    """
    Check if a file is a BED file.
    """
    return filename.endswith(".bed")


def remap_columns(x):
    """
    Rename the first three columns in a BED
    file to be: chr, start, end.

    Remaining columns are named ucol (for
    unknown column) to make sure they show
    up as later components in all viewers
    """
    if x == 0:
        return "chr"
    elif x == 1:
        return "start"
    elif x == 2:
        return "end"
    else:
        return "ucol" + str(x)


@data_factory("BED data loader", is_bed, priority=999)
def read_bed(file_name):
    """
    Read a BED file into a glue Data object.

    The first three columns will be named: chr, start, end

    The remaining columns are named ucol (for unknown column)
    and can be renamed by the user.

    Parameters
    ----------
    file_name: str
        The pathname to the BED file.
    """
    bed_data = pd.read_csv(file_name, sep="\t", header=None)
    bed_data.rename(
        columns=remap_columns, inplace=True
    )  # Rename first three columns as standard

    return Data(
        **{k: bed_data[k] for k in bed_data.columns}, label=Path(file_name).stem
    )
