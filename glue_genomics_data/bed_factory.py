from glue.config import data_factory
from glue.core import Data
import pandas as pd
from pathlib import Path
from qtpy.QtWidgets import QDialog
#from .preprocessors.peak_correlations import PeakCorrelationsPreprocessor


__all__ = ['is_bed', 'read_bed']


def is_bed(filename, **kwargs):
    return filename.endswith('.bed')

def remap_columns(x):
  """
  Custom logic to try and remap columns from a BED file
  """
  if x == 0:
      return 'chr'
  elif x == 1:
      return 'start'
  elif x == 2:
      return 'end'
  else:
      return 'ucol'+str(x)


@data_factory('BED data loader', is_bed, priority=999)
def read_bed(file_name):
    """
    Read a bed file into glue.
    """
    bed_data = pd.read_csv(file_name, sep='\t', header=None)
    bed_data.rename(columns=remap_columns,inplace=True) #Rename first three columns as standard
    
    return Data(**{k:bed_data[k] for k in bed_data.columns},
                        label=Path(file_name).stem)