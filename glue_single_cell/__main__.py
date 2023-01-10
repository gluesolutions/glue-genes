import numpy as np
import pandas as pd
from scipy import sparse

from glue.core.component_id import ComponentID
from glue.core.data import BaseCartesianData
from glue.utils import view_shape
from .data import DataAnnData

# We now create a data object using the above class,
# and launch a a glue session

from glue.core import DataCollection
from glue.app.qt.application import GlueApplication

import anndata

def main(argv):
    #def load_data():
    #df = pd.read_csv('sparse_log_counts_sub_plus_labels2.csv')
    #sparse_coo = sparse.coo_matrix((df['x'].values,
    #                               (pd.factorize(df['cell_id'])[0],
    #                                pd.factorize(df['gene'])[0])))
    #sparse_csc = sparse_coo.tocsc()
    #adata = anndata.AnnData(X=sparse_csc)
    #adata.filename = 'test_backfile_small.h5ad'
    
    adata = anndata.read_h5ad('def_sparse.h5ad',backed='r')
    d = DataAnnData(adata,label="test_name")
    dc = DataCollection([d])
    ga = GlueApplication(dc)
    ga.start()


if __name__ == '__main__':
    import sys
    main(sys.argv)