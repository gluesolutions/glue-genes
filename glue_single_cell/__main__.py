import anndata
from glue.app.qt.application import GlueApplication
from glue.core import DataCollection

from .data import DataAnnData


def main(argv):

    adata = anndata.read_h5ad("def_sparse.h5ad", backed="r")
    d = DataAnnData(adata, label="test_name")
    dc = DataCollection([d])
    ga = GlueApplication(dc)
    ga.start()


if __name__ == "__main__":
    import sys

    main(sys.argv)
