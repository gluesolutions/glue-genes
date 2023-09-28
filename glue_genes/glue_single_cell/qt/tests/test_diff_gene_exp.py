import os
import pytest

from glue_qt.app import GlueApplication
from glue.core import data_factories as df

from glue_genes.glue_single_cell.anndata_factory import read_anndata

from ..diff_gene_exp import get_gene_diff_exp

DATA = os.path.join(os.path.dirname(__file__), "data")


class TestDiffGeneExp(object):
    def get_data(self, try_backed=False):
        data = df.load_data(
            os.path.join(DATA, "test_data.h5ad"),
            factory=read_anndata,
            skip_dialog=True,
            try_backed=try_backed,
        )
        return data

    def setup_method(self):
        self.app = GlueApplication()
        self.session = self.app.session
        self.hub = self.session.hub

        self.dc = self.session.data_collection
        self.dc.append(self.get_data())
        assert len(self.dc) == 3
        self.dd = self.dc[0]
        self.var_data = self.dc[1]
        self.obs_data = self.dc[2]

        self.subset1 = self.dc.new_subset_group(
            label="B cells", subset_state=self.obs_data.id["cell_type"] == "B"
        )
        self.subset2 = self.dc.new_subset_group(
            label="T cells", subset_state=self.obs_data.id["cell_type"] == "T"
        )
        self.subset_single = self.dc.new_subset_group(
            label="SinglePoint", subset_state=self.obs_data.id["Pixel Axis 0 [x]"] == 0
        )
        assert len(self.dd.meta['obs_data'].subsets[0]['cell_type']) == 15
        assert len(self.dd.meta['obs_data'].subsets[1]['cell_type']) == 48
        assert len(self.dd.meta['obs_data'].subsets[2]['cell_type']) == 1

    def test_diff_gene_exp(self):
        df = get_gene_diff_exp(
            self.subset1, self.subset2, self.dd
        )
        assert len(df) == 2000

        assert df["pvals_adj"].min() < 0.5
        assert df["pvals_adj"].min() > 0.4

        assert df["pvals_adj"].max() > 0.95

    def test_diff_gene_exp_versus_rest(self):
        df = get_gene_diff_exp(self.subset2, None, self.dd)
        assert len(df) == 2000
        assert df["pvals_adj"].min() < 0.25
        assert df["pvals_adj"].min() > 0.05

        assert df["pvals_adj"].max() > 0.95

    def test_diff_gene_exp_fails(self):
        with pytest.raises(ValueError):
            _ = get_gene_diff_exp(self.subset_single, None, self.dd)


class TestDiffGeneExpBacked(TestDiffGeneExp):
    def get_data(self, try_backed=True):
        data = df.load_data(
            os.path.join(DATA, "test_data.h5ad"),
            factory=read_anndata,
            skip_dialog=True,
            try_backed=try_backed,
        )
        return data
