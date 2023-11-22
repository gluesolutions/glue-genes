import os
import pytest

from glue_qt.app import GlueApplication
from glue.core import data_factories as df

from glue_genes.glue_single_cell.anndata_factory import read_anndata

from ..diff_gene_exp import get_gene_diff_exp, DiffGeneExpDialog

DATA = os.path.join(os.path.dirname(__file__), "data")


# If SyncComponents kept their sync state, we could use this to check
def has_insync_component(data):
    for comp in data.main_components:
        if 'sync' in comp.label:
            return True
    return False


def has_outofsync_component(data):
    for comp in data.main_components:
        if 'out-of-sync' in comp.label:
            return True
    return False


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
        # The first ten genes in the table
        self.gene_subset = self.dc.new_subset_group(
            label="Genes", subset_state=self.var_data.id["Pixel Axis 0 [x]"] < 10
        )

        # A subset of genes on the X array will be over all observations, but is
        # not a valid subset of observations for differential gene expression
        self.gene_subset_from_X_array = self.dc.new_subset_group(
            label="Genes from X array", subset_state=self.dd.id["Pixel Axis 1 [x]"] < 10
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

    def test_dge_fails_too_small(self):
        with pytest.raises(ValueError, match='less than 2 observations'):
            _ = get_gene_diff_exp(self.subset_single, None, self.dd)

    def test_dge_fails_bad_subset(self):
        with pytest.raises(ValueError,  match='Failed to generate a mask over observations'):
            _ = get_gene_diff_exp(self.gene_subset, None, self.dd)

        with pytest.raises(ValueError, match='includes all the observations'):
            _ = get_gene_diff_exp(self.gene_subset_from_X_array, None, self.dd)

    def test_dge_fails_bad_subsets_two_subsets(self):
        with pytest.raises(ValueError, match='Failed to generate a mask over observations'):
            _ = get_gene_diff_exp(self.subset1, self.gene_subset, self.dd)

        with pytest.raises(ValueError, match='less than 2 observations'):
            _ = get_gene_diff_exp(self.subset_single, self.subset2, self.dd)

    def test_full_with_dialog(self):

        d = DiffGeneExpDialog(self.dc)
        d.state.subset1 = self.subset1
        d.state.subset2 = self.subset2
        d._apply(show_dialog=False)
        assert len(d.state.data.listeners) == 1
        assert has_insync_component(self.var_data)

        d.state.subset1.subset_state = self.obs_data.id["Pixel Axis 0 [x]"] == 0
        assert has_outofsync_component(self.var_data)

        d.state.subset1.subset_state = self.obs_data.id["cell_type"] == "B"

        assert not has_outofsync_component(self.var_data)
        assert has_insync_component(self.var_data)

        for subset in d.state.subset1.subsets:
            subset.delete()
        assert len(d.state.data.listeners) == 0
        assert not has_insync_component(self.var_data)


class TestDiffGeneExpBacked(TestDiffGeneExp):
    def get_data(self, try_backed=True):
        data = df.load_data(
            os.path.join(DATA, "test_data.h5ad"),
            factory=read_anndata,
            skip_dialog=True,
            try_backed=try_backed,
        )
        return data
