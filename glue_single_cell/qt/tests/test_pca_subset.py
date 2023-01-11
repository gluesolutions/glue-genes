"""
To test this properly we need a way to read in some test files
since otherwise we have to duplicate or abstract out all the
linking codes

Which implies we need some (small) test AnnData files

And a way to deal with skipping the dialog in the anndata_factory

Known problem:

You can select a subset that is not defined on the same dataset
and that produces an IncompatibleAttribute exception


"""
import os
from unittest.mock import patch

import numpy as np
from glue.app.qt import GlueApplication
from glue.core import Data
from glue.core import data_factories as df
from glue.core.data_collection import DataCollection
from glue.core.link_helpers import JoinLink
from glue.core.state import GlueUnSerializer

from glue_single_cell.anndata_factory import read_anndata

from ..pca_subset import (PCASubsetDialog, apply_data_arr,
                          do_calculation_over_gene_subset)

DATA = os.path.join(os.path.dirname(__file__), "data")


class TestCellSummarySession(object):
    def test_session_save_and_restore(self, tmpdir):
        self.app = GlueApplication()
        self.session = self.app.session
        self.hub = self.session.hub

        self.dc = self.session.data_collection
        d1 = df.load_data(
            os.path.join(DATA, "test_data.h5ad"), factory=read_anndata, skip_dialog=True
        )
        d2 = df.load_data(
            os.path.join(DATA, "test_other_data.h5ad"),
            factory=read_anndata,
            skip_dialog=True,
        )

        d3 = Data(
            gene_id=["Gene_2", "Gene_3", "Gene_5", "Gene_8", "Gene_11", "Gene_12"],
            qtl=[1, 2, 3, 4, 5, 5],
            label="qtl",
        )
        self.dc.append(d1)
        self.dc.append(d2)
        self.dc.append(d3)

        assert len(self.dc[2].components) == 5

        d1_var = self.dc[1]

        s = self.dc.new_subset_group()

        # s = d1_var.new_subset()
        s.subset_state = d1_var.id["gene_stuff_0"] > 0

        sumdiag = PCASubsetDialog(self.dc)
        sumdiag.state.data = self.dc[2]
        sumdiag.state.genesubset = s
        sumdiag.state.do_means = True
        sumdiag.state.do_pca = False
        sumdiag.state.do_module = False

        with patch("glue_single_cell.qt.pca_subset.dialog") as fakedialog:  # noqa: F841
            sumdiag._apply()

        assert len(self.dc[0].listeners) == 1
        assert len(self.dc[2].components) == 6
        assert np.sum(self.dc[2]["Subset 1_Means_0"]) > 99
        assert np.sum(self.dc[2]["Subset 1_Means_0"]) < 100

        sumdiag.state.genesubset.subset_state = d1_var.id["gene_stuff_0"] > 0.6

        assert np.sum(self.dc[2]["Subset 1_Means_0"]) > 99
        assert (
            not np.sum(self.dc[2]["Subset 1_Means_0"]) < 100
        )  # Check that updating subset changes the result

        filename = tmpdir.join("test_anndata_load_session.glu").strpath
        self.session.application.save_session(filename)

        with open(filename, "r") as f:
            session = f.read()

        state = GlueUnSerializer.loads(session)

        ga = state.object("__main__")

        dc = ga.session.data_collection

        assert len(dc) == 7

        assert len(dc[0].listeners) == 1
        assert len(dc[2].components) == 6

        dc.subset_groups[0].subset_state = dc[1].id["gene_stuff_0"] > 0

        assert np.sum(dc[2]["Subset 1_Means_0"]) > 99
        assert np.sum(dc[2]["Subset 1_Means_0"]) < 100
        ga.close()


class TestCellSummary(object):
    def setup_method(self, method):

        d1 = df.load_data(
            os.path.join(DATA, "test_data.h5ad"), factory=read_anndata, skip_dialog=True
        )
        d2 = df.load_data(
            os.path.join(DATA, "test_other_data.h5ad"),
            factory=read_anndata,
            skip_dialog=True,
        )

        d3 = Data(
            gene_id=["Gene_2", "Gene_3", "Gene_5", "Gene_8", "Gene_11", "Gene_12"],
            qtl=[1, 2, 3, 4, 5, 5],
            label="qtl",
        )
        self.dc = DataCollection([d1, d2, d3])

        self.dc.append(d1)
        self.dc.append(d2)
        # TODO: Implicit indexing of the returned datasets here to get the var
        #       array is not ideal
        mylink = JoinLink(
            cids1=[d1[1].id["var_names"]],
            cids2=[d2[1].id["var_names"]],
            data1=d1[1],
            data2=d2[1],
        )
        self.dc.add_link(mylink)

        qtllink = JoinLink(
            cids1=[d1[1].id["var_names"]],
            cids2=[d3.id["gene_id"]],
            data1=d1[1],
            data2=d3,
        )
        self.dc.add_link(qtllink)

    def test_data_load_and_links(self):
        assert len(self.dc) == 7  # Each datafile produces 3 datasets
        d1_var = self.dc[1]
        d2_var = self.dc[4]

        assert d1_var.components == [
            "Pixel Axis 0 [x]",
            "var_names",
            "gene_stuff_0",
            "gene_stuff_1",
            "gene_stuff_2",
            "gene_stuff_3",
            "gene_stuff_4",
        ]

        s = d1_var.new_subset()
        s.subset_state = d1_var.id["gene_stuff_0"] > 0
        s = d2_var.new_subset()
        s.subset_state = d1_var.id["gene_stuff_0"] > 0
        assert s.to_mask().sum() > 0

    def do_calculation_on_multiple_datasets(self, **kwargs):
        d1_adata = self.dc[0]
        d1_var = self.dc[1]

        s = d1_var.new_subset()
        s.subset_state = d1_var.id["gene_stuff_0"] > 0

        data_arr = do_calculation_over_gene_subset(d1_adata, s, **kwargs)
        assert len(data_arr) == 100

        # Now use the same subset on the other dataset
        d2_adata = self.dc[3]
        d2_var = self.dc[4]

        s = d2_var.new_subset()
        s.subset_state = d1_var.id["gene_stuff_0"] > 0

        data_arr = do_calculation_over_gene_subset(d2_adata, s, **kwargs)
        assert len(data_arr) == 500

    def test_calculation(self, **kwargs):
        self.do_calculation_on_multiple_datasets(calculation="Means")
        # If Module and backed mode we get this error
        # https://github.com/scverse/scanpy/issues/2153
        self.do_calculation_on_multiple_datasets(calculation="Module")

    def test_calculation_through_qtl(self):
        """
        Test that defining a subset on a linked third dataset works
        """
        d1_adata = self.dc[0]
        d1_var = self.dc[1]
        qtl = self.dc[6]

        s = d1_var.new_subset()
        s.subset_state = qtl.id["qtl"] > 2

        data_arr = do_calculation_over_gene_subset(d1_adata, s, calculation="Means")
        assert len(data_arr) == 100

        # Now use the same subset on the other dataset
        d2_adata = self.dc[3]
        d2_var = self.dc[4]

        s = d2_var.new_subset()
        s.subset_state = qtl.id["qtl"] > 2

        data_arr = do_calculation_over_gene_subset(d2_adata, s, calculation="Means")
        assert len(data_arr) == 500

    def test_calculation_through_qtl_subsets(self):
        """
        Test that defining a subset on a linked third dataset works
        """
        d1_adata = self.dc[0]
        qtl = self.dc[6]

        state = qtl.id["qtl"] > 2
        subset_group = self.dc.new_subset_group("QTL genes", state)
        # This is what we do in the dialog
        for subset in subset_group.subsets:
            if subset.data == self.dc[0].meta["var_data"]:
                genesubset = subset

        data_arr = do_calculation_over_gene_subset(
            d1_adata, genesubset, calculation="Means"
        )
        assert len(data_arr) == 100

        # Now use the same subset on the other dataset
        d2_adata = self.dc[3]

        for subset in subset_group.subsets:
            if subset.data == self.dc[3].meta["var_data"]:
                genesubset = subset

        data_arr = do_calculation_over_gene_subset(
            d2_adata, genesubset, calculation="Means"
        )
        assert len(data_arr) == 500

    def test_wrong_subset(self):
        """
        Test that defining a subset on cells does not crash things
        (and also that subset definitions on obs to not propagate
        to var through the X array).
        """
        d1_adata = self.dc[0]
        d2_obs = self.dc[5]
        state = d2_obs.id["cell_type"] == "T"

        subset_group = self.dc.new_subset_group("T cells", state)

        # This is what we do in the dialog
        for subset in subset_group.subsets:
            if subset.data == self.dc[0].meta["var_data"]:
                genesubset = subset
        with patch("glue_single_cell.qt.pca_subset.dialog") as fakedialog:  # noqa: F841
            data_arr = do_calculation_over_gene_subset(
                d1_adata, genesubset, calculation="Means"
            )
        assert data_arr is None

    def test_adding_to_dataset(self):
        d1_adata = self.dc[0]
        d1_var = self.dc[1]
        d1_obs = self.dc[2]  # This is target_dataset

        assert len(d1_obs.components) == 5

        s = d1_var.new_subset()
        s.subset_state = d1_var.id["gene_stuff_0"] > 0

        data_arr = do_calculation_over_gene_subset(d1_adata, s, calculation="Means")
        apply_data_arr(d1_obs, data_arr, "d1", key="Means")
        assert len(d1_obs.components) == 6


class TestCellSummaryBacked(TestCellSummary):
    def setup_method(self, method):

        d1 = df.load_data(
            os.path.join(DATA, "test_data.h5ad"),
            factory=read_anndata,
            skip_dialog=True,
            try_backed=True,
        )
        d2 = df.load_data(
            os.path.join(DATA, "test_other_data.h5ad"),
            factory=read_anndata,
            skip_dialog=True,
            try_backed=True,
        )

        d3 = Data(
            gene_id=["Gene_2", "Gene_3", "Gene_5", "Gene_8", "Gene_11", "Gene_12"],
            qtl=[1, 2, 3, 4, 5, 5],
            label="qtl",
        )
        # self.app = GlueApplication()
        # self.session = self.app.session
        # self.hub = self.session.hub
        self.dc = DataCollection([d1, d2, d3])

        self.dc.append(d1)
        self.dc.append(d2)
        # TODO: Implicit indexing of the returned datasets here to get the var
        #       array is not ideal
        mylink = JoinLink(
            cids1=[d1[1].id["var_names"]],
            cids2=[d2[1].id["var_names"]],
            data1=d1[1],
            data2=d2[1],
        )
        self.dc.add_link(mylink)

        qtllink = JoinLink(
            cids1=[d1[1].id["var_names"]],
            cids2=[d3.id["gene_id"]],
            data1=d1[1],
            data2=d3,
        )
        self.dc.add_link(qtllink)
