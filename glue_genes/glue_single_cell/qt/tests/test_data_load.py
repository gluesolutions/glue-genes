import os

import anndata as ad
from glue_qt.app import GlueApplication
from glue.core import data_factories as df
from glue.core.state import GlueUnSerializer
from numpy.testing import assert_array_equal

from glue_genes.glue_single_cell.anndata_factory import read_anndata

DATA = os.path.join(os.path.dirname(__file__), "data")


class TestData(object):
    def get_data(self, try_backed=False):
        original_anndata = ad.read_h5ad(os.path.join(DATA, "test_data.h5ad"))
        # Calling self.app.load_data() here does NOT skip_dialog
        data = df.load_data(
            os.path.join(DATA, "test_data.h5ad"),
            factory=read_anndata,
            skip_dialog=True,
            try_backed=try_backed,
        )
        return (original_anndata, data)

    def test_anndata_access(self):
        self.app = GlueApplication()
        self.session = self.app.session
        self.hub = self.session.hub

        self.data_collection = self.session.data_collection
        original_anndata, data = self.get_data()
        self.data_collection.append(data)

        assert len(self.data_collection) == 3
        dd = self.data_collection[0]

        # The new anndata object is not exactly the same, since we get the obsm and varm arrays mixed
        assert_array_equal(original_anndata.obs["cell_type"], dd.Xdata.obs["cell_type"])
        # assert original_anndata.var['gene_stuff'] == dd.Xdata.var['gene_stuff']

    def test_save_restore(self, tmpdir):

        self.app = GlueApplication()
        self.session = self.app.session
        self.hub = self.session.hub
        self.data_collection = self.session.data_collection

        orig, data = self.get_data()
        self.data_collection.append(data)
        assert len(self.data_collection) == 3

        filename = tmpdir.join("test_anndata_save_restore_session.glu").strpath
        self.session.application.save_session(filename)
        with open(filename, "r") as f:
            session = f.read()

        state = GlueUnSerializer.loads(session)
        ga = state.object("__main__")
        dc = ga.session.data_collection
        assert dc[0].Xdata is not None


class TestDataBacked(TestData):
    def get_data(self, try_backed=True):
        original_anndata = ad.read_h5ad(os.path.join(DATA, "test_data.h5ad"))
        # Calling self.app.load_data() here does NOT skip_dialog
        data = df.load_data(
            os.path.join(DATA, "test_data.h5ad"),
            factory=read_anndata,
            skip_dialog=True,
            try_backed=try_backed,
        )
        return (original_anndata, data)


class TestAnnDataLoader(object):
    def test_session_save_and_restore(self, tmpdir):
        self.app = GlueApplication()
        self.session = self.app.session
        self.hub = self.session.hub

        self.data_collection = self.session.data_collection

        # Calling self.app.load_data() here does NOT skip_dialog
        data = df.load_data(
            os.path.join(DATA, "test_data.h5ad"), factory=read_anndata, skip_dialog=True
        )
        self.data_collection.append(data)

        assert len(self.data_collection) == 3

        filename = tmpdir.join("test_anndata_load_session.glu").strpath
        self.session.application.save_session(filename)
        with open(filename, "r") as f:
            session = f.read()

        state = GlueUnSerializer.loads(session)
        ga = state.object("__main__")
        dc = ga.session.data_collection

        assert len(dc) == 3
        ga.close()

    def test_session_save_and_restore_with_options(self, tmpdir):
        self.app = GlueApplication()
        self.session = self.app.session
        self.hub = self.session.hub
        self.data_collection = self.session.data_collection

        # Calling self.app.load_data() here does NOT skip_dialog
        data = df.load_data(
            os.path.join(DATA, "test_data.h5ad"),
            factory=read_anndata,
            skip_dialog=True,
            try_backed=False,
            skip_components=["gene_stuff"],
        )
        self.data_collection.append(data)

        assert len(self.data_collection) == 3
        filename = tmpdir.join("test_anndata_kwargs.glu").strpath

        self.session.application.save_session(filename)
        with open(filename, "r") as f:
            session = f.read()

        state = GlueUnSerializer.loads(session)

        ga = state.object("__main__")

        dc = ga.session.data_collection

        assert len(dc) == 3
        assert len(dc[1].components) == 2
        ga.close()

    def test_session_save_and_restore_with_two_datasets(self, tmpdir):
        self.app = GlueApplication()
        self.session = self.app.session
        self.hub = self.session.hub
        self.data_collection = self.session.data_collection

        # Calling self.app.load_data() here does NOT skip_dialog
        data = df.load_data(
            os.path.join(DATA, "test_data.h5ad"),
            factory=read_anndata,
            skip_dialog=True,
            try_backed=False,
            skip_components=["gene_stuff"],
        )
        self.data_collection.append(data)
        data2 = df.load_data(
            os.path.join(DATA, "test_other_data.h5ad"),
            factory=read_anndata,
            skip_dialog=True,
            try_backed=False,
        )
        self.data_collection.append(data2)

        assert len(self.data_collection) == 6

        filename = tmpdir.join("test_anndata_twodata.glu").strpath

        self.session.application.save_session(filename)
        with open(filename, "r") as f:
            session = f.read()

        state = GlueUnSerializer.loads(session)

        ga = state.object("__main__")

        dc = ga.session.data_collection

        assert len(dc) == 6
        assert len(dc[1].components) == 2

        assert len(dc[4].components) > 3

        ga.close()
