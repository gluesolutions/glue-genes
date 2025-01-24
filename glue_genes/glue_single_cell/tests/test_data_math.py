import os

import anndata
import context  # noqa F401
import numpy as np
import pytest
from glue.core import data_factories as df
from numpy.random import RandomState
from numpy.testing import assert_almost_equal, assert_equal
from scipy.sparse import csc_matrix, find

from glue_genes.glue_single_cell.anndata_factory import read_anndata

SPARSE_BACKED_OBS_NUM = 500
SPARSE_BACKED_VAR_NUM = 700


@pytest.fixture
def data_sparse_backed(tmpdir):
    rs = RandomState(12345)
    C = rs.rand(SPARSE_BACKED_OBS_NUM, SPARSE_BACKED_VAR_NUM)
    C[C < 0.90] = 0
    C = csc_matrix(C)
    new_adata = anndata.AnnData(X=C, dtype="float64")
    file_location = os.path.join(tmpdir, "test_dataset_float64.h5ad")
    new_adata.write(filename=file_location)
    d = df.load_data(
        file_location, factory=read_anndata, skip_dialog=True, try_backed=True
    )
    yield C, d


@pytest.fixture
def data_sparse_inmemory(tmpdir):
    rs = RandomState(12345)
    C = rs.rand(SPARSE_BACKED_OBS_NUM, SPARSE_BACKED_VAR_NUM)
    C[C < 0.90] = 0
    C = csc_matrix(C)
    new_adata = anndata.AnnData(X=C, dtype="float64")
    file_location = os.path.join(tmpdir, "test_dataset_float64.h5ad")
    new_adata.write(filename=file_location)
    d = df.load_data(
        file_location, factory=read_anndata, skip_dialog=True, try_backed=False
    )
    yield C, d


def test_data_setup_sparse_backed(data_sparse_backed):
    C, d = data_sparse_backed
    assert len(d[0]._components) == 5
    assert d[0].sparse is True
    assert d[0].backed is True
    assert isinstance(
        d[0].Xdata, anndata.AnnData
    )  # As a reference to the full AnnData object


def test_data_setup_sparse_inmemory(data_sparse_inmemory):
    C, d = data_sparse_inmemory
    assert len(d[0]._components) == 5
    assert d[0].sparse is True
    assert d[0].backed is False
    assert isinstance(d[0].Xdata, anndata.AnnData)


def test_get_data_view_sparse_backed(data_sparse_backed):
    C, d = data_sparse_backed
    views = (
        np.s_[:],
        np.s_[0, :],
        np.s_[:, 0],
        np.s_[:, 10:30],
        np.s_[37:87, :],
    )
    for view in views:
        for comp in d[0].components:
            if comp in d[0]._pixel_component_ids:
                pass
            else:
                comp_data = d[0].get_data(comp, view)
                assert comp_data.size == C[view].size
                assert_equal(comp_data, C[view].data)


def test_get_data_view_sparse_inmemory(data_sparse_inmemory):
    C, d = data_sparse_inmemory
    views = (
        np.s_[:],
        np.s_[0, :],
        np.s_[:, 0],
        np.s_[:, 10:30],
        np.s_[37:87, :],
    )
    for view in views:
        for comp in d[0].main_components:
            comp_data = d[0].get_data(comp, view, keep_sparse=True)
            assert comp_data.size == C[view].size
            assert_equal(comp_data.toarray(), C[view].toarray())


def test_get_data_sparse_backed(data_sparse_backed):
    C, d = data_sparse_backed
    for comp in d[0].main_components:
        comp_data = d[0].get_data(comp)
        assert len(comp_data) == C.size


def test_get_data_sparse_inmemory(data_sparse_inmemory):
    C, d = data_sparse_inmemory
    for comp in d[0].main_components:
        comp_data = d[0].get_data(comp, keep_sparse=True)
        assert comp_data.toarray().shape == C.toarray().shape


def test_make_histogram_sparse_backed(data_sparse_backed):
    C, d = data_sparse_backed
    xmax, ymax = C.shape
    hist_range = ((0, xmax), (0, ymax))
    bins = (
        11,
        11,
    )  # Histogram logic only works if range//bins = integer!! Maybe we don't care too much.
    x, y, w = find(C)
    correct_histogram = np.histogram2d(
        x=x, y=y, bins=bins, range=hist_range, weights=w
    )[0]
    histogram = d[0].compute_histogram(
        cids=["Pixel Axis 0 [y]", "Pixel Axis 1 [x]"],
        bins=bins,
        range=hist_range,
        weights=["X"],
    )
    assert_almost_equal(correct_histogram, histogram)
    assert np.allclose(correct_histogram, histogram)


def test_make_histogram_sparse_inmemory(data_sparse_inmemory):
    C, d = data_sparse_inmemory
    xmax, ymax = C.shape
    hist_range = ((0, xmax), (0, ymax))
    bins = (
        11,
        11,
    )  # Histogram logic only works if range//bins = integer!! Maybe we don't care too much.
    x, y, w = find(C)
    correct_histogram = np.histogram2d(
        x=x, y=y, bins=bins, range=hist_range, weights=w
    )[0]
    histogram = d[0].compute_histogram(
        cids=["Pixel Axis 0 [y]", "Pixel Axis 1 [x]"],
        bins=bins,
        range=hist_range,
        weights=["X"],
    )
    assert_almost_equal(
        correct_histogram, histogram, decimal=5
    )  # Lower tolerance required here
    assert np.allclose(correct_histogram, histogram)
