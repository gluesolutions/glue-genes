import dask.array as da
import numpy as np
import pytest
from glue.core import Data
from glue.core.fixed_resolution_buffer import compute_fixed_resolution_buffer
from glue.core.roi import RangeROI, RectangularROI
from glue.core.subset import RangeSubsetState, RoiSubsetState
from numpy.testing import assert_equal

from glue_genes.glue_genomics_data.multires_data import MultiResolutionData


def downsample_xy(a, zN=4):
    zN = zN
    BSZ = (2, 2)
    m, n = a.shape[1:]
    return a.reshape(zN, m // BSZ[0], BSZ[0], n // BSZ[1], BSZ[1]).mean(axis=(2, 4))


a = np.arange(4096).reshape((4, 32, 32))


def test_downsample():
    """
    https://stackoverflow.com/questions/44527579/whats-the-best-way-to-downsample-a-numpy-array
    """
    a = np.arange(4096).reshape((4, 32, 32))
    # We make small chunks to simulate the typical case
    # where we have very large dask arrays with lots of chunks
    a_d = da.from_array(a, chunks=2)
    b = downsample_xy(a, zN=4)
    b_d = da.from_array(b, chunks=2)
    assert b.shape == (4, 16, 16)

    assert np.max(b_d).compute() < np.max(a_d).compute()
    assert np.min(b_d).compute() > np.min(a_d).compute()

    assert np.max(b) < np.max(a)
    assert np.min(b) > np.min(a)

    down_data = MultiResolutionData(
        x=a,
        all_resolutions=[{"x": a}, {"x": b}],
        label="Downsampled data",
        reduced_dims=[1, 2],
    )
    est_min = down_data.compute_statistic(
        "minimum", down_data.main_components[0], random_subset=1000
    )
    est_max = down_data.compute_statistic(
        "maximum", down_data.main_components[0], random_subset=1000
    )
    # assert est_min == 2
    assert np.isclose(np.min(b), est_min, atol=10)
    assert np.isclose(np.max(b), est_max, atol=10)

    # assert np.isclose(est_min, 0.5) == 0
    # assert pytest.approx(est_max, 0.5) == 4096

    down_data_dask = MultiResolutionData(
        x=a_d,
        all_resolutions=[{"x": a_d}, {"x": b_d}],
        label="Downsampled data",
        reduced_dims=[1, 2],
    )
    est_min = down_data_dask.compute_statistic(
        "minimum", down_data_dask.main_components[0], random_subset=1000
    )
    est_max = down_data_dask.compute_statistic(
        "maximum", down_data_dask.main_components[0], random_subset=1000
    )
    # assert est_min == 2
    assert np.isclose(np.min(b_d).compute(), est_min)
    assert np.isclose(np.max(b_d).compute(), est_max)


x1 = np.eye(8, 8)
x2 = np.eye(4, 4)
x3 = np.eye(2, 2)
base_data = MultiResolutionData(
    x=x1,
    all_resolutions=[{"x": x1}, {"x": x2}, {"x": x3}],
    label="Test data",
    reduced_dims=[0, 1],
)

x1 = da.asarray(np.eye(8, 8))
x2 = da.asarray(np.eye(4, 4))
x3 = da.asarray(np.eye(2, 2))
dask_data = MultiResolutionData(
    x=x1,
    all_resolutions=[{"x": x1}, {"x": x2}, {"x": x3}],
    label="Test data",
    reduced_dims=[0, 1],
)


class TestMultiResolutionData:
    def setup_class(self):
        self.data = base_data
        self.dask_data = dask_data

    def test_basic_init(self):
        d2 = self.data._reduced_res_data_sets[0]
        assert len(self.data._reduced_res_data_sets) == 2
        assert d2.scale_factor == 2
        assert_equal(
            d2[d2._downsampled_pixel_cids[0]].data,
            np.array([[0, 0, 0, 0], [2, 2, 2, 2], [4, 4, 4, 4], [6, 6, 6, 6]]),
        )
        assert d2[d2._downsampled_pixel_cids[0]].data.shape == d2.shape

    def test_cid_lookup(self):
        red_data = self.data._reduced_res_data_sets[0]
        for red_comp, parent_comp in zip(
            red_data.main_components, self.data.main_components
        ):
            assert red_data.convert_reduced_to_full_cid(red_comp) == parent_comp

    def test_get_mask(self):
        subset = self.data.new_subset()
        roi = RectangularROI(xmin=-1, xmax=3, ymin=-1, ymax=5)
        subset.subset_state = RoiSubsetState(
            xatt=self.data.pixel_component_ids[1],
            yatt=self.data.pixel_component_ids[0],
            roi=roi,
        )
        d3 = self.data._reduced_res_data_sets[1]
        red_mask = d3.get_mask(subset.subset_state)
        assert_equal(red_mask, np.array([[True, False], [True, False]]))
        roi = RectangularROI(xmin=-1, xmax=3, ymin=-1, ymax=9)
        subset.roi = roi
        assert_equal(red_mask, np.array([[True, False], [True, False]]))

        range_roi = RangeROI("x", -1, 3)
        subset.roi = range_roi
        assert_equal(red_mask, np.array([[True, False], [True, False]]))

    @pytest.mark.parametrize("data", [(base_data), (dask_data)])
    def test_fixed_resolution_buffer(self, data):
        assert_equal(
            data.compute_fixed_resolution_buffer(
                [(0, 7, 8), (0, 7, 8)],
                target_data=data,
                target_cid=data.main_components[0],
            ),
            np.eye(8, 8),
        )
        assert_equal(
            data.compute_fixed_resolution_buffer(
                [(0, 7, 4), (0, 7, 4)],
                target_data=data,
                target_cid=data.main_components[0],
            ),
            np.eye(4, 4),
        )

    @pytest.mark.parametrize("data", [(base_data), (dask_data)])
    def test_fixed_resolution_buffer_with_subset(self, data):
        subset = data.new_subset()
        subset.subset_state = RangeSubsetState(-1, 0.5, data.main_components[0])
        assert_equal(
            data.compute_fixed_resolution_buffer(
                [(0, 7, 8), (0, 7, 8)],
                target_data=data,
                subset_state=subset.subset_state,
            ),
            ~np.eye(8, 8).astype(bool),
        )


def test_frb_with_multi_resolution_data():
    x = da.arange(512).reshape((8, 8, 8))
    y = da.arange(64).reshape((4, 4, 4))

    data1 = MultiResolutionData(x=x, all_resolutions=[{"x": x}, {"x": y}], label="d1")

    assert_equal(
        data1.compute_fixed_resolution_buffer(
            [1, (-1, 1, 3), (0, 3, 4)],
            target_data=data1,
            target_cid=data1.main_components[0],
        ),
        np.array(
            [
                [-np.inf, -np.inf, -np.inf, -np.inf],
                [64, 65, 66, 67],
                [72, 73, 74, 75],
            ]
        ),
    )

    assert_equal(
        data1.compute_fixed_resolution_buffer(
            [1, (-10, 10, 1), (0, 20, 1)],
            target_data=data1,
            target_cid=data1.main_components[0],
        ),
        np.array([[-np.inf]]),
    )


def test_frb_with_sliced_multi_resolution_data():
    x = da.arange(512).reshape((1, 1, 8, 8, 8))
    y = da.arange(64).reshape((1, 1, 4, 4, 4))

    x = x[0, 0, ...]
    y = y[0, 0, ...]

    data1 = MultiResolutionData(x=x, all_resolutions=[{"x": x}, {"x": y}], label="d1")

    assert_equal(
        data1.compute_fixed_resolution_buffer(
            [1, (-1, 1, 3), (0, 3, 4)],
            target_data=data1,
            target_cid=data1.main_components[0],
        ),
        np.array(
            [
                [-np.inf, -np.inf, -np.inf, -np.inf],
                [64, 65, 66, 67],
                [72, 73, 74, 75],
            ]
        ),
    )

    assert_equal(
        data1.compute_fixed_resolution_buffer(
            [1, (-10, 10, 1), (0, 20, 1)],
            target_data=data1,
            target_cid=data1.main_components[0],
        ),
        np.array([[-np.inf]]),
    )


def test_frb_with_dask():
    data1 = Data(x=da.arange(12).reshape((1, 4, 3)), label="d1")

    assert_equal(
        compute_fixed_resolution_buffer(
            data1,
            target_data=data1,
            bounds=[(-1, 1, 3), (0, 3, 4), 1],
            target_cid=data1.main_components[0],
        ),
        np.array(
            [
                [-np.inf, -np.inf, -np.inf, -np.inf],
                [1, 4, 7, 10],
                [-np.inf, -np.inf, -np.inf, -np.inf],
            ]
        ),
    )


def test_frb_with_sliced_dask():
    x = da.arange(12).reshape((1, 1, 4, 3))
    y = x[:, 0, ...]
    data1 = Data(y=y, label="d1")

    assert_equal(
        compute_fixed_resolution_buffer(
            data1,
            target_data=data1,
            bounds=[(-1, 1, 3), (0, 3, 4), 1],
            target_cid=data1.main_components[0],
        ),
        np.array(
            [
                [-np.inf, -np.inf, -np.inf, -np.inf],
                [1, 4, 7, 10],
                [-np.inf, -np.inf, -np.inf, -np.inf],
            ]
        ),
    )
