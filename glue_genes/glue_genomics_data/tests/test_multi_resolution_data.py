import dask.array as da
import numpy as np
from glue.core import Data
from glue.core.fixed_resolution_buffer import compute_fixed_resolution_buffer
from numpy.testing import assert_equal

from glue_genes.glue_genomics_data.multires_data import MultiResolutionData


class TestMultiResolutionData:
    def setup_class(self):
        x1 = np.eye(64, 64)
        x2 = np.eye(32, 32)
        x3 = np.eye(16, 16)
        self.data = MultiResolutionData(
            x=x1, all_resolutions=[{"x": x1}, {"x": x2}, {"x": x3}], label="Test data"
        )

    def test_basic_init(self):
        d2 = self.data._reduced_res_data_sets[0]
        assert len(self.data._reduced_res_data_sets) == 2
        # assert d2._cid_to_parent_cid == {}

    def test_cid_lookup(self):
        red_data = self.data._reduced_res_data_sets[0]
        for red_comp, parent_comp in zip(
            red_data.main_components, self.data.main_components
        ):
            assert red_data._translate_reduced_cid_to_full_cid(red_comp) == parent_comp


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
