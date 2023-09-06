"""
Generate some test OME-Zarr data and test that it can be loaded into glue
using the ome-zarr data factory.
"""
import pytest
import numpy as np
import zarr
from ome_zarr.io import parse_url
from ome_zarr.writer import write_image
from glue.core.data_factories import load_data
from glue_genes.glue_genomics_data.ome_zarr_factory import read_ome_zarr

MEAN_VAL = 10
SIZE_XY = 128


@pytest.fixture(scope="session")
def test_zyx(tmp_path_factory):
    path = tmp_path_factory.mktemp("test_zyx.zarr", numbered=False)

    size_z = 10
    rng = np.random.default_rng(0)
    data = rng.poisson(MEAN_VAL, size=(size_z, SIZE_XY, SIZE_XY)).astype(np.uint8)

    # write the image data
    store = parse_url(path, mode="w").store
    root = zarr.group(store=store)
    write_image(
        image=data,
        group=root,
        axes="zyx",
        storage_options=dict(chunks=(1, SIZE_XY, SIZE_XY)),
    )
    return path


@pytest.fixture(scope="session")
def test_zyx_trivial(tmp_path_factory):
    path = tmp_path_factory.mktemp("test_zyx_trivial.zarr", numbered=False)

    size_z = 1
    rng = np.random.default_rng(0)
    data = rng.poisson(MEAN_VAL, size=(size_z, SIZE_XY, SIZE_XY)).astype(np.uint8)

    # write the image data
    store = parse_url(path, mode="w").store
    root = zarr.group(store=store)
    write_image(
        image=data,
        group=root,
        axes="zyx",
        storage_options=dict(chunks=(1, SIZE_XY, SIZE_XY)),
    )
    return path


@pytest.fixture(scope="session")
def test_tczyx(tmp_path_factory):
    path = tmp_path_factory.mktemp("test_tczyx.zarr", numbered=False)

    size_z = 10
    size_t = 2
    size_c = 3
    rng = np.random.default_rng(0)
    data = rng.poisson(
        MEAN_VAL, size=(size_t, size_c, size_z, SIZE_XY, SIZE_XY)
    ).astype(np.uint8)

    # write the image data
    store = parse_url(path, mode="w").store
    root = zarr.group(store=store)
    write_image(
        image=data,
        group=root,
        axes="tczyx",
        storage_options=dict(chunks=(1, 1, 1, SIZE_XY, SIZE_XY)),
    )
    return path


@pytest.fixture(scope="session")
def test_tczyx_z_trivial(tmp_path_factory):
    path = tmp_path_factory.mktemp("test_tczyx_z_trivial.zarr", numbered=False)

    size_z = 1
    size_t = 2
    size_c = 3
    rng = np.random.default_rng(0)
    data = rng.poisson(
        MEAN_VAL, size=(size_t, size_c, size_z, SIZE_XY, SIZE_XY)
    ).astype(np.uint8)

    # write the image data
    store = parse_url(path, mode="w").store
    root = zarr.group(store=store)
    write_image(
        image=data,
        group=root,
        axes="tczyx",
        storage_options=dict(chunks=(1, 1, 1, SIZE_XY, SIZE_XY)),
    )
    return path


def test_zarr_with_no_channels(test_zyx):
    # Not sure why we have to explicit about factory here
    d = load_data(test_zyx, factory=read_ome_zarr)
    assert d.shape == (10, SIZE_XY, SIZE_XY)
    assert len(d.main_components) == 1
    assert d.main_components[0].label == "value"


def test_zarr_with_no_channels_trivial(test_zyx_trivial):
    d = load_data(test_zyx_trivial, factory=read_ome_zarr)
    assert d.shape == (SIZE_XY, SIZE_XY)
    assert len(d.main_components) == 1
    assert d.main_components[0].label == "value"


def test_zarr_with_all_dimensions(test_tczyx):
    d = load_data(test_tczyx, factory=read_ome_zarr)
    assert d.shape == (2, 10, SIZE_XY, SIZE_XY)
    assert len(d.main_components) == 3
    assert d.main_components[0].label == "ch0"


def test_zarr_with_all_dimensions_z_trivial(test_tczyx_z_trivial):
    d = load_data(test_tczyx_z_trivial, factory=read_ome_zarr)
    assert d.shape == (2, SIZE_XY, SIZE_XY)
    assert len(d.main_components) == 3
    assert d.main_components[0].label == "ch0"
