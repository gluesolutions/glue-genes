import pytest
from glue.core.data_factories import load_data
from glue_genes.glue_genomics_data.tifffile_factory import read_tifffile_as_zarr
from tifffile import imwrite
import numpy as np


SIZE_XY = 256


@pytest.fixture(scope="session")
def basic_tiff(tmp_path_factory):
    path = tmp_path_factory.mktemp("data", numbered=False) / "basic_tiff.tif"

    data = np.random.randint(0, 255, (SIZE_XY, SIZE_XY, 3), "uint8")
    imwrite(path, data, photometric="rgb")
    return path


def test_basic_tiff(basic_tiff):
    d = load_data(basic_tiff, factory=read_tifffile_as_zarr)
    assert d.shape == (SIZE_XY, SIZE_XY)
    assert len(d.main_components) == 3
