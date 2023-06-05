def setup():
    from .bed_factory import read_bed  # noqa
    from .bigwig_factory import read_bigwig  # noqa
    from .spaceranger_factory import read_spaceranger_directory  # noqa
    from .spaceranger_factory import spatialtrans_autolink  # noqa
    from .zarr_factory import read_zarr  # noqa
