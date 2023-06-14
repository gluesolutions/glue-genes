def setup():
    from glue.config import qt_client

    from .anndata_factory import read_anndata  # noqa
    from .menubar_plugin import diff_gene_exp_plugin  # noqa
    from .menubar_plugin import enrichrpy_plugin  # noqa
    from .menubar_plugin import summarize_gene_subset_exp_plugin  # noqa
    from .qtl_viewer.viewer import QTLViewer  # noqa

    qt_client.add(QTLViewer)

    d3_dozen = [
        "#1E77B3",
        "#FF7E0F",
        "#2BA02A",
        "#D62628",
        "#9367BC",
        "#8C564B",
        "#E277C1",
        "#7E7E7E",
        "#BCBC21",
        "#17BDCF",
        "#3A0182",
        "#004201",
    ]
    d3_twenty = [
        "#1E77B3",
        "#FF7E0F",
        "#2BA02A",
        "#D62628",
        "#9367BC",
        "#8C564B",
        "#E277C1",
        "#7E7E7E",
        "#BCBC21",
        "#17BDCF",
        "#3A0182",
        "#004201",
        "#0DFFA8",
        "#5D003F",
        "#BCBCFF",
        "#D8AFA1",
        "#B80080",
        "#004D52",
        "#6B6400",
        "#7C0100",
        "#6027FE",
    ]

    import matplotlib.cm as cm
    from glue.config import colormaps
    from matplotlib.colors import ListedColormap

    # Add some categorical colormaps

    colormaps.add("d3_dozen", ListedColormap(d3_dozen, name="d3_dozen"))
    colormaps.add("d3_twenty", ListedColormap(d3_twenty, name="d3_twenty"))

    colormaps.add("Pastel1", cm.Pastel1)
    colormaps.add("Paired", cm.Paired)
    colormaps.add("tab20", cm.tab20)
