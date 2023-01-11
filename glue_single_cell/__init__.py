def setup():
    from .anndata_factory import read_anndata  # noqa
    from .anndata_factory import setup_anndata  # noqa
    from .menubar_plugin import diff_gene_exp_plugin  # noqa
    from .menubar_plugin import pca_subset_exp_plugin  # noqa
    from .menubar_plugin import enrichrpy_plugin  # noqa
    from .qtl_viewer.viewer import QTLViewer  # noqa

    from glue.config import qt_client

    qt_client.add(QTLViewer)

    # Add colormaps we can use when we use the pca_subset_exp_plugin
    # to maps gene expression to a summary over cells. A summary generated
    # from a red subset can then be displayed using the 'Reds' colormap
    # This does not happen automatically... at least not yet

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

    from glue.config import colormaps
    import matplotlib.cm as cm
    from matplotlib.colors import ListedColormap

    colormaps.add("Reds", cm.Reds)
    colormaps.add("Greens", cm.Greens)
    colormaps.add("Blues", cm.Blues)
    colormaps.add("Purples", cm.Purples)
    colormaps.add("Oranges", cm.Oranges)
    # Add some categorical colormaps

    colormaps.add("d3_dozen", ListedColormap(d3_dozen, name="d3_dozen"))
    colormaps.add("d3_twenty", ListedColormap(d3_twenty, name="d3_twenty"))

    colormaps.add("Pastel1", cm.Pastel1)
    colormaps.add("Paired", cm.Paired)
    colormaps.add("tab20", cm.tab20)
