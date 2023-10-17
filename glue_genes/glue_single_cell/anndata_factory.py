from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from glue.config import data_factory, autolinker
from glue.core import Data
from glue.core.link_helpers import JoinLink
from glue.core.component import ExtendedComponent
from glue.core.component_id import PixelComponentID
from glue.core.data_region import RegionData
from glue_qt.utils import set_cursor_cm
from qtpy.QtCore import Qt
from shapely.geometry import Point

from .data import DataAnnData
from .qt.load_data import LoadDataDialog

__all__ = [
    "read_anndata",
    "join_anndata_on_keys",
]


def df_to_data(obj, label=None, skip_components=[]):
    result = Data(label=label)
    for c in obj.columns:
        if c not in skip_components:
            result.add_component(obj[c], str(c))
    return result


def is_anndata(filename, **kwargs):
    return filename.endswith(".h5ad") or filename.endswith(".loom")


def join_anndata_on_keys(datasets):
    """
    Use join_on_key to stitch the various components on an anndata
    dataset back together. We join on the pixel ids and indices
    because it is far faster to do this than to join on the string
    names for genes and cells.

    TODO: Test the assumption that matching on pixel ids is really
    faster than matching on string names. Components names for
    cells (obs) and genes (vars) are standard? in AnnData objects
    so we could join on these, which is more intuitive in the UI.
    """
    varset = {d for d in datasets if d.meta["join_on_var"] is True}
    obsset = {d for d in datasets if d.meta["join_on_obs"] is True}

    for dataset in datasets:
        if dataset.meta["anndatatype"] == "X Array":
            for d in obsset:
                if d.meta["anndatatype"] != "X Array":
                    dataset.join_on_key(d, "Pixel Axis 0 [y]", "Pixel Axis 0 [x]")
            for d in varset:
                if d.meta["anndatatype"] != "X Array":
                    dataset.join_on_key(d, "Pixel Axis 1 [x]", "Pixel Axis 0 [x]")
        elif dataset.meta["anndatatype"] == "obs Array":
            for d in obsset:
                # Do not join to self or X Array again
                if (d.meta["anndatatype"] != "X Array") and (
                    d.meta["anndatatype"] != "obs Array"
                ):
                    dataset.join_on_key(d, "Pixel Axis 0 [x]", "Pixel Axis 0 [x]")
        elif dataset.meta["anndatatype"] == "var Array":
            for d in varset:
                # Do not join to self or X Array again
                if (d.meta["anndatatype"] != "X Array") and (
                    d.meta["anndatatype"] != "var Array"
                ):
                    dataset.join_on_key(d, "Pixel Axis 0 [x]", "Pixel Axis 0 [x]")
    return datasets


@data_factory("AnnData Loader", is_anndata, priority=999)
def read_anndata(
    file_name,
    skip_dialog=False,
    skip_components=[],
    subsample=False,
    subsample_factor=1,
    try_backed=False,
):
    """
    Use AnnData to read a single-cell type data file into three linked glue Data objects.

    AnnData objects are mapped into three glue Data objects that are linked by JoinLink

    Currently supports .h5ad and .loom files, but reading in
    backed mode is only supported for .h5ad because of underlying
    limitations with AnnData.

    Parameters
    ----------
    file_name: str
        The file to read
    skip_dialog: bool, optional
        Whether to skip the GUI module (useful for scripts and testing)
    skip_components: list, optional
        The names of columns in obs/var/obsm/varm to NOT load
    subsample: bool, optional
        Whether to subsample on obs to reduce file size
    subsample_factor: int, optional
        If specified, and `subsample`, reduce the size of obs by this factor
    try_backed: bool, optional
        Attempt to use disk-based access to the data. If AnnData fails to load the file this
        way it will be loaded into memory.

    """
    basename = Path(file_name).stem

    if not skip_dialog:
        with set_cursor_cm(Qt.ArrowCursor):
            load_dialog = LoadDataDialog(filename=file_name)
            if load_dialog.exec_():
                skip_components = load_dialog.skip_components
                subsample = load_dialog.subsample
                try_backed = load_dialog.try_backed
                subsample_factor = load_dialog.subsample_factor
            else:
                return []

    if try_backed:
        try:
            adata = sc.read(file_name, sparse=True, backed="r")
            backed = True
        except OSError:
            adata = sc.read(file_name, sparse=True, backed=False)
            backed = False
    else:
        adata = sc.read(file_name, sparse=True, backed=False)
        backed = False

    if "spatial" in adata.uns_keys():
        make_spatial_components = True
    else:
        make_spatial_components = False

    return translate_adata_to_DataAnnData(
        adata,
        subsample=subsample,
        backed=backed,
        basename=basename,
        file_name=file_name,
        skip_components=skip_components,
        subsample_factor=subsample_factor,
        try_backed=try_backed,
        make_spatial_components=make_spatial_components,
    )


def translate_adata_to_DataAnnData(
    adata,
    subsample=False,
    backed=False,
    basename="",
    file_name="",
    skip_components=[],
    subsample_factor=1,
    try_backed=False,
    skip_joins=False,
    make_spatial_components=False,
):
    list_of_data_objs = []

    adata.var_names_make_unique()

    if subsample:
        adata = sc.pp.subsample(
            adata, fraction=subsample_factor, copy=True, random_state=0
        )

    if backed:
        XData = DataAnnData(
            Xarray=adata.X, full_anndata_obj=adata, backed=backed, label=f"{basename}_X"
        )
    else:
        XData = DataAnnData(Xarray=adata.X, backed=backed, label=f"{basename}_X")

    XData.meta["orig_filename"] = basename

    if make_spatial_components:
        library_id = list(adata.uns["spatial"].keys())[0]
        basename = (
            library_id  # This is often a nicer name, but maybe this should be optional?
        )
        # Want radius
        # We should not assume this much about the structure of spatial
        # but this is okay for now...
        scale_fac = adata.uns["spatial"][library_id]["scalefactors"]
        # hi_res_scale_fac = scale_fac["tissue_hires_scalef"]
        # Want radius
        spot_size = scale_fac["spot_diameter_fullres"] / 2.0  # * hi_res_scale_fac / 2.0

    XData.meta["full_filename"] = file_name
    XData.meta["Xdata"] = XData.uuid
    XData.meta["anndatatype"] = "X Array"
    XData.meta["join_on_obs"] = True
    XData.meta["join_on_var"] = True

    # This meta-data is attached to the DataAnnData object so that
    # We can pass it to LoadLog on save/restore
    XData.meta["loadlog_skip_components"] = skip_components
    XData.meta["loadlog_subsample"] = subsample
    XData.meta["loadlog_subsample_factor"] = subsample_factor
    XData.meta["loadlog_try_backed"] = try_backed

    # uns is unstructured data on the AnnData object
    # We just store it in metadata so we can recreate
    # the AnnData object
    XData.meta["uns"] = adata.uns
    # remove rank_genes_groups from uns
    # rank_genes_groups produces a structured
    # array of objects (not strings) that is
    # really hard to sererialize.

    if "rank_genes_groups" in XData.meta["uns"]:
        del XData.meta["uns"]["rank_genes_groups"]
    list_of_data_objs.append(XData)

    # The var array is all components of the same length
    # and is stored by AnnData as a Pandas DataFrame
    try:
        var = adata.var
        var_data = df_to_data(
            var, label=f"{basename}_vars", skip_components=skip_components
        )
        var_data.add_component(adata.var_names, "var_names")
        var_data.meta["Xdata"] = XData.uuid
        var_data.meta["anndatatype"] = "var Array"
        var_data.meta["join_on_obs"] = False
        var_data.meta["join_on_var"] = True
        XData.meta["var_data"] = var_data
        list_of_data_objs.append(var_data)
    except:  # noqa E722
        pass

    for key in adata.varm_keys():
        if key not in skip_components:
            data_arr = adata.varm[key]
            # Sometimes this is a dataframe with names, and sometimes a simple np.array
            if isinstance(data_arr, pd.DataFrame):
                data_to_add = {
                    f"{key}_{col}": data_arr[col].values for col in data_arr.columns
                }
            else:
                data_to_add = {f"{key}_{i}": k for i, k in enumerate(data_arr.T)}
            for comp_name, comp in data_to_add.items():
                var_data.add_component(comp, comp_name)

    # The obs array is all components of the same length
    # and is stored by AnnData as a Pandas DataFrame
    try:
        obs = adata.obs
        obs_data = df_to_data(
            obs, label=f"{basename}_obs", skip_components=skip_components
        )
        obs_data.add_component(adata.obs_names, "obs_names")
        obs_data.meta["Xdata"] = XData.uuid
        obs_data.meta["anndatatype"] = "obs Array"
        obs_data.meta["join_on_obs"] = True
        obs_data.meta["join_on_var"] = False
        XData.meta["obs_data"] = obs_data
        list_of_data_objs.append(obs_data)
    except:  # noqa E722
        pass

    for key in adata.obsm_keys():
        if key not in skip_components:
            data_arr = adata.obsm[key]
            # Sometimes this is a dataframe with names, and sometimes a simple np.array
            if isinstance(data_arr, pd.DataFrame):
                data_to_add = {
                    f"{key}_{col}": data_arr[col].values for col in data_arr.columns
                }
            else:
                data_to_add = {f"{key}_{i}": k for i, k in enumerate(data_arr.T)}
            for comp_name, comp in data_to_add.items():
                obs_data.add_component(comp, comp_name)

    if make_spatial_components:
        obs_data = list_of_data_objs.pop()

        library_id = list(adata.uns["spatial"].keys())[0]
        image_data = adata.uns["spatial"][library_id]["images"]["hires"]
        scale_fac = adata.uns["spatial"][library_id]["scalefactors"]
        hi_res_scale_fac = scale_fac["tissue_hires_scalef"]

        new_spatial_0 = obs_data["spatial_0"]
        new_spatial_1 = (
            int(image_data.shape[0] / hi_res_scale_fac) - obs_data["spatial_1"]
        )  # flip coordinates
        spatial_0_id = obs_data.id["spatial_0"]
        spatial_1_id = obs_data.id["spatial_1"]
        obs_data.update_components({spatial_1_id: new_spatial_1})
        obs_data.update_components({spatial_0_id: new_spatial_0})

        # We need to cast the obs Data object into a RegionData object
        spots = []
        for x, y in zip(obs_data["spatial_0"], obs_data["spatial_1"]):
            spots.append(Point(x, y).buffer(spot_size))
        spot_arr = np.array(spots)

        obs_data_new = RegionData(label=obs_data.label)
        for compid in obs_data.components:
            if not isinstance(compid, PixelComponentID):
                # Use same names (with .label) but NOT same ComponentIDs!
                obs_data_new.add_component(obs_data.get_component(compid), compid.label)

        spot_comp = ExtendedComponent(
            spot_arr,
            center_comp_ids=[
                obs_data_new.id["spatial_0"],
                obs_data_new.id["spatial_1"],
            ],
        )

        obs_data_new.add_component(spot_comp, label="spots")
        obs_data_new.meta = obs_data.meta
        XData.meta["obs_data"] = obs_data_new
        list_of_data_objs.append(obs_data_new)

    # obs_data.meta['xarray_data'] = Xdata
    # var_data.meta['xarray_data'] = Xdata

    return join_anndata_on_keys(list_of_data_objs)


@autolinker("AnnData")
def anndata_autolink(data_collection):
    """
    This sets up automatic links between the components of an Anndata
    dataset, specifically join_on_key links between the X and obs/var
    datasets. This is pretty straightforward because our data loader
    already sets up the _key_joins in the dataset and this is "just"
    adding them as full GUI links.
    """
    anndata_datasets = [
        data for data in data_collection if isinstance(data, DataAnnData)
    ]
    if len(anndata_datasets) < 1:
        return []

    gui_links = []
    for data in anndata_datasets:
        for other, joins in data._key_joins.items():
            cid, cid_other = joins
            gui_link = JoinLink(
                cids1=[cid[0]], cids2=[cid_other[0]], data1=data, data2=other
            )
            if gui_link not in data_collection._link_manager._external_links:
                gui_links.append(gui_link)
    return gui_links
