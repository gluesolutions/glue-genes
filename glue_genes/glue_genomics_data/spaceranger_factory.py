from pathlib import Path

import numpy as np
from glue.config import autolinker, data_factory
from glue.core import data_factories as df
from glue.core.component import ExtendedComponent
from glue.core.component_id import PixelComponentID
from glue.core.data import RegionData
from glue.core.link_helpers import LinkSame
from shapely.geometry import Point
from squidpy.read import visium

from glue_genes.glue_single_cell.anndata_factory import (
    join_anndata_on_keys,
    translate_adata_to_DataAnnData,
)

__all__ = ["read_spaceranger_directory"]


def is_space_ranger(filename, **kwargs):
    """
    Check if the directory is a valid spaceranger dataset
    This is very inefficient, but not sure what else to do
    """
    try:
        _ = visium(filename)
    except Exception:
        return False
    else:
        return True


@data_factory(label="Spaceranger directory", identifier=is_space_ranger)
def read_spaceranger_directory(filename, **kwargs):
    """
    Read in the data as an anndata object
    Read in the hires image
    Convert pixel locations
    Modify the obs file to have the locations
    """

    spatial_path = Path(filename) / "spatial"
    adata_obj = visium(filename, load_images=True)
    library_id = list(adata_obj.uns["spatial"].keys())[0]

    # We do the spatial components by hand
    adata_objs = translate_adata_to_DataAnnData(
        adata_obj,
        basename=library_id,
        file_name=filename,
        skip_joins=True,
        make_spatial_components=False,
    )
    hi_res_image = list(spatial_path.glob("*hires*"))
    if hi_res_image:
        # print("Trying to read some image data")
        # print(hi_res_image[0])

        image_data = df.load_data(hi_res_image[0])  # this duplicates the read from disk

        image_data.label = library_id + "_image"
        # print(adata_objs)

        scale_fac = adata_obj.uns["spatial"][library_id]["scalefactors"]
        hi_res_scale_fac = scale_fac["tissue_hires_scalef"]
        # Want radius
        spot_size = scale_fac["spot_diameter_fullres"] * hi_res_scale_fac / 2.0

        obs = adata_objs[2]

        new_spatial_0 = obs["spatial_0"] * hi_res_scale_fac
        new_spatial_1 = (
            image_data.shape[0] - obs["spatial_1"] * hi_res_scale_fac
        )  # flip
        spatial_0_id = obs.id["spatial_0"]
        spatial_1_id = obs.id["spatial_1"]

        obs.update_components({spatial_1_id: new_spatial_1})
        obs.update_components({spatial_0_id: new_spatial_0})

        # We need to cast the obs Data object into a RegionData object
        spots = []
        for x, y in zip(obs["spatial_0"], obs["spatial_1"]):
            spots.append(Point(x, y).buffer(spot_size))
        spot_arr = np.array(spots)

        obs_data_new = RegionData(label=obs.label)
        for compid in obs.components:
            if not isinstance(compid, PixelComponentID):
                # Use same names (with .label) but NOT same ComponentIDs!
                obs_data_new.add_component(obs.get_component(compid), compid.label)
        spot_comp = ExtendedComponent(
            spot_arr,
            center_comp_ids=[
                obs_data_new.id["spatial_0"],
                obs_data_new.id["spatial_1"],
            ],
        )

        obs_data_new.add_component(spot_comp, label="spots")
        obs_data_new.meta = obs.meta
        obs_data_new.meta["link_spatial_trans"] = image_data.uuid

        adata_objs[0]._key_joins = {}
        adata_objs[1]._key_joins = {}
        adata_objs[0].meta["obs_data"] = obs_data_new

    return join_anndata_on_keys([adata_objs[0]] + [adata_objs[1]] + [obs_data_new]) + [
        image_data
    ]


@autolinker("Spatial Transcriptomics")
def spatialtrans_autolink(data_collection):
    """
    This is to set up automatic links between spatial
    transcriptomics datasets when importing data.

    TODO: This is nicer than setting up a Listener with a custom
    start-up action and we should probably do it for non-spatial
    DataAnnData objects as well.
    """

    spatial_datasets = [
        data for data in data_collection if "link_spatial_trans" in data.meta
    ]

    if len(spatial_datasets) < 1:
        return []

    # We will skip over any datasets that already have links
    existing = set()
    for link in data_collection.external_links:
        existing.add((link.data1, link.data2))
    print(f"Existing links: {existing}")

    all_links = []
    for i1, data1 in enumerate(spatial_datasets):
        for data2 in data_collection:
            if (data1, data2) not in existing and (
                data2.uuid == data1.meta["link_spatial_trans"]
            ):
                try:
                    link1 = LinkSame(
                        data1.id["spatial_1"], data2.pixel_component_ids[0]
                    )
                    link2 = LinkSame(
                        data1.id["spatial_0"], data2.pixel_component_ids[1]
                    )
                except ValueError:
                    continue
                all_links.append(link1)
                all_links.append(link2)
    return all_links
