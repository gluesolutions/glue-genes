"""
A simple geojson parser for glue that does not require geopandas.
This is designed to work primarily with QuPath annotation files.

QuPath GeoJSON defines (0,0) as the upper left corner of the image,
while glue defines (0,0) as the lower left corner of the image.

We need to tackle the transformation of ExtendedComponents
in order to deal with this.

We cannot fix this in the data loader as there is no reference
to the size of the image (which would be necessary to do the transform).
"""

import json
from shapely.geometry import shape
import numpy as np
from glue.core.data import RegionData
from glue.core.component import ExtendedComponent
from glue.config import data_factory


def is_geojson(filename):
    return filename.endswith(".geojson")


@data_factory("GeoJSON data loader", is_geojson, priority=100)
def read_geojson(filename):
    features = json.load(open(filename))
    data = RegionData()

    geometries = []
    properties = []
    center_x = []
    center_y = []
    for feature in features["features"]:
        geom = shape(feature["geometry"])
        rep_point = geom.point_on_surface()
        center_x.append(rep_point.x)
        center_y.append(rep_point.y)
        geometries.append(geom)
        properties.append(feature["properties"])
    data.add_component(center_x, "x")
    data.add_component(center_y, "y")

    geom_comp = ExtendedComponent(
        np.asarray(geometries),
        center_comp_ids=[
            data.id["x"],
            data.id["y"],
        ],
    )
    data.add_component(geom_comp, label="regions")

    property_keys = {k for d in properties for k in d}
    for prop_key in property_keys:
        new_col = [x.get(prop_key, None) for x in properties]
        try:
            data.add_component(new_col, prop_key)
        except ValueError as e:
            print(
                f"The property {prop_key} was not added to the Data object because: {e}"
            )
            pass
    return data
