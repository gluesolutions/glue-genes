import numpy as np
from glue.core.component_id import ComponentID
from glue.core.component_link import ComponentLink
from glue.core.data import Data
from glue.core.exceptions import IncompatibleAttribute
from glue.core.fixed_resolution_buffer import compute_fixed_resolution_buffer
from glue.core.joins import get_mask_with_key_joins


class ReducedResolutionData(Data):
    """
    A simple data class that represents reduced resolution versions of
    data in a MultiResolutionData object

    Parameters
    ----------
    parent : `~.MultiResolutionData`
        The parent MulitResolutionData that this object is a reducted-resolution
        version of.
    scale_factor : `int`
        The scale factor between this data and parent. Currently only a single
        value, since it must be the same in all dimensions

    """

    def __init__(self, label="", coords=None, parent=None, scale_factor=1, **kwargs):
        super().__init__(label=label, coords=coords, **kwargs)
        self._parent_data = parent
        self._cid_to_parent_cid = {}
        self._parent_cid_to_cid = {}
        self.scale_factor = scale_factor

        # Construct a list of original pixel component IDs
        self._parent_pixel_cids = []
        for idim in range(self._parent_data.ndim):
            self._parent_pixel_cids.append(self._parent_data.pixel_component_ids[idim])

        # Construct a list of original world component IDs
        self._parent_world_cids = []
        if len(self._parent_data.world_component_ids) > 0:
            idim_new = 0
            for idim in range(self._parent_data.ndim):
                if self._indices[idim] is None:
                    self._cid_to_parent_cid[
                        self.world_component_ids[idim_new]
                    ] = self._parent_data.world_component_ids[idim]
                    idim_new += 1
        # Yeah, but this isn't quite right. We actually have different
        # ComponentIDs and we need to set this dictionary up when we init
        # The MultiResolutionData object
        # _ = self.main_components  # Just to trigger this code to run

    # @property
    # def main_components(self):
    #    main = []
    #    for cid in self._parent_data.main_components:
    #        if cid not in self._parent_cid_to_cid:
    #            cid_new = ComponentID(label=cid.label, parent=self)
    #            self._parent_cid_to_cid[cid] = cid_new
    #            self._cid_to_parent_cid[cid_new] = cid
    #        main.append(self._parent_cid_to_cid[cid])
    #    return main

    def _translate_full_cid_to_reduced_cid(self, cid):
        """
        This translates the full resolution cid to the reduced resolution cid
        """
        if cid in self._parent_pixel_cids:
            cid = self.pixel_component_ids[cid.axis]
        elif cid in self._parent_cid_to_cid:
            cid = self._parent_cid_to_cid[cid]
        return cid

    def _translate_reduced_cid_to_full_cid(self, cid):
        """
        This translates the reduced resolution cid to the full resolution cid
        """
        if cid in self.pixel_component_ids:
            cid = self._parent_pixel_cids[cid.axis]
        elif cid in self._cid_to_parent_cid:
            cid = self._cid_to_parent_cid[cid]
        return cid

    # def get_data(self, cid, view=None):
    #    """
    #    Sometimes
    #    """
    #    cid = self._translate_reduced_cid_to_full_cid(cid)
    #    return self._parent_data.get_data(cid, view=view)

    def get_mask(self, subset_state, view=None):
        """
        We need to either translate subset_state rois
        or use pixel/world components in the original
        coordinate system.

        The following basically works for polygon selections
        but will break for range selections, and anything
        with a categorical. Also probably any composite subset.

        We need to do something more complicated for those.

        This is adjusting the subset_state.roi into the coordinates
        of the ReducedResolutionData. Probably a better thing is
        to maintain two set of coordinates for our ReducedResolutionData
        perhaps in a dictionary and then translate any attributes that are
        in the subset_state.attributes (will this always be sufficient?)
        to use these other attributes. See IndexedData

        """

        # import pdb

        # pdb.set_trace()

        new_atts = []
        for att in subset_state.attributes:
            for attribute in self._components:
                if (
                    attribute.label == att.label
                ):  # This is a weak check, we should do better by maintaining our own lookup. See IndexedData
                    new_atts.append(attribute)
        # import pdb

        def scale_roi(x):
            return x / self.scale_factor

        # pdb.set_trace()
        # subset_state_reduced.attributes = tuple(new_atts)
        subset_state_reduced = subset_state.copy()

        try:
            subset_state_reduced.xatt = new_atts[0]
        except IndexError:
            pass
        try:
            subset_state_reduced.yatt = new_atts[1]
        except IndexError:
            pass
        subset_state_reduced.roi = subset_state.roi.transformed(
            xfunc=scale_roi, yfunc=scale_roi
        )
        print(f"{self.scale_factor=}")
        print(subset_state.roi)
        print(subset_state_reduced.roi)
        # TODO: We should probably chain these pretransforms in case we already have one
        # subset_state_reduced.pretransform = scale_down
        try:
            array = subset_state_reduced.to_mask(self, view=view)
            return array
        except IncompatibleAttribute:
            return get_mask_with_key_joins(
                self, self._key_joins, subset_state_reduced, view=view
            )


class MultiResolutionData(Data):
    """
    A data class that supports multi-resolution n-dimensional datasets.
    If these multiple resolutions are accessed as DaskArrays this provides
    a memory-efficient way to browse very high-resolution images.

    The basic scheme is to use the highest resolution dataset for setting
    up components but to store the lower-resolution versions as a list of
    Data objects. Then, whenever we need a chunk of the data, we calculate
    whether we can use one of the lower-resolution images instead. The
    lower-resolution versions are not part of the DataCollection and so
    we can't use the general linking framework, but as they are all different
    views of the same data we can do this manually.

    Currently this class assumes OME-Zarr datasets where the multi-resolution
    datasets are stored from highest to lowest resolution with each lower resolution
    version downsampled by a factor of two.

    """

    def __init__(self, label="", coords=None, all_resolutions=[], **kwargs):
        super(MultiResolutionData, self).__init__(label=label, coords=coords, **kwargs)
        self.scale = 2

        if len(all_resolutions) > 1:
            self._reduced_res_data_sets = [
                ReducedResolutionData(
                    **x, parent=self, scale_factor=self.scale ** (i + 1)
                )
                for i, x in enumerate(all_resolutions[1:])
            ]
        else:
            self._reduced_res_data_sets = []
        # Assumes that data is listed from highest resolution to lowest resolution
        # This is fixed 2x, but we could trivially derive this from the data
        # Our assumption is that there is a consistent downscaling across
        # all relevant dimensions. This probably does not have to be true.
        self._resolutions = np.array(self.scale ** np.arange(len(all_resolutions)))
        self.setup_cid_lookups()

    def setup_cid_lookups(self):
        for reduced_data in self._reduced_res_data_sets:
            for r_compid, compid in zip(
                reduced_data.main_components, self.main_components
            ):
                reduced_data._cid_to_parent_cid[r_compid] = compid
                reduced_data._parent_cid_to_cid[compid] = r_compid

    def compute_fixed_resolution_buffer(self, *args, **kwargs):
        """
        Get a fixed_resolution_buffer from the smallest resolution dataset that we can
        glue flips x/y relative to ome-zarr, so we read in the data as
        tczxy

        """
        # This might only be getting the *first* view, and sometimes full_view could be integers?

        full_view = args[0]  # This should be the only thing in args
        # view is t,z,x,y where t and z are optional
        if len(full_view) == 4:
            t = full_view[0]
            z = full_view[1]
            y = full_view[3]
            x = full_view[2]
        elif len(full_view) == 3:
            z = full_view[0]
            y = full_view[2]
            x = full_view[1]
        else:
            y = full_view[1]
            x = full_view[0]

        # Gets the index of the smallest Data object that has sufficient resolution for the request.
        # If we request a higher resolution than we have present in the data,
        # this just gets the highest resolution available

        if isinstance(x, tuple):
            x_res = abs(x[1] - x[0]) / x[2]
            x_res_mask = np.ma.masked_greater(self._resolutions, x_res)
            xx = np.ma.argmax(x_res_mask)
        else:  # This generally should not happen
            xx = 0
        if isinstance(y, tuple):
            y_res = abs(y[1] - y[0]) / y[2]
            y_res_mask = np.ma.masked_greater(self._resolutions, y_res)
            yy = np.ma.argmax(y_res_mask)
        else:
            yy = 0

        print(f"{xx=}")
        print(f"{yy=}")
        b = min(xx, yy)  # Use the highest resolution needed for either x or y

        if b == 0:
            print(f"{full_view=}")
            frb = compute_fixed_resolution_buffer(self, full_view, **kwargs)
            print(f"{frb.shape}")
            return frb
        else:
            print(f"{full_view=}")
            # Is the cache working properly?
            reduced_data = self._reduced_res_data_sets[b - 1]
            target_cid = kwargs.pop("target_cid", None)

            for r in range(b):  # for each reduction we resize full_view by 2
                # TODO: Simplify this logic for different number of dimensions
                if len(full_view) == 2:
                    old_x, old_y = full_view
                    full_view = [
                        (
                            old_x[0] / self.scale,
                            old_x[1] / self.scale,
                            max(1, int(old_x[2])),
                        ),
                        (
                            old_y[0] / self.scale,
                            old_y[1] / self.scale,
                            max(1, int(old_y[2])),
                        ),
                    ]
                elif len(full_view) == 3:
                    z, old_x, old_y = full_view
                    full_view = [
                        z,
                        (
                            old_x[0] / self.scale,
                            old_x[1] / self.scale,
                            max(1, int(old_x[2])),
                        ),
                        (
                            old_y[0] / self.scale,
                            old_y[1] / self.scale,
                            max(1, int(old_y[2])),
                        ),
                    ]
                elif len(full_view) == 4:
                    t, z, old_x, old_y = full_view
                    full_view = [
                        t,
                        z,
                        (
                            old_x[0] / self.scale,
                            old_x[1] / self.scale,
                            max(1, int(old_x[2])),
                        ),
                        (
                            old_y[0] / self.scale,
                            old_y[1] / self.scale,
                            max(1, int(old_y[2])),
                        ),
                    ]

            print(f"{full_view=}")
            print(reduced_data._parent_cid_to_cid)
            # This is not doing it for some reason...
            kwargs["target_cid"] = reduced_data._translate_full_cid_to_reduced_cid(
                target_cid
            )
            target_data = kwargs.pop("target_data", None)
            kwargs["target_data"] = reduced_data
            print(full_view)
            print(kwargs)
            frb = compute_fixed_resolution_buffer(reduced_data, full_view, **kwargs)
            print(f"{frb.shape}")
            return frb
