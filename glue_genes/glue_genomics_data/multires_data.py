import numpy as np
from glue.core.component_link import ComponentLink
from glue.core.data import Data
from glue.core.exceptions import IncompatibleAttribute
from glue.core.fixed_resolution_buffer import compute_fixed_resolution_buffer
from glue.core.joins import get_mask_with_key_joins


class ReducedResolutionData(Data):
    """
    A simple data class that represents reduced resolution versions of
    data in a MultiResolutionData object
    """

    def __init__(self, label="", coords=None, parent=None, scale_factor=1, **kwargs):
        super().__init__(label=label, coords=coords, **kwargs)
        self.parent = parent
        self.scale_factor = scale_factor

    def get_data(self, cid, view=None):
        if isinstance(cid, ComponentLink):
            return cid.compute(self, view)

        if cid in self._components:
            comp = self._components[cid]
        elif cid in self.parent._components:
            comp = self.parent._components[cid]
        elif cid in self._externally_derivable_components:
            comp = self._externally_derivable_components[cid]
        elif cid in self.parent._externally_derivable_components:
            comp = self.parent._externally_derivable_components[cid]
        else:
            raise IncompatibleAttribute(cid)

        if view is not None:
            result = comp[view]
        else:
            result = comp.data

        return result

    def get_mask(self, subset_state, view=None):
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
            for new_comp, old_comp in zip(reduced_data.components, self.components):
                if target_cid == old_comp:
                    kwargs["target_cid"] = new_comp
                    break
            target_data = kwargs.pop("target_data", None)
            kwargs["target_data"] = reduced_data
            print(full_view)
            print(kwargs)
            frb = compute_fixed_resolution_buffer(reduced_data, full_view, **kwargs)
            print(f"{frb.shape}")
            return frb
