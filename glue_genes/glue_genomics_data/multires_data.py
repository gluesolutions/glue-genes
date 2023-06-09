import numpy as np
from glue.core.component import CoordinateComponent
from glue.core.component_id import PixelComponentID
from glue.core.component_link import ComponentLink
from glue.core.data import Data, pixel_label
from glue.core.exceptions import IncompatibleAttribute
from glue.core.fixed_resolution_buffer import compute_fixed_resolution_buffer
from glue.core.joins import get_mask_with_key_joins


class ReducedCoordinateComponent(CoordinateComponent):
    def __init__(self, data, axis, world=False, stride=1):
        super().__init__(data, axis, world=world)
        self.stride = stride

    def _calculate(self, view=None):
        if self.world:
            return super()._calculate(view=view)
        else:
            slices = [
                slice(0, s, self.stride) for s in np.array(self.shape) * self.stride
            ]  # Use stride here to get downsampled pixel coordinates
            grids = np.broadcast_arrays(*np.ogrid[slices])
            if view is not None:
                grids = [g[view] for g in grids]
            return grids[self.axis]


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

    TODO: Check that any of this works for world coordinate components
    """

    def __init__(
        self,
        label="",
        coords=None,
        parent=None,
        scale_factor=1,
        reduced_dims=[],
        **kwargs,
    ):
        super().__init__(label=label, coords=coords, **kwargs)
        self._parent_data = parent
        self._cid_to_parent_cid = {}
        self._parent_cid_to_cid = {}
        self.scale_factor = scale_factor
        self.reduced_dims = reduced_dims
        # def scale_pix(x):
        #    return x / self.scale_factor

        # Construct a list of original pixel component IDs
        self._parent_pixel_cids = []
        self._downsampled_pixel_cids = []
        for idim in range(self._parent_data.ndim):
            self._parent_pixel_cids.append(self._parent_data.pixel_component_ids[idim])
            # For ease of referencing we produce a ReducedCoordinateComponent
            # for every dimension, but only axes in self.reduced_dims are actually
            # reduced
            if idim in self.reduced_dims:
                comp = ReducedCoordinateComponent(
                    self, idim, world=False, stride=self.scale_factor
                )
            else:
                comp = ReducedCoordinateComponent(self, idim, world=False, stride=1)

            label = pixel_label(idim, self._parent_data.ndim)
            cid = PixelComponentID(idim, "Reduced Pixel Axis %s" % label, parent=self)
            self.add_component(comp, cid)
            self._downsampled_pixel_cids.append(cid)

        # Construct a list of original world component IDs
        # TODO: test if this works at all
        self._parent_world_cids = []
        if len(self._parent_data.world_component_ids) > 0:
            idim_new = 0
            for idim in range(self._parent_data.ndim):
                if self._indices[idim] is None:
                    self._cid_to_parent_cid[
                        self.world_component_ids[idim_new]
                    ] = self._parent_data.world_component_ids[idim]
                    idim_new += 1

    def convert_full_to_reduced_cid(self, cid, downsample=False):
        """
        This translates the full resolution cid to the reduced resolution cid
        """
        if cid in self._parent_pixel_cids:
            if downsample:
                cid = self._downsampled_pixel_cids[cid.axis]
            else:
                cid = self.pixel_component_ids[cid.axis]
        elif cid in self._parent_cid_to_cid:
            cid = self._parent_cid_to_cid[cid]
        return cid

    def convert_reduced_to_full_cid(self, cid, upsample=True):
        """
        This translates the reduced resolution cid to the full resolution cid
        """
        if upsample:
            if cid in self._downsampled_pixel_cids:
                cid = self._parent_pixel_cids[cid.axis]
        if cid in self.pixel_component_ids:
            cid = self._parent_pixel_cids[cid.axis]
        elif cid in self._cid_to_parent_cid:
            cid = self._cid_to_parent_cid[cid]
        return cid

    def get_data(self, cid, view=None):
        """
        In the case were we are trying to get a pixel cid
        from the original dataset we want to return the
        reduced pixel, otherwise we do the normal thing
        and we could simplify this logic to call the parent
        method
        """

        if isinstance(cid, ComponentLink):
            return cid.compute(self, view)

        if cid in self._components:
            comp = self._components[cid]
        elif cid in self._externally_derivable_components:
            comp = self._externally_derivable_components[cid]
        else:
            try:
                cid = self.convert_full_to_reduced_cid(cid, downsample=True)
                comp = self._components[cid]
            except KeyError:
                raise IncompatibleAttribute(cid)
        if view is not None:
            result = comp[view]
        else:
            result = comp.data

        return result

    def get_mask(self, subset_state, view=None):
        """ """
        try:
            array = subset_state.to_mask(self, view=view)
            return array
        except IncompatibleAttribute:
            array = get_mask_with_key_joins(
                self, self._key_joins, subset_state, view=view
            )
            return array


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

    Since the MultiResolutionData set contains the full data, we do not
    need to override get_mask or get_data.

    """

    def __init__(
        self, label="", coords=None, all_resolutions=[], reduced_dims=[], **kwargs
    ):
        super(MultiResolutionData, self).__init__(label=label, coords=coords, **kwargs)
        self.scale = 2

        if len(all_resolutions) > 1:
            self._reduced_res_data_sets = [
                ReducedResolutionData(
                    **x,
                    parent=self,
                    scale_factor=self.scale ** (i + 1),
                    reduced_dims=reduced_dims,
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
            # print(f"{full_view=}")
            frb = compute_fixed_resolution_buffer(self, full_view, **kwargs)
            # print(f"{frb.shape}")
            return frb
        else:
            # print(f"{full_view=}")
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

            # print(f"{full_view=}")
            # print(reduced_data._parent_cid_to_cid)
            kwargs["target_cid"] = reduced_data.convert_full_to_reduced_cid(target_cid)
            _ = kwargs.pop("target_data", None)
            kwargs["target_data"] = reduced_data
            # print(full_view)
            # print(kwargs)
            frb = compute_fixed_resolution_buffer(reduced_data, full_view, **kwargs)
            # print(f"{frb.shape}")
            return frb

    def compute_statistic(self, statistic, cid, **kwargs):
        """
        random_views_for_dask_array (used if random_subset is set
        and we have a dask array) has a problem if there is structure
        in the image that corresponds to the chunk size. For biology
        imaging data with wells, this is often the case. Here we assume
        that the smallest reduced-resolution dataset will be small
        enough to just compute statistics without using random subset.
        """
        _ = kwargs.pop("random_subset", None)
        kwargs["random_subset"] = None

        for reduced_data in self._reduced_res_data_sets[::-1]:
            print(f"Trying {reduced_data}")
            print(f"With scale factor {reduced_data.scale_factor}")
            if cid:
                cid = reduced_data.convert_full_to_reduced_cid(cid)
            return reduced_data.compute_statistic(statistic, cid, **kwargs)

    #    pass

    def compute_histogram(
        self, cids, weights=None, range=None, bins=None, log=None, subset_state=None
    ):
        # def compute_histogram(self, *args, **kwargs):
        """
        We want to pass the compute_histogram calls to whatever reduced resolution dataset has
        sufficient number of points. I'm not aware of a particular heuristic here. We could just
        start from the smallest dataset and then expand if we think that there aren't enough
        points.
        """

        # Start with smallest dataset and move to larger sets
        # if we don't have a good number of points within the
        # most populous bin
        for reduced_data in self._reduced_res_data_sets[::-1]:
            print(f"Trying {reduced_data}")
            print(f"With scale factor {reduced_data.scale_factor}")

            target_cids = []
            for cid in cids:
                target_cids.append(reduced_data.convert_full_to_reduced_cid(cid))
            if weights:
                weights = reduced_data.convert_full_to_reduced_cid(weights)

            histogram = reduced_data.compute_histogram(
                target_cids,
                weights=weights,
                range=range,
                bins=bins,
                log=log,
                subset_state=subset_state,
            )
            print(histogram)
            if np.max(histogram) > 30:
                return (
                    histogram * reduced_data.scale_factor
                )  # Inflate reported bins by scale factor
