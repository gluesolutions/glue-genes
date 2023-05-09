import numpy as np
from glue.core.data import Data
from glue.core.fixed_resolution_buffer import compute_fixed_resolution_buffer


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

        if len(all_resolutions) > 1:
            self._reduced_res_data_sets = [Data(**x) for x in all_resolutions[1:]]
        else:
            self._reduced_res_data_sets = []
        # Assumes that data is listed from highest resolution to lowest resolution
        # This is fixed 2x, but we could trivially derive this from the data
        # Our assumption is that there is a consistent downscaling across
        # all relevant dimensions. This probably does not have to be true.
        self._resolutions = np.array(2 ** np.arange(len(all_resolutions)))

    def compute_fixed_resolution_buffer(self, *args, **kwargs):
        """
        Get a fixed_resolution_buffer from the smallest resolution dataset that we can
        glue flips x/y relative to ome-zarr, so we read in the data as
        tczxy

        """
        # This might only be getting the *first* view, and sometimes full_view could be integers?

        full_view = args[0]  #  This should be the only thing in args
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
        # If we request a higher resolution than we have present in the data, this just gets the highest resolution available

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

        # print(f"{xx=}")
        # print(f"{yy=}")
        b = min(xx, yy)  #  Use the highest resolution needed for either x or y

        if b == 0:
            # print(f'{full_view=}')
            frb = compute_fixed_resolution_buffer(self, full_view, **kwargs)
            # print(f'{frb.shape}')
            return frb
        else:
            # print(f'{full_view=}')
            # Is the cache working properly?
            reduced_data = self._reduced_res_data_sets[b - 1]
            target_cid = kwargs.pop("target_cid", None)

            for r in range(b):  # for each reduction we resize full_view by 2
                # TODO: Simplify this logic for different number of dimensions
                if len(full_view) == 2:
                    old_x, old_y = full_view
                    full_view = [
                        (old_x[0] / 2, old_x[1] / 2, max(1, int(old_x[2]))),
                        (old_y[0] / 2, old_y[1] / 2, max(1, int(old_y[2]))),
                    ]
                elif len(full_view) == 3:
                    z, old_x, old_y = full_view
                    full_view = [
                        z,
                        (old_x[0] / 2, old_x[1] / 2, max(1, int(old_x[2]))),
                        (old_y[0] / 2, old_y[1] / 2, max(1, int(old_y[2]))),
                    ]
                elif len(full_view) == 4:
                    t, z, old_x, old_y = full_view
                    full_view = [
                        t,
                        z,
                        (old_x[0] / 2, old_x[1] / 2, max(1, int(old_x[2]))),
                        (old_y[0] / 2, old_y[1] / 2, max(1, int(old_y[2]))),
                    ]

            # print(f'{full_view=}')
            for new_comp, old_comp in zip(reduced_data.components, self.components):
                if target_cid == old_comp:
                    kwargs["target_cid"] = new_comp
                    break
            # target_data = kwargs.pop("target_data", None)
            kwargs["target_data"] = reduced_data
            # print(full_view)
            # print(kwargs)
            frb = compute_fixed_resolution_buffer(reduced_data, full_view, **kwargs)
            # print(f"{frb.shape}")
            return frb
