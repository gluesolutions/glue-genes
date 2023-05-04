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

        self._reduced_res_data_sets = [Data(**x) for x in all_resolutions[1:]]

        # Assumes that data is listed from highest resolution to lowest resolution
        # This is fixed 2x, but we could trivially derive this from the data
        # Our assumption is that there is a consistent downscaling across
        # all relevant dimensions. This probably does not have to be true.
        self._resolutions = np.array(2 ** np.arange(len(all_resolutions)))

    def compute_fixed_resolution_buffer(self, *args, **kwargs):
        # This might only be getting the *first* view, and sometimes full_view could be integers?

        full_view = args[0]

        print(full_view)
        x = full_view[0]
        y = full_view[1]

        # Gets the index of the smallest Data object that has sufficient resolution for the request.
        # If we request a higher resolution than we have present in the data, this just gets the highest resolution available

        if isinstance(x, tuple):
            x_res = (
                abs(x[1] - x[0]) / x[2] / 2
            )  # Factor of 2 because we want to err on getting higher resolution
            x_res_mask = np.ma.masked_greater(self._resolutions, x_res)
            xx = np.ma.argmax(x_res_mask)
        else:  # This is theoretically possible, but what should we properly do?
            xx = 0
        if isinstance(y, tuple):
            y_res = abs(y[1] - y[0]) / y[2] / 2
            y_res_mask = np.ma.masked_greater(self._resolutions, y_res)
            yy = np.ma.argmax(y_res_mask)
        else:
            yy = 0

        print(f"{xx=}")
        print(f"{yy=}")
        b = min(xx, yy)
        # The only really remaining problem is whether we can look up the components in these datasets
        # if we haven't added to them the data collection object, do they have the proper links defined?

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
                old_x, old_y = full_view
                full_view = [
                    (old_x[0] / 2, old_x[1] / 2, int(old_x[2] / 2)),
                    (old_y[0] / 2, old_y[1] / 2, int(old_y[2] / 2)),
                ]
            # print(f'{full_view=}')
            for new_comp, old_comp in zip(reduced_data.components, self.components):
                if target_cid == old_comp:
                    kwargs["target_cid"] = new_comp
                    break
            target_data = kwargs.pop("target_data", None)
            kwargs["target_data"] = reduced_data
            print(kwargs)
            frb = compute_fixed_resolution_buffer(reduced_data, full_view, **kwargs)
            print(f"{frb.shape}")
            return frb
