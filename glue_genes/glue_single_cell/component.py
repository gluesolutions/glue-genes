from glue.core.component import Component
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap
from glue.config import colormaps
import warnings

__all__ = ['SyncComponent']


class SyncComponent(Component):
    """
    A component that is automatically updated by something else within glue,
    currently one or more subset definitions.

    TODO: This _could_ be where the Listener lives that keeps it up-to-date. Currently
          this is a GeneSummaryListener added to the Xarray dataset

    Parameters
    ----------
    preferred_cmap : `str` or :class:`~matplotlib.colors.Colormap`, optional
        A colormap to be used as the preferred colormap, by name or instance. Default is `None`.
"""
    def __init__(self, data, units=None, subsets=[]):
        super().__init__(data, units=units)
        if len(subsets) > 0:
            self.cmap_name, self.preferred_cmap = self.get_cmap_from_subsets(subsets)

    def get_cmap_from_subsets(self, subsets):
        if len(subsets) == 1:
            colors = ["white", subsets[0].style.color]
            cmap_name = f"{subsets[0].label}_cmap"
        elif len(subsets) == 2:
            colors = [subsets[0].style.color, "white", subsets[1].style.color]
            cmap_name = f"{subsets[0].label}_{subsets[1].label}_cmap"
        else:
            raise NotImplementedError()
        my_cmap = LinearSegmentedColormap.from_list(cmap_name, colors)

        # my_cmap_r = my_cmap.reversed()
        colormaps.add(cmap_name, my_cmap)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            mpl.colormaps.register(cmap=my_cmap, force=True) # Do we need this?
        return (cmap_name, my_cmap)
