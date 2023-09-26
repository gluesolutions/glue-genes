from glue.core.component import Component
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap
from glue.config import colormaps
import warnings

__all__ = ["SyncComponent"]


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

    def __init__(self, data, units=None, subsets=[], preferred_cmap=None):
        super().__init__(data, units=units)
        self.subsets = subsets
        if preferred_cmap is not None:
            if isinstance(preferred_cmap, str):
                try:
                    self.preferred_cmap = mpl.colormaps.get_cmap(preferred_cmap)
                    self.cmap_name = preferred_cmap
                except ValueError:
                    if len(subsets) > 0:
                        self.cmap_name, self.preferred_cmap = self.get_cmap_from_subsets(subsets)
                    else:
                        self.preferred_cmap = None
                        self.cmap_name = None
            elif isinstance(preferred_cmap, mpl.colors.Colormap):
                self.preferred_cmap = preferred_cmap
                self.cmap_name = preferred_cmap.name
        elif len(subsets) > 0:
            self.cmap_name, self.preferred_cmap = self.get_cmap_from_subsets(subsets)
        else:
            self.preferred_cmap = None
            self.cmap_name = None

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
            mpl.colormaps.register(cmap=my_cmap, force=True)  # Do we need this?
        return (cmap_name, my_cmap)

    def __gluestate__(self, context):
        if isinstance(self.preferred_cmap, str):
            return dict(data=context.id(self._data),
                        units=context.id(self._units),
                        subsets=context.id(self.subsets),
                        preferred_cmap=self.preferred_cmap)
        else:
            return dict(data=context.id(self._data),
                        units=context.id(self._units),
                        subsets=context.id(self.subsets),
                        preferred_cmap=self.cmap_name)

    @classmethod
    def __setgluestate__(cls, rec, context):
        return cls(context.object(rec['data']),
                   units=context.object(rec['units']),
                   subsets=context.object(rec['subsets']),
                   preferred_cmap=rec['preferred_cmap'])
