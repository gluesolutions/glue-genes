from glue.viewers.scatter.layer_artist import ScatterLayerArtist
from glue.viewers.scatter.state import ScatterLayerState
from glue.viewers.scatter.python_export import python_export_scatter_layer
from glue.utils import defer_draw, broadcast_to, ensure_numerical
from glue.core.exceptions import IncompatibleAttribute


from glue.viewers.scatter.layer_artist import (set_mpl_artist_cmap,
                                                DensityMapLimits,
                                                ravel_artists,
                                                InvertedNormalize)
from glue.viewers.scatter.layer_artist import (CMAP_PROPERTIES,
                                               MARKER_PROPERTIES,
                                               LINE_PROPERTIES,
                                               DENSITY_PROPERTIES,
                                               VISUAL_PROPERTIES,
                                               DATA_PROPERTIES)
import numpy as np


DATA_PROPERTIES.update(['lod_att', 'lod_thresh'])


class QTLLayerArtist(ScatterLayerArtist):
    """
    This is not going to work for a density artist
    
    Specifically, _update_data just filters the data list before
    display, which does not work for the density map mode
    """

    _layer_state_cls = ScatterLayerState
    _python_exporter = python_export_scatter_layer
    
    def __init__(self, axes, viewer_state, layer_state=None, layer=None):
    
        super(QTLLayerArtist, self).__init__(axes, viewer_state,
                                                 layer_state=layer_state, layer=layer)
        
        # Watch for changes in the viewer state which would require the
        # layers to be redrawn
        self._viewer_state.add_global_callback(self._update_scatter)
        self.state.add_global_callback(self._update_scatter)
        
        # Scatter density
        self.density_auto_limits = DensityMapLimits()
        self._set_axes(axes)
        self.errorbar_index = 2
        self.vector_index = 3
        self.lod_mask = None
        
        # NOTE: Matplotlib can't deal with NaN values in errorbar correctly, so
        # we need to prefilter values - the following variable is used to store
        # the mask for the values we keep, so that we can apply it to the color
        # See also https://github.com/matplotlib/matplotlib/issues/13799
        self._errorbar_keep = None

    @defer_draw
    def _update_scatter(self, force=False, **kwargs):
    
        if (self._viewer_state.x_att is None or
            self._viewer_state.y_att is None or
                self.state.layer is None):
            return
    
        changed = set() if force else self.pop_changed_properties()
    
        if force or len(changed & DATA_PROPERTIES) > 0:
            self._update_data()
            force = True
    
        if force or len(changed & VISUAL_PROPERTIES) > 0:
            self._update_visual_attributes(changed, force=force)


    @defer_draw
    def _update_data(self):
    
        #print(f"{self._viewer_state.lod_att}")
        # Layer artist has been cleared already
        if len(self.mpl_artists) == 0:
            return
    
        try:
            if not self.state.density_map:
                x = ensure_numerical(self.layer[self._viewer_state.x_att].ravel())
    
        except (IncompatibleAttribute, IndexError):
            # The following includes a call to self.clear()
            self.disable_invalid_attributes(self._viewer_state.x_att)
            return
        else:
            self.enable()
    
        try:
            if not self.state.density_map:
                y = ensure_numerical(self.layer[self._viewer_state.y_att].ravel())
        except (IncompatibleAttribute, IndexError):
            # The following includes a call to self.clear()
            self.disable_invalid_attributes(self._viewer_state.y_att)
            return
        else:
            self.enable()
        
        if self._viewer_state.lod_att is not None and self._viewer_state.lod_thresh is not None:
            #print("Applying lod thresholding...")
            try:
                if not self.state.density_map:
                    lod = ensure_numerical(self.layer[self._viewer_state.lod_att].ravel())
            except (IncompatibleAttribute, IndexError):
                # The following includes a call to self.clear()
                self.disable_invalid_attributes(self._viewer_state.lod_att)
                return
            else:
                self.enable()
    
            self.lod_mask = np.ma.masked_where(lod < self._viewer_state.lod_thresh, lod)
            masked_x = np.ma.masked_where(np.ma.getmask(self.lod_mask), x)
            masked_y = np.ma.masked_where(np.ma.getmask(self.lod_mask), y)
        else:
            masked_x = x
            masked_y = y
    
        if self.state.markers_visible:
    
            if self.state.density_map:
                # We don't use x, y here because we actually make use of the
                # ability of the density artist to call a custom histogram
                # method which is defined on this class and does the data
                # access.
                self.plot_artist.set_data([], [])
                self.scatter_artist.set_offsets(np.zeros((0, 2)))
            else:
                self.density_artist.set_label(None)
                if self._use_plot_artist():
                    # In this case we use Matplotlib's plot function because it has much
                    # better performance than scatter.
                    self.plot_artist.set_data(masked_x, masked_y)
                else:
                    offsets = np.vstack((masked_x, masked_y)).transpose()
                    self.scatter_artist.set_offsets(offsets)
        else:
            self.plot_artist.set_data([], [])
            self.scatter_artist.set_offsets(np.zeros((0, 2)))
            
            
    @defer_draw
    def _update_visual_attributes(self, changed, force=False):
    
        if not self.enabled:
            return
    
        if self.state.markers_visible:
    
            if self.state.density_map:
    
                if self.state.cmap_mode == 'Fixed':
                    if force or 'color' in changed or 'cmap_mode' in changed:
                        self.density_artist.set_color(self.state.color)
                        self.density_artist.set_clim(self.density_auto_limits.min,
                                                     self.density_auto_limits.max)
                elif force or any(prop in changed for prop in CMAP_PROPERTIES):
                    c = ensure_numerical(self.layer[self.state.cmap_att].ravel())
                    if self.lod_mask is not None:
                        c = np.ma.masked_where(np.ma.getmask(self.lod_mask), c)
                    set_mpl_artist_cmap(self.density_artist, c, self.state)
    
                if force or 'stretch' in changed:
                    self.density_artist.set_norm(ImageNormalize(stretch=STRETCHES[self.state.stretch]()))
    
                if force or 'dpi' in changed:
                    self.density_artist.set_dpi(self._viewer_state.dpi)
    
                if force or 'density_contrast' in changed:
                    self.density_auto_limits.contrast = self.state.density_contrast
                    self.density_artist.stale = True
    
            else:
    
                if self._use_plot_artist():
    
                    if force or 'color' in changed or 'fill' in changed:
                        if self.state.fill:
                            self.plot_artist.set_markeredgecolor('none')
                            self.plot_artist.set_markerfacecolor(self.state.color)
                        else:
                            self.plot_artist.set_markeredgecolor(self.state.color)
                            self.plot_artist.set_markerfacecolor('none')
    
                    if force or 'size' in changed or 'size_scaling' in changed:
                        self.plot_artist.set_markersize(self.state.size *
                                                        self.state.size_scaling)
    
                else:
    
                    # TEMPORARY: Matplotlib has a bug that causes set_alpha to
                    # change the colors back: https://github.com/matplotlib/matplotlib/issues/8953
                    if 'alpha' in changed:
                        force = True
    
                    if self.state.cmap_mode == 'Fixed':
                        if force or 'color' in changed or 'cmap_mode' in changed or 'fill' in changed:
                            self.scatter_artist.set_array(None)
                            if self.state.fill:
                                self.scatter_artist.set_facecolors(self.state.color)
                                self.scatter_artist.set_edgecolors('none')
                            else:
                                self.scatter_artist.set_facecolors('none')
                                self.scatter_artist.set_edgecolors(self.state.color)
                    elif force or any(prop in changed for prop in CMAP_PROPERTIES) or 'fill' in changed:
                        self.scatter_artist.set_edgecolors(None)
                        self.scatter_artist.set_facecolors(None)
                        c = ensure_numerical(self.layer[self.state.cmap_att].ravel())
                        if self.lod_mask is not None:
                            c = np.ma.masked_where(np.ma.getmask(self.lod_mask), c)

                        set_mpl_artist_cmap(self.scatter_artist, c, self.state)
                        if self.state.fill:
                            self.scatter_artist.set_edgecolors('none')
                        else:
                            self.scatter_artist.set_facecolors('none')
    
                    if force or any(prop in changed for prop in MARKER_PROPERTIES):
    
                        if self.state.size_mode == 'Fixed':
                            s = self.state.size * self.state.size_scaling
                            s = broadcast_to(s, self.scatter_artist.get_sizes().shape)
                        else:
                            s = ensure_numerical(self.layer[self.state.size_att].ravel())
                            if self.lod_mask is not None:
                                s = np.ma.masked_where(np.ma.getmask(self.lod_mask), s)

                            s = ((s - self.state.size_vmin) /
                                 (self.state.size_vmax - self.state.size_vmin))
                            # The following ensures that the sizes are in the
                            # range 3 to 30 before the final size_scaling.
                            np.clip(s, 0, 1, out=s)
                            s *= 0.95
                            s += 0.05
                            s *= (30 * self.state.size_scaling)
    
                        # Note, we need to square here because for scatter, s is actually
                        # proportional to the marker area, not radius.
                        self.scatter_artist.set_sizes(s ** 2)

        for artist in [self.scatter_artist, self.plot_artist,
                       self.vector_artist, self.line_collection,
                       self.density_artist]:
    
            if artist is None:
                continue
    
            if force or 'alpha' in changed:
                artist.set_alpha(self.state.alpha)
    
            if force or 'zorder' in changed:
                artist.set_zorder(self.state.zorder)
    
            if force or 'visible' in changed:
                # We need to hide the density artist if it is not needed because
                # otherwise it might still show even if there is no data as the
                # neutral/zero color might not be white.
                if artist is self.density_artist:
                    artist.set_visible(self.state.visible and
                                       self.state.density_map and
                                       self.state.markers_visible)
                else:
                    artist.set_visible(self.state.visible)
        if self._use_plot_artist():
            self.scatter_artist.set_visible(False)
        else:
            self.plot_artist.set_visible(False)
        self.redraw()