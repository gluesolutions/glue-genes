from glue.core.roi_pretransforms import ProjectionMplTransform
from glue.core.subset import RangeSubsetState, roi_to_subset_state
from glue.utils import decorate_all_methods, defer_draw, mpl_to_datetime64
from glue.viewers.matplotlib.qt.data_viewer import MatplotlibDataViewer
from glue.viewers.scatter.viewer import MatplotlibScatterMixin

from .layer_artist import QTLLayerArtist
from .qt.layer_style_editor import QTLLayerStyleEditor
from .qt.options_widget import QTLOptionsWidget
from .state import QTLViewerState

__all__ = ["QTLViewer"]


@decorate_all_methods(defer_draw)
class QTLViewer(MatplotlibScatterMixin, MatplotlibDataViewer):
    """
    A dedicated viewer for eQTL-type data.

    Adds LOD threshold filtering and chromosome boundaries on top of
    a 2D Scatter viewer. The LOD threshold filtering happens in the
    LayerArtist, while the chromosome boundaries are set here in
    _update_axes.
    """

    LABEL = "QTL Viewer"
    _layer_style_widget_cls = (
        QTLLayerStyleEditor  # We should trim options that do not make sense
    )
    _state_cls = QTLViewerState
    _options_cls = QTLOptionsWidget  # This has the LOD level widget
    _data_artist_cls = QTLLayerArtist
    _subset_artist_cls = QTLLayerArtist

    tools = ["select:rectangle", "select:xrange", "select:yrange"]

    def __init__(self, session, parent=None, state=None):
        proj = None
        MatplotlibDataViewer.__init__(
            self, session, parent=parent, state=state, projection=proj
        )
        MatplotlibScatterMixin.setup_callbacks(self)
        self.state.add_callback("species", self._update_axes)
        self.state.add_callback("pos_units", self._update_axes)

    def _update_axes(self, *args):

        if (self.state.x_att is not None) and (self.state.y_att is not None):

            chr_bounds = [0]
            chr_label_positions = []
            chr_labels = []

            current_position = 0

            names = self.state.chr_pos[self.state.species]["Names"]
            gridsize = self.state.chr_pos[self.state.species]["GridSize"]
            for i in range(len(names)):
                length = gridsize / self.state.pos_units
                center_position = current_position + length / 2
                chr_label_positions.append(center_position)
                current_position += length
                chr_bounds.append(current_position)
                chr_labels.append(names[i])
            self.axes.set_xticks(chr_label_positions, minor=False)
            self.axes.set_xticklabels(chr_labels, minor=False)
            self.axes.set_xticks(chr_bounds, minor=True)

            self.axes.set_yticks(chr_label_positions, minor=False)
            self.axes.set_yticklabels(chr_labels, minor=False)
            self.axes.set_yticks(chr_bounds, minor=True)

            self.axes.tick_params(which="major", length=0)
            self.axes.tick_params(which="minor", length=10)

            self.axes.grid(visible=True, which="minor")

            self.state.x_axislabel = "Chr: " + self.state.x_att.label

            self.state.y_axislabel = "Chr: " + self.state.y_att.label

        self.axes.figure.canvas.draw_idle()

    def apply_roi(self, roi, override_mode=None):

        # Force redraw to get rid of ROI. We do this because applying the
        # subset state below might end up not having an effect on the viewer,
        # for example there may not be any layers, or the active subset may not
        # be one of the layers. So we just explicitly redraw here to make sure
        # a redraw will happen after this method is called.
        self.redraw()

        if len(self.layers) == 0:
            return

        x_date = "datetime" in self.state.x_kinds
        y_date = "datetime" in self.state.y_kinds

        if x_date or y_date:
            roi = roi.transformed(
                xfunc=mpl_to_datetime64 if x_date else None,
                yfunc=mpl_to_datetime64 if y_date else None,
            )

        use_transform = self.state.plot_mode != "rectilinear"
        subset_state = roi_to_subset_state(
            roi,
            x_att=self.state.x_att,
            x_categories=self.state.x_categories,
            y_att=self.state.y_att,
            y_categories=self.state.y_categories,
            use_pretransform=use_transform,
        )
        if use_transform:
            subset_state.pretransform = ProjectionMplTransform(
                self.state.plot_mode,
                self.axes.get_xlim(),
                self.axes.get_ylim(),
                self.axes.get_xscale(),
                self.axes.get_yscale(),
            )

        if self.state.lod_att is not None and self.state.lod_thresh is not None:
            lod_roi = RangeSubsetState(
                self.state.lod_thresh, 999, att=self.state.lod_att
            )  # hi=999 is a bad hack.
            # We could set this to the max value of self.state.lod_att instead
            subset_state = subset_state & lod_roi

        self.apply_subset_state(subset_state, override_mode=override_mode)
