import os

import numpy as np
from echo.qt import autoconnect_callbacks_to_qt, connect_value
from glue.utils.qt import fix_tab_widget_fontsize, load_ui
from qtpy import QtGui, QtWidgets
from qtpy.QtCore import Qt


class QTLLayerStyleEditor(QtWidgets.QWidget):
    def __init__(self, layer, parent=None):

        super(QTLLayerStyleEditor, self).__init__(parent=parent)

        self.ui = load_ui(
            "layer_style_editor.ui", self, directory=os.path.dirname(__file__)
        )

        connect_kwargs = {
            "alpha": dict(value_range=(0, 1)),
            "size_scaling": dict(value_range=(0.1, 10), log=True),
            "density_contrast": dict(value_range=(0, 1)),
        }
        self._connections = autoconnect_callbacks_to_qt(
            layer.state, self.ui, connect_kwargs
        )

        self._connection_dpi = connect_value(
            layer.state.viewer_state,
            "dpi",
            self.ui.value_dpi,
            value_range=(12, 144),
            log=True,
        )

        fix_tab_widget_fontsize(self.ui.tab_widget)

        self.layer_state = layer.state

        self.layer_state.add_callback("markers_visible", self._update_markers_visible)

        self.layer_state.add_callback("cmap_mode", self._update_cmap_mode)
        self.layer_state.add_callback("size_mode", self._update_size_mode)

        self.layer_state.add_callback("density_map", self._update_size_mode)
        self.layer_state.add_callback("density_map", self._update_warnings)
        self.layer_state.add_callback("density_map", self._update_checkboxes)

        self.layer_state.viewer_state.add_callback("x_att", self._update_checkboxes)
        self.layer_state.viewer_state.add_callback("y_att", self._update_checkboxes)

        self.layer_state.add_callback("layer", self._update_warnings)

        self._update_markers_visible()

        self._update_size_mode()
        self._update_cmap_mode()

        self._update_checkboxes()

        self._update_warnings()

    def _update_warnings(self, *args):

        if self.layer_state.layer is None:
            n_points = 0
        else:
            n_points = np.product(self.layer_state.layer.shape)

        warning = " (may be slow given data size)"

        for combo, threshold in [
            (self.ui.combosel_size_mode, 10000),
            (self.ui.combosel_cmap_mode, 50000),
        ]:

            if n_points > threshold and not self.layer_state.density_map:
                for item in range(combo.count()):
                    text = combo.itemText(item)
                    if text != "Fixed":
                        combo.setItemText(item, text + warning)
                        combo.setItemData(item, QtGui.QBrush(Qt.red), Qt.TextColorRole)
            else:
                for item in range(combo.count()):
                    text = combo.itemText(item)
                    if text != "Fixed":
                        if warning in text:
                            combo.setItemText(item, text.replace(warning, ""))
                            combo.setItemData(item, QtGui.QBrush(), Qt.TextColorRole)

    def _update_size_mode(self, size_mode=None):

        visible = (
            not self.layer_state.density_map
            and not self.layer_state.size_mode == "Fixed"
        )
        self.ui.label_size_attribute.setVisible(visible)
        self.ui.combosel_size_att.setVisible(visible)
        self.ui.label_size_limits.setVisible(visible)
        self.ui.valuetext_size_vmin.setVisible(visible)
        self.ui.valuetext_size_vmax.setVisible(visible)
        self.ui.button_flip_size.setVisible(visible)

        visible = (
            not self.layer_state.density_map and self.layer_state.size_mode == "Fixed"
        )
        self.ui.value_size.setVisible(visible)

        density = self.layer_state.density_map
        self.ui.value_dpi.setVisible(density)
        self.ui.label_dpi.setVisible(density)
        self.ui.label_stretch.setVisible(density)
        self.ui.combosel_stretch.setVisible(density)
        self.ui.value_density_contrast.setVisible(density)
        self.ui.label_contrast.setVisible(density)
        self.ui.combosel_size_mode.setVisible(not density)
        self.ui.value_size_scaling.setVisible(not density)
        self.ui.label_size_mode.setVisible(not density)
        self.ui.label_size_scaling.setVisible(not density)
        self.ui.label_fill.setVisible(not density)
        self.ui.bool_fill.setVisible(not density)

    def _update_markers_visible(self, *args):
        self.ui.combosel_size_mode.setEnabled(self.layer_state.markers_visible)
        self.ui.value_size.setEnabled(self.layer_state.markers_visible)
        self.ui.combosel_size_att.setEnabled(self.layer_state.markers_visible)
        self.ui.valuetext_size_vmin.setEnabled(self.layer_state.markers_visible)
        self.ui.valuetext_size_vmax.setEnabled(self.layer_state.markers_visible)
        self.ui.button_flip_size.setEnabled(self.layer_state.markers_visible)
        self.ui.value_size_scaling.setEnabled(self.layer_state.markers_visible)
        self.ui.value_dpi.setEnabled(self.layer_state.markers_visible)
        self.ui.combosel_stretch.setEnabled(self.layer_state.markers_visible)
        self.ui.label_size_scaling.setEnabled(self.layer_state.markers_visible)
        self.ui.combosel_points_mode.setEnabled(self.layer_state.markers_visible)
        self.ui.value_density_contrast.setEnabled(self.layer_state.markers_visible)

    def _update_checkboxes(self, *args):
        pass

    def _update_cmap_mode(self, cmap_mode=None):

        if self.layer_state.cmap_mode == "Fixed":
            self.ui.label_cmap_attribute.hide()
            self.ui.combosel_cmap_att.hide()
            self.ui.label_cmap_limits.hide()
            self.ui.valuetext_cmap_vmin.hide()
            self.ui.valuetext_cmap_vmax.hide()
            self.ui.button_flip_cmap.hide()
            self.ui.combodata_cmap.hide()
            self.ui.label_colormap.hide()
            self.ui.color_color.show()
        else:
            self.ui.label_cmap_attribute.show()
            self.ui.combosel_cmap_att.show()
            self.ui.label_cmap_limits.show()
            self.ui.valuetext_cmap_vmin.show()
            self.ui.valuetext_cmap_vmax.show()
            self.ui.button_flip_cmap.show()
            self.ui.combodata_cmap.show()
            self.ui.label_colormap.show()
            self.ui.color_color.hide()
