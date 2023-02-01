import os

from echo import keep_in_sync
from echo.qt import autoconnect_callbacks_to_qt
from glue.core.state_objects import CallbackProperty, State
from glue.utils import nonpartial
from glue.utils.decorators import avoid_circular
from glue.utils.qt import fix_tab_widget_fontsize, load_ui
from qtpy import QtWidgets

__all__ = ["SliderState", "SliderLabelWidget", "QTLOptionsWidget"]


class SliderState(State):

    label = CallbackProperty()
    slider_label = CallbackProperty()
    slider_pos = CallbackProperty()


class SliderLabelWidget(QtWidgets.QWidget):
    def __init__(self, label="", lo=0, hi=1, parent=None):
        super().__init__(parent=parent)
        self.state = SliderState()
        self.state.label = label
        self.lo = lo
        self.state.slider_pos = lo
        self.label_fmt = "{:.2f}"

        self.scale = (100) / (hi - lo)  # Slider has 100 increments

        self.ui = load_ui(
            "slider_with_label_widget.ui", self, directory=os.path.dirname(__file__)
        )
        connect_kwargs = {"slider_pos": dict(value_range=(lo, hi))}
        self._connections = autoconnect_callbacks_to_qt(
            self.state, self.ui, connect_kwargs
        )

        self.value_slider_pos.valueChanged.connect(
            nonpartial(self.set_label_from_slider)
        )

        # self.valuetext_slider_label.setMinimumWidth(80)
        self.state.slider_label = self.label_fmt.format(self.value_slider_pos.value())
        self.valuetext_slider_label.editingFinished.connect(
            nonpartial(self.set_slider_from_label)
        )

        self.set_label_from_slider()

    @avoid_circular
    def set_label_from_slider(self):
        value = self.state.slider_pos
        self.state.slider_label = self.label_fmt.format(value)

    @avoid_circular
    def set_slider_from_label(self):
        # Ignore recursive calls - we do this rather than ignore_callback
        # below when setting slider_label, otherwise we might be stopping other
        # subscribers to that event from being correctly updated
        if getattr(self, "_in_set_slider_from_label", False):
            return
        else:
            self._in_set_slider_from_label = True

        text = self.valuetext_slider_label.text()
        value = int((float(text) - self.lo) * self.scale)
        self.value_slider_pos.setValue(value)

        self._in_set_slider_from_label = False


class QTLOptionsWidget(QtWidgets.QWidget):
    def __init__(self, viewer_state, session, parent=None):

        super().__init__(parent=parent)

        self._slider = None  # This should only actually be a single entry, not a list
        self.ui = load_ui(
            "options_widget.ui", self, directory=os.path.dirname(__file__)
        )

        fix_tab_widget_fontsize(self.ui.tab_widget)

        self._connections = autoconnect_callbacks_to_qt(viewer_state, self.ui)

        connect_kwargs = {"alpha": dict(value_range=(0, 1))}
        self._connections_legend = autoconnect_callbacks_to_qt(
            viewer_state.legend, self.ui.legend_editor.ui, connect_kwargs
        )
        self.layout = self.ui.layout_slices
        self.layout.setSpacing(4)
        self.layout.setContentsMargins(0, 3, 0, 3)

        self.viewer_state = viewer_state

        self.viewer_state.add_callback("lod_att", self.update_slider_widget)
        self.update_slider_widget()

    def _clear(self):

        if self._slider is not None:
            # slider = self.layout.takeAt(0)
            self._slider.close()

        self._slider = None

    def update_slider_widget(self, *args):
        if self.viewer_state.lod_att is None:
            return

        self._clear()

        lod_slider = SliderLabelWidget(
            label="Thresh", lo=self.viewer_state.lod_min, hi=self.viewer_state.lod_max
        )
        self._slider = lod_slider
        self.slider_state = self._slider.state
        self._sync_lod_thresh = keep_in_sync(
            self.slider_state, "slider_pos", self.viewer_state, "lod_thresh"
        )
        self.layout.addWidget(self._slider)


if __name__ == "__main__":

    from glue.utils.qt import get_qapp

    app = get_qapp()

    widget = QTLOptionsWidget()
    # widget = SliderLabelWidget(label='LOD Threshold', lo=0.12, hi=13.54)
    # widget.show()

    app.exec_()
