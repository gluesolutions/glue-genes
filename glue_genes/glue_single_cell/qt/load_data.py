import os

import anndata as ad
import h5py
from glue.utils.qt import load_ui
from psutil import virtual_memory
from psutil._common import bytes2human
from qtpy.QtCore import Qt
from qtpy.QtWidgets import QDialog, QListWidgetItem

__all__ = ["LoadDataDialog"]

MEM_WARNING = """The file you are loading is comparable in size \
to available RAM on your system. You should open the file in backed mode.
"""

NO_MEM_WARNING = """Estimated memory usage is much less than system memory."""

LARGE_X_WARNING = """The file you are loading contains a very large number \
of measurements (n_obs x n_vars > 1e8) which may result in poor performance for \
certain interactive features, including calculating differential gene \
expression and calculating the summary over a gene subset. Glue can subsample \
random cells to improve performance or you can cancel and reduce the filesize \
outside of glue. If you can tolerate the above features being slow, you can \
also just proceed to load the file as-is.
"""

NO_LARGE_X_WARNING = """The number of measurements (n_obs x n_vars) in this \
dataset is a good size for interactive features in glue. There is no need \
to subsample this dataset."""

CANCEL_HELP = """Choose cancel to do your own preprocessing of the data to \
reduce its size. We recommend the \
<a href=\"https://github.com/tanaylab/metacells\">Metacells package</a> \
which will reduce the dataset down to a small number of metacells carefully \
chosen to represent the full dataset. Running metacells properly is an \
interactive process and may require significant computational resources. See \
<a href=\"https://metacells.readthedocs.io/en/latest/Metacells_Vignette.html\">\
this tutorial</a> for additional information."""


NO_CANCEL_HELP = """"""

# We do not have the exact compression of a h5ad file
# This is just a conservative guess from some sample files
COMPRESSION_FACTOR = 5


def get_system_memory():
    mem = virtual_memory()
    total_mem = mem.total
    human_mem = bytes2human(total_mem)
    return (total_mem, human_mem)


class LoadDataDialog(QDialog):
    """
    A dialog to allow the user to choose how to load an :class:`~.DataAnnData` object.

    Currently we return parameters directly from this dialog rather than
    the ideal approach which is to encapsulate these parameters in a
    custom State class.

    Parameters
    ----------
    filename : str
                Path to data file to load
    """

    def __init__(self, filename=None, parent=None):

        super(LoadDataDialog, self).__init__(parent=parent)

        self.subsample = False
        self.try_backed = (
            None  # This is the default for scanpy.read to read into memory
        )
        self.skip_components = []
        self.subsample_factor = 1
        trial_read = ad.read(filename, backed="r")
        nobs = trial_read.n_obs
        nvars = trial_read.n_vars
        filesize = os.path.getsize(filename)

        # Some AnnData files in the wild have groups structured in different ways
        try:
            with h5py.File(filename) as h5f:
                compression = h5f["X/data"].compression
        except (AttributeError, KeyError):
            try:
                with h5py.File(filename) as h5f:
                    compression = h5f["X"].compression
            except (AttributeError, KeyError):
                compression = None

        if compression:
            filesize = filesize * COMPRESSION_FACTOR

        human_file_size = bytes2human(filesize)

        self.ui = load_ui(
            "load_data.ui", parent=self, directory=os.path.dirname(__file__)
        )
        # self._connections = autoconnect_callbacks_to_qt(self.state, self.ui)

        self.ui.label_n_obs.setText(str(nobs))
        self.ui.label_n_vars.setText(str(nvars))
        self.ui.label_memory_usage.setText(str(human_file_size))

        system_mem, human_system_memory = get_system_memory()
        self.ui.label_system_memory.setText(str(human_system_memory))

        if filesize > system_mem / 2:  # Make sure we have plenty of RAM
            do_memory_warning = True
        else:
            do_memory_warning = False

        if do_memory_warning:
            self.ui.label_memory_warning.setText(MEM_WARNING)
            self.ui.label_memory_warning.setStyleSheet("color: red")
            self.ui.checkbox_backed.checked = True
        else:
            self.ui.label_memory_warning.setText(NO_MEM_WARNING)
            self.ui.label_memory_warning.setStyleSheet("color: black")

        if nobs * nvars > 1e8:
            do_large_x_warning = True
            minimum_points = 10_000 / nobs  # We always keep at least 10_000 obs
            self.subsample_factor = max(1e8 / (nobs * nvars), minimum_points)
        else:
            do_large_x_warning = False
            self.subsample_factor = 0.1  # If the data is not that large, just do 10%

        if do_large_x_warning:
            self.ui.label_largex_warning.setText(LARGE_X_WARNING)
            self.ui.label_largex_warning.setStyleSheet("color: red")
        else:
            self.ui.label_largex_warning.setText(NO_LARGE_X_WARNING)
            self.ui.label_largex_warning.setStyleSheet("color: black")

        if do_large_x_warning or do_memory_warning:
            self.ui.button_subsample.setDefault(True)
            self.ui.label_cancel_help.setText(CANCEL_HELP)
        else:
            self.ui.button_ok.setDefault(True)
            self.ui.label_cancel_help.setText(NO_CANCEL_HELP)

        obss = trial_read.obs_keys()
        self.populate_list("obs", obss)

        varss = trial_read.var_keys()
        self.populate_list("var", varss)

        obsms = trial_read.obsm_keys()
        self.populate_list("obsm", obsms)

        varms = trial_read.varm_keys()
        self.populate_list("varm", varms)

        self.ui.button_cancel.clicked.connect(self.reject)
        self.ui.button_ok.clicked.connect(self.accept)
        self.ui.button_subsample.clicked.connect(self.do_subsample)

        self.ui.button_select_none.clicked.connect(self.select_none)
        self.ui.button_select_all.clicked.connect(self.select_all)

        self.ui.list_component.itemChanged.connect(self._on_check_change)

        trial_read.file.close()
        del trial_read

    def populate_list(self, data_type, data_set):
        item = QListWidgetItem(data_type)
        item.setFlags(item.flags() & ~Qt.ItemIsSelectable)
        item.setForeground(Qt.gray)
        self.ui.list_component.addItem(item)
        for data_obj in data_set:
            item = QListWidgetItem(data_obj)
            item.setCheckState(Qt.Checked)
            self.ui.list_component.addItem(item)

    def _on_check_change(self, *event):

        any_checked = False

        for idx in range(self.ui.list_component.count()):
            item = self.ui.list_component.item(idx)
            if item.checkState() == Qt.Checked:
                any_checked = True
                break

        self.button_ok.setEnabled(any_checked)

    def select_none(self, *event):
        self._set_all_checked(False)

    def select_all(self, *event):
        self._set_all_checked(True)

    def _set_all_checked(self, check_state):
        for idx in range(self.ui.list_component.count()):
            item = self.ui.list_component.item(idx)
            item.setCheckState(Qt.Checked if check_state else Qt.Unchecked)

    def set_components(self):
        skip_components = []
        for idx in range(self.ui.list_component.count()):
            item = self.ui.list_component.item(idx)
            if (item.checkState() != Qt.Checked) and (item.foreground() != Qt.gray):
                skip_components.append(item.text())
        self.skip_components = skip_components

    def accept(self):
        self.set_components()
        if self.ui.checkbox_backed.isChecked():
            self.try_backed = "r"
        super(LoadDataDialog, self).accept()

    def do_subsample(self):
        self.set_components()
        self.subsample = True
        if self.ui.checkbox_backed.isChecked():
            self.try_backed = "r"
        super(LoadDataDialog, self).accept()


if __name__ == "__main__":

    from glue.utils.qt import get_qapp

    app = get_qapp()

    dialog = LoadDataDialog(filename="test2.h5ad")
    if dialog.exec_():
        print(dialog.components)
        print(dialog.subsample)
    else:
        print("Rejected")
