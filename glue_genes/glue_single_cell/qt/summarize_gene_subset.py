"""
Calculate a summary statistic for cells based on genesubset

This menubar plugin in normally invoked through a dialog
which not only does the initial calculation and stores the
summary is a new component in the target dataset, but sets up
a Listener to keep that new component up-to-date when changes
are made to the gene subset.
"""

import os

import numpy as np
import scanpy as sc
from echo.qt import autoconnect_callbacks_to_qt
from glue.core import HubListener
from glue.core.exceptions import IncompatibleAttribute
from glue.core.message import SubsetDeleteMessage, SubsetUpdateMessage
from glue_qt.utils import load_ui
from qtpy import QtWidgets
from qtpy.QtWidgets import QMessageBox
from glue_genes.glue_single_cell.component import SyncComponent
from glue.core.component_id import ComponentID

from ..state import SummarizeGeneSubsetState

__all__ = ["do_calculation_over_gene_subset", "SummarizeGeneSubsetDialog", "GeneSummaryListener"]


def dialog(title, text, icon):
    if icon == "warn":
        icon = QMessageBox.Warning
    elif icon == "info":
        icon = QMessageBox.Information
    info = QMessageBox()
    info.setIcon(icon)
    info.setText(title)
    info.setInformativeText(text)
    info.setStandardButtons(info.Ok)
    result = info.exec_()  # noqa: F841
    return True


def do_calculation_over_gene_subset(data_with_Xarray, genesubset, calculation="Means"):
    """
    Calculate a summary statistic for cells based on genesubset

    Currently this code is a bit of a mess because the slicing we are
    doing here implicitly? uses the new todense() cast on the X array
    which causes all this fragile code originally deal with sparse
    arrays to fail. Of course, in theory we would like to keep the
    benefit of sparse arrays here. I think the call to get_data()
    is inside the slice opperation? Or possibly in the to_index_list()
    call? I am not sure. It might even be in our call to
    data_with_Xarray.Xdata. But... I don't think so.



    Parameters
    ----------
    data_with_Xarray : :class:`~.DataAnnData`
        The expression matrix (X) used for the calculation
    genesubset : :class:`glue.core.subset.Subset`
        The subset of genes to use for the calculation
    calculation : str
        The type of calculation to perform: [PCA, Module, Means]

    """
    adata = data_with_Xarray.Xdata
    raw = False
    try:
        mask = genesubset.to_index_list()
    except IncompatibleAttribute:
        dialog("Failed", "Failed to generate a mask on the selected subset.", "warn")
        return None
    if mask.sum() == 0:
        return None
    # Slicing the adata object is probably not the fastest thing we can do
    if calculation == "PCA":
        adata_sel = adata[:, mask]  # This will fail if genesubset is not actually over genes
        try:
            adata_sel = adata_sel.to_memory()
        except ValueError:
            pass
        sc.pp.pca(adata_sel, n_comps=10)
        data_arr = adata_sel.obsm["X_pca"]
    elif calculation == "Module":
        adata_sel = adata[:, mask]  # This will fail if genesubset is not actually over genes
        gene_list = list(adata_sel.var_names)
        try:
            # TODO/FIXME: This could crash glue if a really large file is loaded into memory
            # We do not need all genes in the reference dataset. We could subsample down to a
            # "reasonable" number of genes, append our target genes and only load that into memory
            adata_temp = adata.to_memory()
        except ValueError:  # ValueError if data is already in memory
            adata_temp = adata

        try:
            sc.tl.score_genes(adata_temp, gene_list=gene_list)
            data_arr = np.expand_dims(adata_temp.obs["score"], axis=1)
        except ValueError:
            print("No genes found!")
            return None

    elif calculation == "Means":
        if raw:
            adata_sel = adata.raw.X[:, mask]  # This will fail if genesubset is not actually over genes
            if data_with_Xarray.sparse is True:
                data_arr = adata_sel.mean(axis=1)#.A1,)
            else:
                data_arr = adata_sel.mean(axis=1)

        else:
            adata_sel = adata.X[:, mask]
            if data_with_Xarray.sparse is True:
                data_arr = adata_sel.mean(axis=1)#.A1
            else:
                data_arr = adata_sel.mean(axis=1)

    return np.squeeze(np.asarray(data_arr))


def apply_data_arr(target_dataset, data_arr, basename, subset, key="Means"):
    """
    Add appropriately named SyncComponents from data_arr to target_dataset

    If we are doing PCA we will have multiple components to add

    Parameters
    ----------
    target_dataset : :class:`~.Data`
        The dataset to add the components to
    data_arr : :class:`numpy.ndarray`
        The array of data to add, this will be a 1D array for Means and Module
        and a 2D array for PCA
    basename : str
        The name of the subset
    subset : :class:`glue.core.subset.Subset`
        The subset that will be tracked by the SyncComponent
    key : str
        The type of calculation performed: [PCA, Module, Means]
    """
    cids = []
    if data_arr.ndim == 2:
        for i, k in enumerate(data_arr.T):
            component_name = f"{basename}_{key}_{i} (sync)"
            cid = ComponentID(component_name)
            comp = SyncComponent(k, subsets=[subset])
            _ = target_dataset.add_component(comp, cid)
            cids.append(cid)
    else:
        component_name = f"{basename}_{key} (sync)"
        cid = ComponentID(component_name)
        comp = SyncComponent(data_arr, subsets=[subset])
        _ = target_dataset.add_component(comp, cid)
        cids.append(cid)
    return cids


class GeneSummaryListener(HubListener):
    """
    A Listener to keep the new components in target_dataset
    up-to-date with any changes in the genesubset.

    Parameters
    ----------
    genesubset : Subset
        The subset over genes to watch
    basename : str
        The subset label
    key : str
        The kind of summary performed: [Means, PCA, Module]
    data_with_Xarray : :class:`~.DataAnnData`
        The expression matrix (X) used to reference the target_dataset
    comps : list
        A list of the components to keep in sync
    """

    def __init__(self, genesubset, basename, key, data_with_Xarray=None, comps=[]):
        self.genesubset = genesubset
        self.basename = basename
        self.key = key
        self.comps = comps
        self.problem_subset = False
        if data_with_Xarray is not None:
            self.set_circular_refs(data_with_Xarray)

    def set_circular_refs(self, data_with_Xarray):
        self.data_with_Xarray = data_with_Xarray
        self.hub = data_with_Xarray.hub
        self.target_dataset = self.data_with_Xarray.meta["obs_data"]

    def register_to_hub(self, hub=None):
        if hub is not None:
            self.hub = hub
        if self.hub is None:
            self.hub = self.data_with_Xarray.hub
        # hub.subscribe(self, SubsetCreateMessage,
        #              handler=self.update_subset)
        self.hub.subscribe(self, SubsetUpdateMessage, handler=self.update_subset)
        self.hub.subscribe(self, SubsetDeleteMessage, handler=self.delete_subset)

    def update_subset(self, message):
        """
        if the subset is the one we care about then we rerun the calculation.

        If we are just renaming the subset, then we just update
        the component names.

        TODO: If the subset is no longer over a valid set of attributes
              we should... what?
        """
        subset = message.subset
        if subset == self.genesubset:
            if message.attribute == "label":
                if len(self.comps) == 1:
                    self.comps[0].label = f"{subset.label}_{self.key} (sync)"
                else:
                    for i, c in enumerate(self.comps):
                        c.label = f"{subset.label}_{self.key}_{i} (sync)"
            # if subset.attributes == self.genesubset_attributes:
                return
            try:
                new_data = do_calculation_over_gene_subset(
                    self.data_with_Xarray, self.genesubset, calculation=self.key
                )
            except ValueError as e:
                if not self.problem_subset:
                    confirm = dialog(  # noqa: F841
                        "Failed to update differential gene expression",
                        f"{e}",
                        "warn")
                else:
                    pass
                self.problem_subset = True
                return
            self.problem_subset = False
            mapping = {}
            if new_data.ndim == 2:
                for c, k in zip(self.comps, new_data.T):
                    mapping[c] = k
            else:
                mapping[self.comps[0]] = new_data
            self.target_dataset.update_components(mapping)

    def delete_subset(self, message):
        """
        TODO: Remove the attributes from target_dataset
        """
        pass

    def receive_message(self, message):
        pass

    def __gluestate__(self, context):
        return dict(
            genesubset=context.id(self.genesubset),
            basename=context.do(self.basename),
            key=context.do(self.key),
            data_with_Xarray=context.id(self.data_with_Xarray),
            comps=context.id(self.comps),
        )

    @classmethod
    def __setgluestate__(cls, rec, context):
        # target_dataset = context.object(rec['data_with_Xarray'])
        result = cls(
            genesubset=context.object(rec["genesubset"]),
            basename=context.object(rec["basename"]),
            key=context.object(rec["key"]),
            comps=context.object(rec["comps"]),
        )
        yield result
        data_with_Xarray = context.object(rec["data_with_Xarray"])
        result.set_circular_refs(data_with_Xarray)
        # result.register_to_hub()


class SummarizeGeneSubsetDialog(QtWidgets.QDialog):
    def __init__(self, collect, default=None, parent=None):
        super().__init__(parent=parent)

        self.state = SummarizeGeneSubsetState(collect)

        self.ui = load_ui(
            "summarize_gene_subset.ui", self, directory=os.path.dirname(__file__)
        )
        self._connections = autoconnect_callbacks_to_qt(self.state, self.ui)

        self._collect = collect

        if default is not None:
            self.state.data = default

        self.ui.button_ok.clicked.connect(self.accept)
        self.ui.button_cancel.clicked.connect(self.reject)

    def _apply(self):
        """
        Apply the appropriate calculation over the X array for the
        gene susbset, add the result to the target dataset and
        add a GeneSummaryListener to keep that up-to-date.
        """
        genesubset = None

        target_dataset = self.state.data

        for data in self._collect:
            if target_dataset.meta["Xdata"] == data.uuid:
                data_with_Xarray = data

        for subset in self.state.genesubset.subsets:
            if (
                subset.data == data_with_Xarray.meta["var_data"]
            ):  # Find the subset on the genes, assuming we are adding to cell data
                genesubset = subset
        if not genesubset:
            print(
                f"Selected subset {self.state.genesubset.label} does not define genes in {self.state.data.label}"
            )

        basename = genesubset.label
        if self.state.do_means:
            key = "Means"
        elif self.state.do_pca:
            key = "PCA"
        elif self.state.do_module:
            key = "Module"
        data_arr = do_calculation_over_gene_subset(
            data_with_Xarray, genesubset, calculation=key
        )
        new_comps = []
        if data_arr is not None:
            new_comps = apply_data_arr(
                target_dataset, data_arr, basename, genesubset, key=key
            )
            gene_summary_listener = GeneSummaryListener(
                genesubset, basename, key, data_with_Xarray, new_comps
            )
            gene_summary_listener.register_to_hub()
            data_with_Xarray.listeners.append(gene_summary_listener)

        confirm = dialog(  # noqa: F841
            "Adding a new component",
            f"The component:\n"
            f"{new_comps}\n"
            f"has been added to:\n"
            f"{target_dataset.label}\n"
            f"and will be automatically updated when {genesubset.label} is changed.",
            "info",
        )

    @classmethod
    def summarize(cls, collect, default=None, parent=None):
        self = cls(collect, parent=parent, default=default)
        value = self.exec_()

        if value == QtWidgets.QDialog.Accepted:
            self._apply()
