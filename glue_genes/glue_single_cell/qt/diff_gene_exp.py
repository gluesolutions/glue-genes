"""
A menubar plugin to compute differential gene expression for two subsets.

This plugin is normally invoked with a GUI but the core logic can be
invoked as `get_gene_diff_exp`.

TODO: If the subset is over 1 or fewer datapoints then scanpy will crash, so we
        should check for this and disallow this subset.
TODO: Probably the logic here does not all work if subset2 is None/Rest

"""
import os
import numpy as np
import scanpy as sc
import pandas as pd
from echo.qt import autoconnect_callbacks_to_qt

from glue_qt.utils import load_ui
from glue.core.component_id import ComponentID
from glue.core.message import SubsetDeleteMessage, SubsetUpdateMessage
from glue.core import HubListener

from qtpy import QtWidgets

from glue_genes.glue_single_cell.component import SyncComponent
from ..state import DiffGeneExpState
from .summarize_gene_subset import dialog

__all__ = ["get_gene_diff_exp", "DiffGeneExpDialog", "DifferentialGeneExpressionListener"]


def get_gene_diff_exp(subset1, subset2, data):
    """
    Get differential gene expression for two subsets over a dataset

    Uses scanpy `rank_genes_groups` with just one group (subset1)
    and one reference set (subset2) and method=wilcoxon.

    Parameters
    ----------
    subset1 : :class:`~glue.core.subset.Subset`
        The highest ranked genes are expressed more in this subset.
    subset2 : :class:`~glue.core.subset.Subset`
        The reference subset. The lowest ranked genes are expressed more in this subset.
    data : :class:`~.DataAnnData`
        The gene expression (X) matrix connecting genes and cells.
    """

    adata = data.Xdata
    obsdata = data.meta["obs_data"]
    if subset2 is not None:
        m1 = obsdata.get_mask(subset1.subset_state)
        m2 = obsdata.get_mask(subset2.subset_state)
        # Cannot do rank_gene_groups with less than 2 observations
        if (np.sum(m1) < 2) or (np.sum(m2) < 2):
            raise ValueError("Failed to compute differential gene expression. "
                             "Cannot do differential gene expression with less than 2 observations in each subset. "
                             f"Subset 1 has {np.sum(m1)} observations, Subset 2 has {np.sum(m2)} observations.")
        conditions = [m1, m2]
        choices = ["1", "2"]

        adata.obs["glue_subsets"] = np.select(conditions, choices, default="0")
    else:  # Here we compare against the full dataset
        m1 = obsdata.get_mask(subset1.subset_state)
        # Cannot do rank_gene_groups with less than 2 observations
        if (np.sum(m1) < 2):
            raise ValueError("Failed to compute differential gene expression. "
                             "Cannot do differential gene expression with less than 2 observations in each subset. "
                             f"Subset 1 has {np.sum(m1)} observations.")
        conditions = [m1]
        choices = ["1"]
        adata.obs["glue_subsets"] = np.select(conditions, choices, default="2")

    adata_selected = adata[adata.obs["glue_subsets"] != "0", :]
    try:
        adata_selected = (
            adata_selected.to_memory()
        )  # We should check that this is not going to be too large
    except ValueError:
        pass
    sc.tl.rank_genes_groups(
        adata_selected, "glue_subsets", groups=["1"], reference="2", method="wilcoxon"
    )
    df = sc.get.rank_genes_groups_df(adata_selected, None)

    # Note that adata and adata_selected still have the same number of var_names
    # This sorts the results of rank genes based on the original gene order
    dummy = pd.Series(adata.var_names, name='names').to_frame()
    sorted_df = pd.merge(dummy, df, on='names', how='left')

    return sorted_df


class DiffGeneExpDialog(QtWidgets.QDialog):
    def __init__(self, collect, default=None, parent=None):
        super(DiffGeneExpDialog, self).__init__(parent=parent)

        self.state = DiffGeneExpState(collect)

        self.ui = load_ui("diff_gene_exp.ui", self, directory=os.path.dirname(__file__))
        self._connections = autoconnect_callbacks_to_qt(self.state, self.ui)

        self._collect = collect

        if default is not None:
            self.state.data = default

        self.ui.button_ok.clicked.connect(self.accept)
        self.ui.button_cancel.clicked.connect(self.reject)

    def _apply(self):
        """
        Calculate differential gene expression between two selected subsets
        or one subset versus all the rest of the data.

        This assumes that the subsets have been defined over on the obs array

        Note that this copies the relevant subsets into memory as otherwise
        this won't work on anndata in disk-backed mode:

        https://github.com/theislab/scanpy/issues/2147

        (the above is technically for a different scanpy function, but the same problem occurs for rank_genes_groups)
        """
        try:
            df = get_gene_diff_exp(
                self.state.subset1,
                self.state.subset2,
                self.state.data
            )
        except ValueError as e:
            confirm = dialog(  # noqa: F841
                "Failed to compute differential gene expression",
                f"{e}",
                "warn")

        label1 = self.state.subset1.label
        if self.state.subset2 is None:
            label2 = "Rest"
            subsets = [self.state.subset1]
            message = f"{label1}"
        else:
            label2 = self.state.subset2.label
            subsets = [self.state.subset1, self.state.subset2]
            message = f"{label1} or {label2}"

        if df is not None:
            vardata = self.state.data.meta["var_data"]
            component_name = f"z-scores for {label1} vs {label2} (sync)"
            cid1 = ComponentID(component_name)
            comp1 = SyncComponent(df['scores'].values, subsets=subsets)
            _ = vardata.add_component(comp1, cid1)
            component_name = f"Adj p-vals for {label1} vs {label2} (sync)"
            cid2 = ComponentID(component_name)
            comp2 = SyncComponent(df['pvals_adj'].values, subsets=subsets)
            _ = vardata.add_component(comp2, cid2)

            dge_listener = DifferentialGeneExpressionListener(self.state.data,
                                                              self.state.subset1,
                                                              self.state.subset2, comps=[cid1, cid2])
            dge_listener.register_to_hub()
            self.state.data.listeners.append(dge_listener)

            confirm = dialog(  # noqa: F841
                "Adding a new component",
                f"The components:\n"
                f"{comp1} and {comp2}\n"
                f"have been added to:\n"
                f"{vardata.label}\n"
                f"and will be automatically updated when {message} is changed.",
                "info")

    @classmethod
    def calculate_deg(cls, collect, default=None, parent=None):
        self = cls(collect, parent=parent, default=default)
        value = self.exec_()

        if value == QtWidgets.QDialog.Accepted:
            self._apply()


class DifferentialGeneExpressionListener(HubListener):
    """
    A Listener that updates the differential gene expression SynComponent(s)
    when the underlying subsets are changed.

    Parameters
    ----------
    data_with_Xarray : :class:`~.DataAnnData`
        The data object containing the var_data component to be updated
    subset1 : :class:`~glue.core.subset.Subset`
        The subset that is the numerator in the differential gene expression
    subset2 : :class:`~glue.core.subset.Subset`
        The subset that is the denominator in the differential gene expression
    comps : list of :class:`~glue.core.component_id.ComponentID`
        The components to be kept up-to-date.
    """

    def __init__(self, data_with_Xarray, subset1, subset2, comps=[]):
        super().__init__()
        self.data = data_with_Xarray
        self.subset1 = subset1
        self.subset2 = subset2
        self.comps = comps
        self.problem_subset = False
        if data_with_Xarray is not None:
            self.set_circular_refs(data_with_Xarray)

    def set_circular_refs(self, data_with_Xarray):
        self.data_with_Xarray = data_with_Xarray
        self.hub = data_with_Xarray.hub
        self.target_dataset = self.data_with_Xarray.meta["var_data"]

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
        if the subset is the one we care about
        then we rerun the calculation.

        TODO: If the subset is the same subset and has just been
        renamed then we need to update the component name

        TODO: If the subset is no longer over a valid set of attributes
              we should... do what?
        """
        subset = message.subset
        if subset in self.subset2.subsets:
            subset2valid = True
        else:
            subset2valid = False

        if (subset in self.subset1.subsets) or (subset2valid):
            try:
                new_df = get_gene_diff_exp(self.subset1, self.subset2, self.data_with_Xarray)
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
            mapping[self.comps[0]] = new_df['scores'].values
            mapping[self.comps[1]] = new_df['pvals_adj'].values

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
            subset1=context.id(self.subset1),
            subset2=context.id(self.subset2),
            data_with_Xarray=context.id(self.data),
            comps=context.id(self.comps),
        )

    @classmethod
    def __setgluestate__(cls, rec, context):
        result = cls(
            subset1=context.object(rec["subset1"]),
            subset2=context.object(rec["subset2"]),
            comps=context.object(rec["comps"]),
        )
        yield result
        data_with_Xarray = context.object(rec["data_with_Xarray"])
        result.set_circular_refs(data_with_Xarray)
