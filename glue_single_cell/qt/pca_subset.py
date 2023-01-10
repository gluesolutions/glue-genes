import os
import numpy as np
from qtpy import QtWidgets
from echo.qt import autoconnect_callbacks_to_qt


from glue.core.subset import MultiOrState
from glue.utils.qt import load_ui
from glue.core import Data, Hub, HubListener
from glue.core.message import (SubsetMessage,
                               SubsetCreateMessage,
                               SubsetUpdateMessage,
                               SubsetDeleteMessage,
                               )
from glue.core.exceptions import IncompatibleAttribute

from ..state import PCASubsetState
from ..anndata_factory import df_to_data
from qtpy.QtWidgets import QMessageBox

import scanpy as sc
from scipy.sparse import issparse
import time

__all__ = ['PCASubsetDialog','GeneSummaryListener']


def dialog(title, text, icon):
    if icon=='warn':
        icon = QMessageBox.Warning
    elif icon=='info':
        icon = QMessageBox.Information
    info = QMessageBox()
    info.setIcon(icon)
    info.setText(title)
    info.setInformativeText(text)
    info.setStandardButtons(info.Ok)
    result = info.exec_()
    return True


def do_calculation_over_gene_subset(data_with_Xarray, genesubset, calculation = 'Means'):
    """
    """
    adata = data_with_Xarray.Xdata
    raw = False
    try: 
        mask = genesubset.to_index_list()
    except IncompatibleAttribute:
        dialog('Failed', "Failed to generate a mask on the selected subset.", 'warn')
        return None
    if mask.sum() == 0:
        return None
    #Slicing the adata object is probably not the fastest thing we can do
    if calculation == 'PCA':
        adata_sel = adata[:, mask]  # This will fail if genesubset is not actually over genes
        try:
            adata_sel = adata_sel.to_memory()
        except ValueError:
            pass
        sc.pp.pca(adata_sel, n_comps=10)
        data_arr = adata_sel.obsm['X_pca']
    elif calculation == 'Module':
        adata_sel = adata[:, mask]  # This will fail if genesubset is not actually over genes
        gene_list = list(adata_sel.var_names)
        try:
            # TODO/FIXME: This could crash glue if a really large file is loaded into memory
            # We do not need all genes in the reference dataset. We could subsample down to a
            # "reasonable" number of genes, append our target genes and only load that into memory
            adata_temp = adata.to_memory() 
        except ValueError: # ValueError if data is already in memory
            adata_temp = adata
        
        try:
            sc.tl.score_genes(adata_temp, gene_list = gene_list)
            data_arr = np.expand_dims(adata_temp.obs['score'],axis=1)
        except ValueError:
            print("No genes found!")
            return None

    elif calculation == 'Means':
        if raw:
            adata_sel = adata.raw.X[: , mask]  # This will fail if genesubset is not actually over genes
            if data_with_Xarray.sparse == True:
                data_arr = np.expand_dims(adata_sel.mean(axis=1).A1,axis=1)  # Expand to make same dimensionality as PCA
            else:
                data_arr = np.expand_dims(adata_sel.mean(axis=1),axis=1) 

        else:
            adata_sel = adata.X[: , mask]
            if data_with_Xarray.sparse == True:
                data_arr = np.expand_dims(adata_sel.mean(axis=1).A1,axis=1)  # Expand to make same dimensionality as PCA
            else:
                data_arr = np.expand_dims(adata_sel.mean(axis=1),axis=1) 

    return data_arr

def apply_data_arr(target_dataset, data_arr, basename, key='PCA'):
    """
    This is a sort of clunky approach to doing this.

    We generate a Data object from the data_arr that was returned,
    and then add the non-coordinate components of this Data object
    to the target dataset. 
    """
    
    data = Data(**{f'{key}_{i}':k for i,k in enumerate(data_arr.T)},label=f'{basename}_{key}')
    for x in data.components:
        if x not in data.coordinate_components:
            new_comp_name = f'{basename}_{x}'
            target_dataset.add_component(data.get_component(f'{x}'), new_comp_name)
    return new_comp_name #This just gets the last one, which is not quite correct for PCA


class GeneSummaryListener(HubListener):
    """
    A Listener to keep the new components in target_dataset
    up-to-date with any changes in the genesubset. 
    
    SubsetMessage define `subset` and `attribute` (for update?) 
    """
    def __init__(self, genesubset, basename, key, data_with_Xarray=None):
        self.genesubset = genesubset
        self.basename = basename
        self.key = key
        if data_with_Xarray is not None:
            self.set_circular_refs(data_with_Xarray)

    def set_circular_refs(self, data_with_Xarray):
        self.data_with_Xarray = data_with_Xarray
        self.hub = data_with_Xarray.hub
        self.target_dataset = self.data_with_Xarray.meta['obs_data']

    def register_to_hub(self, hub=None):
        if hub is not None:
            self.hub = hub
        if self.hub is None:
            self.hub = self.data_with_Xarray.hub
        #hub.subscribe(self, SubsetCreateMessage,
        #              handler=self.update_subset)
        self.hub.subscribe(self, SubsetUpdateMessage,
                      handler=self.update_subset)
        self.hub.subscribe(self, SubsetDeleteMessage,
                      handler=self.delete_subset)

    def __gluestate__(self, context):
        
        return dict(genesubset = context.id(self.genesubset),
                    basename = context.do(self.basename),
                    key = context.do(self.key),
                    data_with_Xarray = context.id(self.data_with_Xarray)
                )
    
    @classmethod
    def __setgluestate__(cls, rec, context):
        #target_dataset = context.object(rec['data_with_Xarray'])
        result =  cls(genesubset=context.object(rec['genesubset']),
                   basename=context.object(rec['basename']),
                   key=context.object(rec['key']),
                   )
        yield result
        data_with_Xarray = context.object(rec['data_with_Xarray'])
        result.set_circular_refs(data_with_Xarray)
        #result.register_to_hub()

    def update_subset(self, message):
        """
        if the subset is the one we care about
        then we rerun the calculation.

        TODO: If the subset is the same subset and has just been
        renamed then we need to update the component name
        """
        subset = message.subset
        if subset == self.genesubset:
            #if subset.attributes == self.genesubset_attributes:
            new_data = do_calculation_over_gene_subset(self.data_with_Xarray, self.genesubset, calculation = self.key)

            if new_data is not None:
                mapping = {f'{self.basename}_{self.key}_{i}':k for i,k in enumerate(new_data.T)}
                for x in self.target_dataset.components:  # This is to get the right component ids
                    xstr = f'{x.label}'
                    if xstr in mapping.keys():
                        mapping[x] = mapping.pop(xstr)
                self.target_dataset.update_components(mapping)
                

    def delete_subset(self, message):
        """
        Remove the attributes from target_dataset
        """
        pass
        
    def receive_message(self, message):
        pass

class PCASubsetDialog(QtWidgets.QDialog):

    def __init__(self, collect, default=None, parent=None):

        super(PCASubsetDialog, self).__init__(parent=parent)

        self.state = PCASubsetState(collect)

        self.ui = load_ui('pca_subset.ui', self,
                          directory=os.path.dirname(__file__))
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
            if target_dataset.meta['Xdata'] == data.uuid:
                data_with_Xarray = data
        
        for subset in self.state.genesubset.subsets:
            if subset.data == data_with_Xarray.meta['var_data']:  #  Find the subset on the genes, assuming we are adding to cell data
                genesubset = subset
                genesubset_attributes = subset.attributes
        if not genesubset:
            print(f"Selected subset {self.state.genesubset.label} does not seem to define genes in for {self.state.data.label}")

        basename = genesubset.label
        if self.state.do_means:
            key = 'Means'
        elif self.state.do_pca:
            key = 'PCA'
        elif self.state.do_module:
            key = "Module"
        data_arr = do_calculation_over_gene_subset(data_with_Xarray, genesubset, calculation = key)

        if data_arr is not None:
        
            new_comp_name = apply_data_arr(target_dataset, data_arr, basename, key=key)
            gene_summary_listener = GeneSummaryListener(genesubset, basename, key, data_with_Xarray)
            gene_summary_listener.register_to_hub()
            data_with_Xarray.listeners.append(gene_summary_listener)

        confirm = dialog('Adding a new component',
                        f'The component:\n'
                        f'{new_comp_name}\n'
                        f'has been added to:\n'
                        f'{target_dataset.label}\n'   
                        f'and will be automatically updated when {genesubset.label} is changed.',
                        'info')

    @classmethod
    def summarize(cls, collect, default=None, parent=None):
        self = cls(collect, parent=parent, default=default)
        value = self.exec_()

        if value == QtWidgets.QDialog.Accepted:
            self._apply()
