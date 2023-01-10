import os
from qtpy import QtWidgets
from echo.qt import autoconnect_callbacks_to_qt
import pandas as pd

from glue.utils.qt import load_ui

from ..anndata_factory import df_to_data
from ..state import EnrichpyState
from .pca_subset import dialog

import enrichrpy.enrichr

__all__ = ['EnrichpyDialog']

def convert_genes_to_list(row):
    arr = row.split(';')
    l = [x for x in arr]
    return l

def convert_to_mouse_ids(row):
    """
    Mouse Entrez IDs are like Sox17, while enrichr returns SOX17
    This function will restore these gene IDs.
    
    Ideally we should figure out when we need to capitalize the
    gene list and only restore if we need to. 
    """
    return [x.capitalize() for x in row]

class EnrichpyDialog(QtWidgets.QDialog):

    def __init__(self, collect, default=None, parent=None):

        super(EnrichpyDialog, self).__init__(parent=parent)

        self.state = EnrichpyState(collect)

        self.ui = load_ui('enrichr.ui', self,
                          directory=os.path.dirname(__file__))
        self._connections = autoconnect_callbacks_to_qt(self.state, self.ui)

        self._collect = collect

        if default is not None:
            self.state.data = default

        self.ui.button_ok.clicked.connect(self.accept)
        self.ui.button_cancel.clicked.connect(self.reject)

    def _apply(self, do_dialog=True):
        """
        """
        if self.state.subset is not None:
            for subset in self.state.subset.subsets:
                if subset.data == self.state.data:
                    gene_subset = subset
        gene_list = gene_subset[self.state.gene_att]
        gene_list_upper = [x for x in gene_list] #[x.upper() for x in gene_list] Maybe just for humans?
        output_df = enrichrpy.enrichr.get_pathway_enrichment(gene_list_upper, gene_set_library=self.state.gene_set)
        output_df['Overlapping genes'] = output_df['Overlapping genes'].apply(lambda x: '|'.join(x))
        new_name = f'{self.state.gene_set} for {self.state.subset.label}'
        results_data = df_to_data(output_df,label=new_name)
        # Output.results has Pandas Object types which don't play nice
        # with np.load in session files (unless we change to have allow_pickle=True)
        # so we just force them to be strings in glue
        for component in results_data.components:
            try:
                if results_data[component].dtype == 'O':
                    results_data[component] = results_data[component].astype('str')
            except AttributeError:
                pass
        
        self._collect.append(results_data)
        
        if do_dialog:
            confirm = dialog('New dataset created',
                    f'The Dataset:\n'
                    f'{new_name}\n'
                    f'has been created.',
                    'info')

        # Theoretically we could then join on the Genes column somehow...
        # But the Genes column is Gene1;Gene2; etc. so it is not trivial and probably it is not interesting to do this.
        
        
    @classmethod
    def enrich(cls, collect, default=None, parent=None):
        self = cls(collect, parent=parent, default=default)
        value = self.exec_()

        if value == QtWidgets.QDialog.Accepted:
            self._apply()
