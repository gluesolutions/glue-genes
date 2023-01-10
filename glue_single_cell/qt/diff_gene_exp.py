import os
import numpy as np
from qtpy import QtWidgets
from echo.qt import autoconnect_callbacks_to_qt

from glue.core.subset import CategorySubsetState 
from glue.utils.qt import load_ui
from glue.core.link_helpers import JoinLink
from glue.core.data import Data

from ..state import DiffGeneExpState
from .pca_subset import dialog

import scanpy as sc

__all__ = ['DiffGeneExpDialog']


def get_gene_list_diff_exp(subset1, subset2, data, n_genes=50):

    adata = data.Xdata
    obsdata = data.meta['obs_data']
    conditions = [
        (obsdata.get_mask(subset1.subset_state)),
        (obsdata.get_mask(subset2.subset_state))]
    
    choices = ['1','2']
    
    adata.obs['glue_subsets'] = np.select(conditions, choices, default='0')
    adata_selected = adata[adata.obs['glue_subsets'] != '0', :]
    try:
        adata_selected = adata_selected.to_memory()  # We should check that this is not going to be too large
    except ValueError:
        pass
    sc.tl.rank_genes_groups(adata_selected, 'glue_subsets', groups=['1'], reference='2', method='wilcoxon')
    
    gene_list = [x[0] for x in adata_selected.uns['rank_genes_groups']['names']]
    scores = [x[0] for x in adata_selected.uns['rank_genes_groups']['scores']]    
    pvals_adj = [x[0] for x in adata_selected.uns['rank_genes_groups']['pvals_adj']]
    pvals = [x[0] for x in adata_selected.uns['rank_genes_groups']['pvals']]

    dge_data = Data(gene_names=gene_list, scores=scores, pvals=pvals, pvals_adj = pvals_adj)

    # For the starting subset we return just the top N genes
    gene_list = gene_list[0:n_genes]
    
    return gene_list, dge_data

class DiffGeneExpDialog(QtWidgets.QDialog):

    def __init__(self, collect, default=None, parent=None):

        super(DiffGeneExpDialog, self).__init__(parent=parent)

        self.state = DiffGeneExpState(collect)

        self.ui = load_ui('diff_gene_exp.ui', self,
                          directory=os.path.dirname(__file__))
        self._connections = autoconnect_callbacks_to_qt(self.state, self.ui)

        self._collect = collect

        if default is not None:
            self.state.data = default

        self.ui.button_ok.clicked.connect(self.accept)
        self.ui.button_cancel.clicked.connect(self.reject)

    def _apply(self):
        """
        Calculate differential gene expression between the two selected
        subsets and create a new subset_group with the genes that are
        differentially expressed. 
        
        This assumes that these subsets have been defined properly on the 
        obs array
        
        Note that this copies the relevant subsets into memory as otherwise
        this won't work on anndata in disk-backed mode:
        
        https://github.com/theislab/scanpy/issues/2147
        
        (the above is technically for a different scanpy function, but the same problem occurs for rank_genes_groups)
        """
        gene_list, dge_data = get_gene_list_diff_exp(self.state.subset1, self.state.subset2, self.state.data, n_genes=self.state.num_genes)
        new_name = f'DGE between {self.state.subset1.label} and {self.state.subset2.label}'
        dge_data.label = new_name
        self.state.data_collection.append(dge_data)

        vardata = self.state.data.meta['var_data']

        genelink = JoinLink(cids1 = [vardata.id['var_names']], 
                          cids2 = [dge_data.id['gene_names']], data1 = vardata, data2 = dge_data)
        self.state.data_collection.add_link(genelink)

        all_indices = []
        for gene in gene_list: # There is probably a more efficient way to get the codes for certain categories
            matching_indices = np.where(vardata[self.state.gene_att] == gene)
            all_indices.extend(list(matching_indices[0]))
        gene_codes = vardata[self.state.gene_att][all_indices].codes
        gene_state = CategorySubsetState(att=vardata.id[self.state.gene_att],
                                         categories=gene_codes)
        new_name = f'DGE between {self.state.subset1.label} and {self.state.subset2.label}'

        self.state.data_collection.new_subset_group(new_name, gene_state)

        confirm = dialog('New subset created',
                f'The subset:\n'
                f'{new_name}\n'
                f'has been created.',
                'info')


    @classmethod
    def create_subset(cls, collect, default=None, parent=None):
        self = cls(collect, parent=parent, default=default)
        value = self.exec_()

        if value == QtWidgets.QDialog.Accepted:
            self._apply()
