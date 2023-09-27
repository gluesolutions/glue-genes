from glue.config import menubar_plugin

from .qt.diff_gene_exp import DiffGeneExpDialog
from .qt.enrichr import EnrichpyDialog
from .qt.summarize_gene_subset import SummarizeGeneSubsetDialog


@menubar_plugin("Calculate Differential Gene Expression")
def diff_gene_exp_plugin(session, data_collection):
    DiffGeneExpDialog.calculate_deg(data_collection, default=None, parent=None)
    return


@menubar_plugin("Calculate Summary Over Gene Subset")
def summarize_gene_subset_exp_plugin(session, data_collection):
    SummarizeGeneSubsetDialog.summarize(data_collection, default=None, parent=None)
    return


@menubar_plugin("Enrich Gene Set via Enrichr")
def enrichrpy_plugin(session, data_collection):
    EnrichpyDialog.enrich(data_collection, default=None, parent=None)
    return
