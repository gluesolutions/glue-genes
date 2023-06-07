from glue_genes.glue_single_cell.component import SyncComponent
from glue.core.data import Data
from glue.core.component_id import ComponentID
import numpy as np
from glue.config import colormaps


def test_synccomponent():

    n_cmaps_orig = len(colormaps)
    data = Data(label="Test Data")
    sub = data.new_subset()
    assert sub.style.color == '#e31a1c'
    cid = ComponentID("SyncComponent (sync)")
    comp = SyncComponent(np.random.random((2, 3)), subsets=[sub])
    cid2 = data.add_component(comp, cid)
    assert cid2 == cid
    n_cmaps_new = len(colormaps)
    assert n_cmaps_new == n_cmaps_orig+1
    assert data._components[cid] == comp

    assert comp.preferred_cmap == colormaps.members[-1][1]
