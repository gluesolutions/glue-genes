from glue_genes.glue_single_cell.component import SyncComponent
from glue.core.data import Data
from glue.core.component_id import ComponentID
import numpy as np
from glue.config import colormaps
import pytest
from matplotlib.colors import ColorConverter


def test_synccomponent():
    n_cmaps_orig = len(colormaps)
    data = Data(label="Test Data")
    sub = data.new_subset()
    sub.style.color == "#e31a1c"
    cid = ComponentID("SyncComponent (sync)")
    comp = SyncComponent(np.random.random((2, 3)), subsets=[sub])
    cid2 = data.add_component(comp, cid)
    assert cid2 == cid
    n_cmaps_new = len(colormaps)
    assert n_cmaps_new == n_cmaps_orig + 1
    assert data._components[cid] == comp
    assert comp.preferred_cmap == colormaps.members[-1][1]


def test_synccomponent_no_subsets():
    data = Data(label="Test Data")
    cid = ComponentID("SyncComponent (sync)")
    comp = SyncComponent(np.random.random((2, 3)))
    _ = data.add_component(comp, cid)
    assert comp.preferred_cmap is None


def test_synccomponent_two_subsets():
    data = Data(label="Test Data")
    cid = ComponentID("SyncComponent (sync)")
    sub1 = data.new_subset()
    sub2 = data.new_subset()
    sub1.style.color == "#e31a1c"
    sub2.style.color == "#1f78b4"
    comp = SyncComponent(np.random.random((2, 3)), subsets=[sub1, sub2])
    _ = data.add_component(comp, cid)
    assert comp.preferred_cmap == colormaps.members[-1][1]
    assert comp.cmap_name == f"{sub1.label}_{sub2.label}_cmap"


def test_synccomponent_three_subsets():
    data = Data(label="Test Data")
    sub1 = data.new_subset()
    sub2 = data.new_subset()
    sub3 = data.new_subset()
    sub1.style.color == "#e31a1c"
    sub2.style.color == "#1f78b4"
    sub3.style.color == "#33a02c"
    with pytest.raises(NotImplementedError):
        _ = SyncComponent(np.random.random((2, 3)), subsets=[sub1, sub2, sub3])


def test_synccomponent_cmap():
    data = Data(label="Test Data")
    cid = ComponentID("SyncComponent (sync)")
    sub1 = data.new_subset()
    sub2 = data.new_subset()
    sub1.style.color == "#e31a1c"
    sub2.style.color == "#1f78b4"
    comp = SyncComponent(np.random.random((2, 3)), subsets=[sub1, sub2])
    _ = data.add_component(comp, cid)
    assert comp.preferred_cmap == colormaps.members[-1][1]
    assert comp.cmap_name == f"{sub1.label}_{sub2.label}_cmap"
    assert comp.preferred_cmap(0) == ColorConverter().to_rgba(sub1.style.color)
    assert comp.preferred_cmap(255) == ColorConverter().to_rgba(sub2.style.color)
