from glue.core.data_combo_helper import ComboHelper, ComponentIDComboHelper
from glue.core.state_objects import StateAttributeLimitsHelper
from glue.viewers.matplotlib.state import DeferredDrawCallbackProperty as DDCProperty
from glue.viewers.matplotlib.state import (
    DeferredDrawSelectionCallbackProperty as DDSCProperty,
)
from glue.viewers.scatter.state import ScatterViewerState

__all__ = ["QTLViewerState"]

CHR_POSITIONS = {
    # http://www.informatics.jax.org/mgihome/other/mouse_facts1.shtml
    "Mouse": {
        "Names": [
            "1",
            "2",
            "3",
            "4",
            "5",
            "6",
            "7",
            "8",
            "9",
            "10",
            "11",
            "12",
            "13",
            "14",
            "15",
            "16",
            "17",
            "18",
            "19",
            "X",
            "Y",
            "MT",
        ],
        "Lengths": [
            195,
            182,
            160,
            157,
            152,
            150,
            145,
            130,
            124,
            131,
            122,
            120,
            121,
            125,
            104,
            98,
            95,
            91,
            61,
            169,
            91,
            0.01,
        ],
        "GridSize": 200_000_000,  # Genome position seems to use a fixed grid
        # so that there is emtpy space in later chromosomes
    },
    "Human": {  # http://www.insilicase.com/Web/Chromlen.aspx
        "Names": [
            "1",
            "2",
            "3",
            "4",
            "5",
            "6",
            "7",
            "8",
            "9",
            "10",
            "11",
            "12",
            "13",
            "14",
            "15",
            "16",
            "17",
            "18",
            "19",
            "20",
            "21",
            "22",
            "22",
            "X",
            "Y",
            "MT",
        ],
        "GridSize": 250_000_000,  # This is a guess
    },
}
UNITS_LOOKUP = {1: "bp", 1000: "kb", 1_000_000: "Mb", 1_000_000_000: "Gb"}


class QTLViewerState(ScatterViewerState):
    """
    A state class that includes all the attributes for a QTL Viewer

    This is a fairly minimal small set of extra attributes on top
    of the base :class:`~glue.viewers.scatter.state.ScatterViewerState`.
    """

    species = DDSCProperty(
        0, docstring="The species (for showing chromosome boundaries)"
    )
    pos_units = DDSCProperty(0, docstring="Units for gene and marker position")

    lod_att = DDSCProperty(
        docstring="The attribute giving the LOD score ", default_index=2
    )
    lod_thresh = DDCProperty(-99, docstring="The LOD threshold for display and subsets")

    lod_min = DDCProperty(docstring="Lower limit of the visible x range")
    lod_max = DDCProperty(docstring="Upper limit of the visible x range")

    def __init__(self, **kwargs):
        super().__init__()

        self.lod_lim_helper = StateAttributeLimitsHelper(
            self,
            attribute="lod_att",
            lower="lod_min",
            upper="lod_max",
            margin=0,
        )

        self.lod_att_helper = ComponentIDComboHelper(
            self, "lod_att", numeric=True, categorical=False
        )

        self.species_helper = ComboHelper(self, "species")
        self.pos_units_helper = ComboHelper(self, "pos_units")

        self.chr_pos = CHR_POSITIONS

        self.species_helper.choices = list(self.chr_pos.keys())
        try:
            self.species_helper.selection = self.species_helper.choices[0]
        except IndexError:
            pass
        self.species_helper.display = str

        def display_unit_names(unit):
            return UNITS_LOOKUP[unit]

        self.pos_units_helper.choices = [1, 1000, 1_000_000, 1_000_000_000]
        try:
            self.pos_units_helper.selection = self.pos_units_helper.choices[0]
        except IndexError:
            pass
        self.pos_units_helper.display = display_unit_names

        # self._adjust_lod_thresh()
        self.update_from_dict(kwargs)
        # Setup callback only after update_from_dict
        # otherwise session files do not save lod_thresh
        self.add_callback("lod_att", self._adjust_lod_thresh)

    def _layers_changed(self, *args):

        layers_data = self.layers_data
        layers_data_cache = getattr(self, "_layers_data_cache", [])

        if layers_data == layers_data_cache:
            return

        self.x_att_helper.set_multiple_data(self.layers_data)
        self.y_att_helper.set_multiple_data(self.layers_data)
        self.lod_att_helper.set_multiple_data(self.layers_data)

        self._layers_data_cache = layers_data

    def _adjust_lod_thresh(self, *args):
        """
        When we change lod_add, we'll set the threhold to the
        minimum value in order to show all the data.
        """
        if self.lod_att is None:
            return
        self.lod_thresh = self.lod_min
