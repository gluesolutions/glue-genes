from glue.app.qt import GlueApplication
from glue.core import Data
from glue.core.roi import RectangularROI
from glue.core.state import GlueUnSerializer
from glue.utils.qt import process_events
from numpy.testing import assert_equal

from ..viewer import QTLViewer


class TestScatterViewer(object):
    def setup_method(self, method):

        self.data = Data(
            label="d1",
            x=[1125339623, 3025778038, 3025778038, 497597140],
            y=[1125562888, 3026087386, 259628224, 498042878],
            z=[41, 28, 11, 9],
            z2=[0.1, 0.3, 0.4, 0.5],
        )

        self.app = GlueApplication()
        self.session = self.app.session
        self.hub = self.session.hub

        self.data_collection = self.session.data_collection
        self.data_collection.append(self.data)

        self.viewer = self.app.new_data_viewer(QTLViewer)

    # def teardown_method(self, method):
    #    self.viewer.close()
    #    self.viewer = None
    #    self.app.close()
    #    self.app = None

    def test_basic(self):

        viewer_state = self.viewer.state

        # Check defaults when we add data
        self.viewer.add_data(self.data)

        assert viewer_state.x_att is self.data.id["x"]
        assert viewer_state.y_att is self.data.id["y"]
        assert viewer_state.lod_att is self.data.id["z"]

        assert len(viewer_state.layers) == 1

    def test_lod_limits(self):
        viewer_state = self.viewer.state
        self.viewer.add_data(self.data)
        assert viewer_state.lod_min == min(self.data["z"])
        assert viewer_state.lod_max == max(self.data["z"])

        assert viewer_state.lod_thresh == viewer_state.lod_min

        viewer_state.lod_att = self.data.id["z2"]
        assert viewer_state.lod_thresh == viewer_state.lod_min

    def test_lod_roi(self):
        """
        Currently lod_thresh creates a lod_mask inside
        the artist only. It should probably happen
        in a state variable that the artist accesses.

        This would, probably, eliminate the need for calling
        the apply_roi thing again?

        NO! Changing lod_thresh does NOT change the subset
        selection unless the user is explicit
        """

        viewer_state = self.viewer.state
        self.viewer.add_data(self.data)

        # Set lod_threshold to include just two points
        viewer_state.lod_thresh = 20
        roi = RectangularROI(0, 4000000000, 0, 4000000000)
        self.viewer.apply_roi(roi)
        assert len(self.viewer.layers) == 2
        sub_data = self.viewer.layers[1].plot_artist.get_data()
        assert len(sub_data[0]) == 2

        # Set lod_threshold to include just one point
        viewer_state.lod_thresh = 40
        roi = RectangularROI(0, 4000000000, 0, 4000000000)
        self.viewer.apply_roi(roi)
        assert len(self.viewer.layers) == 2
        sub_data = self.viewer.layers[1].plot_artist.get_data()
        assert len(sub_data[0]) == 1  # This is a masked array

    def test_lod_roi_with_cmap_mode(self):
        # Regression test to make sure LOD works with cmap_mode
        # The data is set correctly, but then _update_visual_attributes
        # does not run (correctly?) because we set c and s as well, which
        # means things do not update properly.

        # This seems like a code of code duplication to get a reduced view
        # I can see a potential here for a limitedscatterviewer artist/state
        # that I can use in both QTL viewer and small multiples
        viewer_state = self.viewer.state
        self.viewer.add_data(self.data)
        self.viewer.layers[0].state.cmap_mode = "Linear"
        self.viewer.layers[0].state.cmap_att = self.data.id["z"]

        # Set lod_threshold to include just two points
        viewer_state.lod_thresh = 20
        roi = RectangularROI(0, 4000000000, 0, 4000000000)
        self.viewer.apply_roi(roi)
        self.viewer.layers[1].state.cmap_mode = "Linear"
        self.viewer.layers[1].state.cmap_att = self.data.id["z"]

        assert len(self.viewer.layers) == 2
        sub_data = self.viewer.layers[1].scatter_artist.get_offsets()
        assert len(sub_data[0]) == 2

        # Set lod_threshold to include just one point
        viewer_state.lod_thresh = 40
        roi = RectangularROI(0, 4000000000, 0, 4000000000)
        self.viewer.apply_roi(roi)
        self.viewer.layers[1].state.cmap_mode = "Linear"
        self.viewer.layers[1].state.cmap_att = self.data.id["z"]

        assert len(self.viewer.layers) == 2
        sub_data = self.viewer.layers[1].scatter_artist.get_offsets()

        assert_equal(
            sub_data, [[1125339623, 1125562888]]
        )  # This is not a masked array?

    def test_session_save_and_restore(self, tmpdir):
        viewer_state = self.viewer.state
        self.viewer.add_data(self.data)
        viewer_state.lod_att = self.data.id["z2"]
        viewer_state.lod_thresh = 0.3
        process_events()
        filename = tmpdir.join("test_qtl_session.glu").strpath

        self.session.application.save_session(filename)

        with open(filename, "r") as f:
            session = f.read()

        state = GlueUnSerializer.loads(session)

        ga = state.object("__main__")

        dc = ga.session.data_collection

        viewer = ga.viewers[0][0]
        assert viewer.state.lod_att is dc[0].id["z2"]
        assert viewer.state.lod_thresh == 0.3
        ga.close()
