from glue_qt.app import GlueApplication
from glue.core import data_factories as df
from glue.core.data_region import RegionData
from glue.core.component import ExtendedComponent
from glue.core.state import GlueUnSerializer
from glue_genes.glue_genomics_data.spaceranger_factory import read_spaceranger_directory
import os
import pytest

DATA = os.path.join(os.path.dirname(__file__), "test_spatial_data")


class TestSpacerangerData(object):
    def get_data(self):
        test_data = df.load_data(
            DATA,
            factory=read_spaceranger_directory,
        )
        return test_data

    def test_spaceranger_access(self):
        self.app = GlueApplication()
        self.session = self.app.session
        self.hub = self.session.hub

        self.data_collection = self.session.data_collection
        test_data = self.get_data()
        self.data_collection.append(test_data)

        assert len(self.data_collection) == 4

    # The filter is to remove a harmless warning about np.save and meta on dtypes
    # that is triggered by the test dataset but not real-world datasets
    @pytest.mark.filterwarnings("ignore::UserWarning")
    def test_save_and_restore(self, tmpdir):
        self.app = GlueApplication()
        self.session = self.app.session
        self.hub = self.session.hub

        self.data_collection = self.session.data_collection
        test_data = self.get_data()
        self.data_collection.append(test_data)

        assert len(self.data_collection) == 4

        # Save the session

        filename = tmpdir.join("test_spaceranger_save_restore_session.glu").strpath
        self.session.application.save_session(filename)
        self.app.close()

        with open(filename, "r") as f:
            session = f.read()

        state = GlueUnSerializer.loads(session)
        ga = state.object("__main__")
        dc = ga.session.data_collection

        assert len(dc) == 4
        assert dc[0].label == "custom_X"
        assert dc[1].label == "custom_vars"
        assert dc[2].label == "custom_obs"
        assert dc[3].label == "custom_image"

        assert isinstance(dc[2], RegionData)
        assert isinstance(dc[2].get_component(dc[2].components[7]), ExtendedComponent)
        # Make sure component is restored correctly
        assert dc[2].center_x_id == dc[2].components[5]
        assert dc[2].center_y_id == dc[2].components[6]
