import glue_qt.viewers.scatter.data_viewer as old
from .data_viewer import ScatterViewer  # noqa


def setup():
    from glue_qt.config import qt_client
    qt_client._members = [x for x in qt_client._members if x != old.ScatterViewer]
    qt_client.add(ScatterViewer)
