from glue.core.link_helpers import LinkTwoWay, MultiLink
from glue.core.link_helpers import link_helper
from glue.core.link_helpers import BaseMultiLink

class BaseLinearScaleLink(BaseMultiLink):
    """
    This class can be used to simply add a scale and shift
    transform between (e.g.) images and annotations as:

        @link_helper(category='Images')
        class Fullres_to_Hires(BaseLinearScaleLink):
            description = 'Link Fullres data to HiRes data'
            labels1 = ['FullRes']
            labels2 = ['HiRes']
            display = "Fullres <-> HiRes"
            scale = 0.053972367
            shift = 0

    """
    def forwards(self, input):
        return self.scale * input + self.shift

    def backwards(self, input):
        return (input - self.shift) / self.scale


@link_helper(category="Images")
class FlipLink(MultiLink):

    display = "Flip y axis"
    description = ("Join two components while flipping the y axis. This"
                   "is sometimes necessary when image coordinates are"
                   "specified with (0,0) in the upper left corner.")

    labels1 = ["Image x coords", "Image y coords"]
    labels2 = ["Region/annotation x coords", "Region/annotation y coords"]

    def __init__(self, *args, cids1=None, cids2=None, data1=None, data2=None):
        self.data1 = data1
        image_height, _ = data1.shape
        def flip_y(coords):
            return image_height - coords

        self.data2 = data2
        image_x, image_y = cids1
        annot_x, annot_y = cids2
        x_to_x = LinkSame(image_x, annot_x)
        # AttributeError: 'LinkTwoWay' object has no attribute 'inverse'
        y_to_y = LinkTwoWay(image_y, annot_y, forwards=flip_y, backwards=flip_y)

        self._links = [x_to_x, y_to_y]
