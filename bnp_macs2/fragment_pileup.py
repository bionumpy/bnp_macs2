import logging
from bionumpy.datatypes import Bed6
from bionumpy.arithmetics.geometry import Geometry
logger = logging.getLogger(__name__)


def get_fragment_pileup(reads: Bed6, fragment_length: int, geometry: Geometry):
    logging.info("Getting fragment pileup")
    fragments = geometry.extend_to_size(reads, fragment_length)
    return geometry.get_pileup(fragments)
