import logging
from bionumpy.datatypes import Bed6
from bionumpy.arithmetics.geometry import Geometry
import dataclasses
import numpy as np
logger = logging.getLogger(__name__)


def get_average_pileup(reads: Bed6, window_size: int, geometry: Geometry):
    logging.info(f"Getting average pileup for window_size {window_size}")
    windows = dataclasses.replace(reads, 
                                  start=reads.start-window_size//2,
                                  stop=reads.start+window_size//2)
    clipped = geometry.clip(windows)
    return geometry.get_pileup(clipped)/window_size


def get_control_pileup(reads, window_sizes, read_rate, geometry):
    pileup = float(read_rate)
    for window_size in window_sizes:
        avg_pileup = get_average_pileup(reads, window_size, geometry)
        pileup = np.maximum(pileup, avg_pileup)
    return pileup
