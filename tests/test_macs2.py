import numpy as np
from bnp_macs2.fragment_pileup import get_fragment_pileup
from bionumpy import Bed6
from bionumpy.arithmetics.geometry import Geometry
import pytest


@pytest.fixture
def intervals():
    return Bed6.from_entry_tuples(
        [('chr1', 10, 20, '.', '.', '-'),
         ('chr1', 11, 22, '.', '.', '+'),
         ('chr1', 40, 60, '.', '.', '-'),
         ('chr1', 15, 35, '.', '.', '+')])


@pytest.fixture
def geometry():
    return Geometry({'chr1': 100})


def dense_fragment_pileup(intervals, fragment_length, size):
    pileup = np.zeros(size, dtype=int)
    intervals.strand = intervals.strand.ravel()
    for interval in intervals:
        if interval.strand == '+':
            pileup[interval.start:interval.start+fragment_length] += 1
        else:
            pileup[interval.stop-fragment_length:interval.stop] += 1

    return pileup


def test_get_fragment_pileup(intervals: Bed6, geometry: Geometry):
    true_pileup = {name: dense_fragment_pileup(intervals, 20, geometry.chrom_size(name)) for name in geometry.names()}
    fragment_pileup = get_fragment_pileup(intervals, 20, geometry)
    np.testing.assert_equal(true_pileup, fragment_pileup.to_dict())
