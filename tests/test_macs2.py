import numpy as np
from bnp_macs2.fragment_pileup import get_fragment_pileup
from bnp_macs2.control_pileup import get_average_pileup, get_control_pileup
from bnp_macs2.call_peaks import call_peaks
from bnp_macs2.cli import Macs2Params, macs2
from bionumpy import Bed6, str_equal
from bionumpy.datatypes import Interval
from bionumpy.arithmetics.geometry import Geometry, GenomicTrack
from bionumpy.arithmetics.intervals import GenomicRunLengthArray
from bionumpy.arithmetics.global_offset import GlobalOffset
from bionumpy.util.testing import assert_bnpdataclass_equal
import pytest


@pytest.fixture
def intervals():
    return Bed6.from_entry_tuples(
        [('chr1', 10, 20, '.', '.', '-'),
         ('chr1', 11, 22, '.', '.', '+'),
         ('chr1', 40, 60, '.', '.', '-'),
         ('chr1', 15, 35, '.', '.', '+'),
         ('chr2', 15, 35, '.', '.', '+')])

@pytest.fixture
def chrom_sizes():
    return {'chr1': 100, 'chr2': 60}


@pytest.fixture
def geometry(chrom_sizes):
    return Geometry(chrom_sizes)# {'chr1': 100, 'chr2': 60})

@pytest.fixture
def pileup(chrom_sizes):
    global_offset = GlobalOffset(chrom_sizes)
    dense_pileup = np.full(160, 0.06, dtype=float)
    dense_pileup[10:30] = 0.04
    dense_pileup[39:60] = 0.03
    dense_pileup[71:79] = 0.02
    dense_pileup[110:119] = 0.01
    rle = GenomicRunLengthArray.from_array(dense_pileup)
    return GenomicTrack.from_global_data(rle, global_offset)

@pytest.fixture
def peaks():
    return Interval.from_entry_tuples([
        ('chr1', 10, 60)])

def dense_fragment_pileup(intervals, fragment_length, size):
    pileup = np.zeros(size, dtype=int)
    intervals.strand = intervals.strand.ravel()
    for interval in intervals:
        if interval.strand == '+':
            pileup[interval.start:interval.start+fragment_length] += 1
        else:
            pileup[interval.stop-fragment_length:interval.stop] += 1

    return pileup


def dense_average_pileup(intervals, window_size, size):
    pileup = np.zeros(size, dtype=int)
    for i in range(size):
        for interval in intervals:
            if (i-window_size//2) < interval.start <= (i+window_size//2):
                pileup[i] += 1
    return pileup/window_size


def dense_control_pileup(intervals, window_sizes, read_rate, size):
    r = np.maximum(*[dense_average_pileup(intervals, window_size, size) for window_size in window_sizes])
    return np.maximum(r, read_rate)


def test_get_fragment_pileup(intervals: Bed6, geometry: Geometry):
    true_pileup = {name: dense_fragment_pileup(intervals[str_equal(intervals.chromosome, name)], 20, geometry.chrom_size(name)) for name in geometry.names()}
    fragment_pileup = get_fragment_pileup(intervals, 20, geometry)
    np.testing.assert_equal(true_pileup, fragment_pileup.to_dict())


@pytest.mark.parametrize('window_size', [10, 20])
def test_get_average_pileup(intervals, geometry, window_size):
    true_pileup = {name: dense_average_pileup(intervals[str_equal(intervals.chromosome, name)], window_size, geometry.chrom_size(name)) for name in geometry.names()}
    avg_pileup = get_average_pileup(intervals, window_size, geometry)
    np.testing.assert_equal(true_pileup, avg_pileup.to_dict())


def test_get_control_pileup(intervals, geometry):
    true_pileup = {name: dense_control_pileup(intervals[str_equal(intervals.chromosome, name)], [10, 20], 0.1, geometry.chrom_size(name)) for name in geometry.names()}
    control_pileup = get_control_pileup(intervals, [10, 20], 0.1, geometry)
    np.testing.assert_equal(true_pileup, control_pileup.to_dict())    


def test_call_peaks(pileup, geometry, peaks):
    called_peaks = call_peaks(np.log(pileup), 0.05, 20, 10, geometry)
    called_peaks.chromosome = called_peaks.chromosome.encoding.decode(called_peaks.chromosome)
    assert_bnpdataclass_equal(called_peaks, peaks)


def testmacs2_acceptance(intervals, geometry):
    params = Macs2Params(fragment_length=20,
                         max_gap=10,
                         p_value_cutoff=0.05,
                         n_reads = len(intervals))
    macs2(intervals, geometry, params)
