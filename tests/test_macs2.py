import numpy as np
from numpy.testing import assert_equal
from bnp_macs2.cli import Macs2Params, Macs2
from bnp_macs2.listener import DebugListnerStream
from bionumpy import Bed6, str_equal
from bionumpy.datatypes import Interval
from bionumpy.genomic_data import GenomicTrack, GenomicIntervals
from bionumpy.genomic_data.geometry import Geometry, StreamedGeometry
from bionumpy.arithmetics.intervals import GenomicRunLengthArray
from bionumpy.genomic_data.global_offset import GlobalOffset
from bionumpy.util.testing import assert_bnpdataclass_equal
from bionumpy.streams import NpDataclassStream
import pytest


@pytest.fixture
def debug_listner():
    return DebugListnerStream()


@pytest.fixture
def intervals():
    return Bed6.from_entry_tuples(
        [('chr1', 10, 20, '.', '.', '-'),
         ('chr1', 11, 22, '.', '.', '+'),
         ('chr1', 40, 60, '.', '.', '-'),
         ('chr1', 15, 35, '.', '.', '+'),
         ('chr2', 15, 35, '.', '.', '+')])


@pytest.fixture
def params(chrom_sizes, intervals):
    return Macs2Params(fragment_length=20,
                       max_gap=10,
                       p_value_cutoff=0.05,
                       n_reads=len(intervals),
                       effective_genome_size=sum(chrom_sizes.values()))


@pytest.fixture
def macs2_obj(params, debug_listner):
    return Macs2(params, debug_listner)


@pytest.fixture
def chrom_sizes():
    return {'chr1': 100, 'chr2': 60}


@pytest.fixture
def genomic_intervals(intervals, chrom_sizes):
    return GenomicIntervals.from_intervals(intervals, chrom_sizes)


@pytest.fixture
def geometry(chrom_sizes):
    return Geometry(chrom_sizes)# {'chr1': 100, 'chr2': 60})


@pytest.fixture
def streamed_geometry(chrom_sizes):
    return StreamedGeometry(chrom_sizes)# {'chr1': 100, 'chr2': 60})


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


def test_get_fragment_pileup(intervals: Bed6, macs2_obj, geometry, genomic_intervals):
    true_pileup = {name: dense_fragment_pileup(intervals[str_equal(intervals.chromosome, name)], 20, geometry.chrom_size(name)) for name in geometry.names()}

    fragment_pileup = macs2_obj.get_fragment_pileup(genomic_intervals) # , geometry)
    np.testing.assert_equal(true_pileup, fragment_pileup.to_dict())


@pytest.mark.parametrize('window_size', [10, 20])
def test_get_average_pileup(intervals, geometry, window_size, macs2_obj, genomic_intervals):
    true_pileup = {name: dense_average_pileup(intervals[str_equal(intervals.chromosome, name)], window_size, geometry.chrom_size(name)) for name in geometry.names()}
    avg_pileup = macs2_obj._get_average_pileup(genomic_intervals, window_size)# , geometry)
    np.testing.assert_equal(true_pileup, avg_pileup.to_dict())


def test_get_control_pileup(intervals, geometry, macs2_obj, genomic_intervals):
    read_rate = len(intervals)/macs2_obj.params.effective_genome_size
    true_pileup = {name: dense_control_pileup(intervals[str_equal(intervals.chromosome, name)], [10, 20], read_rate, geometry.chrom_size(name))*macs2_obj.params.fragment_length for name in geometry.names()}
    control_pileup = macs2_obj.get_control_pileup(genomic_intervals, [10, 20])# , 0.1, geometry)
    np.testing.assert_equal(true_pileup, control_pileup.to_dict())


def test_call_peaks(pileup, peaks, macs2_obj):
    called_peaks = macs2_obj.call_peaks(np.log(pileup)).get_data()
    print(called_peaks.chromosome.encoding)
    # called_peaks.chromosome = called_peaks.chromosome.encoding.decode(called_peaks.chromosome)
    assert_bnpdataclass_equal(called_peaks, peaks)


def testmacs2_acceptance(intervals, chrom_sizes, macs2_obj):
    genomic_intervals = GenomicIntervals.from_intervals(intervals, chrom_sizes)
    macs2_obj.run(genomic_intervals)


def testmacs2_acceptance_stream(intervals, chrom_sizes, macs2_obj: Macs2):
    stream = NpDataclassStream(iter([intervals]))
    genomic_intervals = GenomicIntervals.from_interval_stream(stream, chrom_sizes)
    peaks = macs2_obj.run(genomic_intervals)
    genomic_intervals = GenomicIntervals.from_intervals(intervals, chrom_sizes)
    real_peaks = macs2_obj.run(genomic_intervals)
    assert_equal(peaks.start, real_peaks.start)
