"""Console script for bnp_macs2."""
from typing import List
import dataclasses
from scipy.special import pdtrc
import numpy as np
import logging

from bionumpy.datatypes import Interval, Bed6, NarrowPeak
from bionumpy.arithmetics.geometry import Geometry, GenomicTrack

from .listener import Listner

logger = logging.getLogger(__name__)


@dataclasses.dataclass
class Macs2Params:
    fragment_length: int = 150
    n_reads: int = None
    p_value_cutoff: float = 0.001
    max_gap: int = 30
    write_bdg: bool = False


def logsf(count: float, mu: float) -> float:
    return np.log(pdtrc(count, mu))


class register:
    def __init__(self, name):
        self._name = name

    def __call__(self, func):
        def new_func(obj, *args, **kwargs):
            result = func(obj, *args, **kwargs)
            if hasattr(obj, '_listner') and obj._listner is not None:
                getattr(obj._listner, self._name)(result)
            return result
        return new_func


class Macs2:
    def __init__(self, geometry: Geometry, params: Macs2Params, listner: Listner=None):
        self._geometry = geometry
        self._params = params
        self._listner = listner

    def dedup(self, intervals):
        s = self._geometry.sort(intervals)
        filtered = s[np.insert(s.start[1:] != s.start[:-1], 0, True)]
        logger.info(f'{len(filtered)} reads remaining of {len(s)}')
        return filtered

    @register('peaks')
    def run(self, intervals: Interval) -> Interval:
        fragment_pileup = self.get_fragment_pileup(intervals)
        print(fragment_pileup.sum())
        control = self.get_control_pileup(intervals, [10000])
        p_scores = logsf(fragment_pileup, control)
        if self._listner is not None:
            self._listner.p_scores(p_scores)
        print(p_scores)
        peaks = self.call_peaks(p_scores)
        return self.get_narrow_peak(peaks, np.log10(np.e)*-p_scores)

    @register('treat_pileup')
    def get_fragment_pileup(self, reads: Bed6) -> GenomicTrack:
        fragments = self._geometry.extend_to_size(reads, self._params.fragment_length)
        return self._geometry.get_pileup(fragments)

    def _get_average_pileup(self, reads: Bed6, window_size: int) -> GenomicTrack:
        windows = dataclasses.replace(reads,
                                      start=reads.start-window_size//2,
                                      stop=reads.start+window_size//2)
        clipped = self._geometry.clip(windows)
        return self._geometry.get_pileup(clipped)/window_size

    @register('control_lambda')
    def get_control_pileup(self, reads: Bed6, window_sizes: List[int]) -> GenomicTrack:
        pileup = float(self._params.n_reads/self._geometry.size())
        print(self._params.n_reads, self._geometry.size())
        for window_size in window_sizes:
            avg_pileup = self._get_average_pileup(reads, window_size)
            print(avg_pileup)
            pileup = np.maximum(pileup, avg_pileup)
        print(pileup*self._params.fragment_length)
        return pileup*self._params.fragment_length

    # @register('peaks')
    def call_peaks(self, log_p_values):
        peaks = log_p_values < np.log(self._params.p_value_cutoff)
        peaks = peaks.to_intervals()
        assert len(peaks) > 0
        peaks = self._geometry.merge_intervals(peaks, distance=self._params.max_gap)
        peaks = peaks[(peaks.stop-peaks.start) >= self._params.fragment_length]
        return peaks

    def get_narrow_peak(self, peaks: Interval, p_values: GenomicTrack):
        peak_signals = p_values.get_intervals(peaks, stranded=False)
        max_values = peak_signals.max(axis=-1)
        mean_values = peak_signals.mean(axis=-1)
        print(mean_values)
        return NarrowPeak(
            peaks.chromosome,
            peaks.start,
            peaks.stop,
            [f'peak_{i+1}' for i in range(len(peaks))],
            (max_values*10).astype(int),
            ['.']*len(peaks),
            mean_values,
            max_values,
            max_values,
            [0]*len(peaks))

