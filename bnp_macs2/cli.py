"""Console script for bnp_macs2."""
import typer
from typing import List
import dataclasses
from scipy.special import pdtrc
import numpy as np
import logging

from bionumpy.datatypes import Interval, Bed6
from bionumpy.arithmetics.geometry import Geometry, GenomicTrack
#from .control_pileup import get_control_pileup
#from .fragment_pileup import get_fragment_pileup
#from .call_peaks import call_peaks
# import matplotlib.pyplot as plt
import bionumpy as bnp

logging.basicConfig(level=logging.INFO)


@dataclasses.dataclass
class Macs2Params:
    fragment_length: int = 150
    n_reads: int = None
    p_value_cutoff: float = 0.001
    max_gap: int = 30
    write_bdg: bool = False


def logsf(count: float, mu: float) -> float:
    return np.log(pdtrc(count, mu))



def macs2(intervals: Interval, geometry: Geometry, params: Macs2Params):
    if not len(intervals):
        return Interval.empty()
    p_values = get_p_values(intervals, geometry,
                            params.fragment_length,
                            params.n_reads/geometry.size())
    return call_peaks(p_values, params.p_value_cutoff, params.fragment_length, params.max_gap, geometry)


class Macs2:
    def __init__(self, geometry: Geometry, params: Macs2Params):
        self._geometry = geometry
        self._params = params

    def run(self, intervals: Interval):
        p_values = self.get_p_values(intervals)
        return self.call_peaks(p_values)

    def get_fragment_pileup(self, reads: Bed6) -> GenomicTrack:
        logging.info("Getting fragment pileup")
        fragments = self._geometry.extend_to_size(reads, self._params.fragment_length)
        return self._geometry.get_pileup(fragments)

    def _get_average_pileup(self, reads: Bed6, window_size: int) -> GenomicTrack:
        logging.info(f"Getting average pileup for window_size {window_size}")
        windows = dataclasses.replace(reads,
                                      start=reads.start-window_size//2,
                                      stop=reads.start+window_size//2)
        clipped = self._geometry.clip(windows)
        return self._geometry.get_pileup(clipped)/window_size

    def get_control_pileup(self, reads: Bed6, window_sizes: List[int]) -> GenomicTrack:
        logging.info('Getting control lambda')
        pileup = float(self._params.n_reads/self._geometry.size())
        for window_size in window_sizes:
            avg_pileup = self._get_average_pileup(reads, window_size)
            pileup = np.maximum(pileup, avg_pileup)
        return pileup

    def get_p_values(self, intervals: Interval):
        fragment_pileup = self.get_fragment_pileup(intervals)
        control = self._params.fragment_length*self.get_control_pileup(intervals, [1000, 10000])
        p_values = logsf(fragment_pileup, control)
        return p_values

    def call_peaks(self, log_p_values):
        peaks = log_p_values < np.log(self._params.p_value_cutoff)
        peaks = peaks.to_intervals()
        peaks = self._geometry.merge_intervals(peaks, distance=self._params.max_gap)
        peaks = peaks[(peaks.stop-peaks.start) >= self._params.fragment_length]
        return peaks


def main(filename: str, genome_file: str, fragment_length: int = 150, p_value_cutoff: float = 0.001, outfilename: str = None):
    genome = bnp.open(genome_file, buffer_type=bnp.io.files.ChromosomeSizeBuffer).read()
    chrom_sizes = {str(name): size for name, size in zip(genome.name, genome.size)}
    intervals = bnp.open(filename, buffer_type=bnp.io.delimited_buffers.Bed6Buffer).read()
    tag_size = np.median(intervals.stop-intervals.start)
    params = Macs2Params(
        fragment_length=fragment_length,
        p_value_cutoff=p_value_cutoff,
        max_gap=int(tag_size),
        n_reads=bnp.count_entries(filename))
    m = Macs2(Geometry(chrom_sizes), params)
    result = m.run(intervals)
    if outfilename is not None:
        with bnp.open(outfilename, 'w') as f:
            f.write(result)
    return result


def test():
    #res = main("example_data/small_interval.bed", "example_data/small_genome.fa.fai")
    res = main('example_data/simulated_chip_seq.bed', 'example_data/simulated.chrom.sizes', fragment_length=200)
    for chunk in res:
        print(chunk)


def big():
    with bnp.open("/home/knut/Data/peaks.bed", "w") as outfile:
        for data in main("/home/knut/Data/subset.bed", "/home/knut/Data/hg38.chrom.sizes"):
            if len(data):
                outfile.write(data)


if __name__ == "__main__":

    typer.run(main)


def run():
    typer.run(main)


if __name__ == "__main__":
    typer.run(main)
